import torch
from convex_hulls import ConvexHulls
import numpy as np
from numpy.linalg import inv

# define pi in torch
pi = torch.acos(torch.zeros(1)).item() * 2


class BilevelOptimizer(object):
    def __init__(self, hull0: ConvexHulls, target: ConvexHulls):
        self.beta = torch.Tensor([pi / 4]).requires_grad_(True)
        self.phi = torch.Tensor([pi / 4]).requires_grad_(True)
        self.theta = torch.Tensor([1, 1, 1]).requires_grad_(True)
        self.hull0 = hull0
        self.hull1 = target
        self.distance, self.closest_pos0, self.closest_pos1 = hull0.distance_between_convex_hulls(target)
        self.d = torch.Tensor([self.distance / 10]).requires_grad_(True)
        self.n = self._get_n()
        self.tmp_params = None

    # calculate the value of objective function
    def obj(self, mode="Normal"):
        # using the parameters from the optimizer
        if mode == "Normal":
            points0 = torch.Tensor(self.hull0.points)
            points1 = points0 + self.theta
            points2 = torch.Tensor(self.hull1.points)
            # vertices of two convex hulls
            v1 = points1[self.hull0.vertices, :]
            v2 = points2[self.hull1.vertices, :]
            # value of objective function
            objective = torch.sum(torch.log(torch.matmul(v2, self.n) - self.d)) + torch.sum(
                torch.log(self.d - torch.matmul(v1, self.n)))
            # objective *= -1
            return objective
        # using the tmp variables during line search
        if mode == "Temp":
            beta, phi, theta, d = self.tmp_params['beta'], \
                                  self.tmp_params['phi'], \
                                  self.tmp_params['theta'], \
                                  self.tmp_params['d']
            n = torch.cat([torch.cos(beta) * torch.cos(phi), torch.sin(beta) * torch.cos(phi), torch.sin(phi)]).reshape(
                (3, 1))
            points0 = torch.Tensor(self.hull0.points)
            #     ch1 = ConvexHull(points0 + theta)
            #     points1 = torch.Tensor(ch1.points)
            points1 = points0 + theta
            points2 = torch.Tensor(self.hull1.points)
            # vertices of two convex hulls
            #     v1 = points1[ch1.vertices, :].T
            v1 = points1[self.hull0.vertices, :]
            v2 = points2[self.hull1.vertices, :]
            # value of objective function
            objective = torch.sum(torch.log(torch.matmul(v2, n) - d)) + torch.sum(
                torch.log(d - torch.matmul(v1, n)))

            # objective *= -1
            return objective

    def line_search(self, niters: int = 500, tol: float = 1e-4, scale=0.7, c1=0.7, c2=0.8):
        result = self.obj()
        result.backward()
        for i in range(niters):
            old_obj = result
            print('last obj', old_obj)
            s = 1
            self._get_tmp_params(s)
            result_temp = self.obj(mode='Temp')
            print('finding feasible step size...')
            # keep objective value from Nan
            while torch.isnan(result_temp):
                s *= scale
                self._get_tmp_params(s)
                result_temp = self.obj(mode="Temp")
            print('s1 = ', s)
            print('finding proper step size to minimize the objective function...')

            # Armijo Condition
            while result_temp > old_obj - c1 * s * torch.sum(
                    torch.Tensor(
                        [self.beta.grad, self.phi.grad, self.d.grad, torch.sum(self.theta.grad)])) or torch.isnan(
                result_temp):
                tmp_s = s
                s *= scale
                if torch.abs(torch.Tensor([s - tmp_s])) <= tol:
                    break

            print('s2 = ', s)
            result_temp.backward()

            # Curvature condition

            self._update_params(s)
            self._reset_n()
            result = self.obj()
            result.backward()
            if torch.abs(result - old_obj) < tol:
                print("Converged")
                break
            print('Next Round')

    # reset beta and phi and n
    def _reset_n(self):
        rot_y_angle = -self.phi.detach()
        rot_z_angle = self.beta.detach()

        rot_y = torch.Tensor([[torch.cos(rot_y_angle), 0, torch.sin(rot_y_angle)],
                              [0, 1, 0],
                              [-torch.sin(rot_y_angle), 0, torch.cos(rot_y_angle)]]).requires_grad_(True)
        rot_z = torch.Tensor([[torch.cos(rot_z_angle), -torch.sin(rot_z_angle), 0],
                              [torch.sin(rot_z_angle), torch.cos(rot_z_angle), 0],
                              [0, 0, 1]]).requires_grad_(True)
        self.beta.data.zero_()
        self.phi.data.zero_()
        self.n = rot_z @ rot_y @ self._get_n()

    def _get_n(self) -> torch.Tensor:
        return torch.stack([torch.cos(self.beta) * torch.cos(self.phi), torch.sin(self.beta) * torch.cos(self.phi),
                            torch.sin(self.phi)]).reshape(3, 1).requires_grad_(True)

    def _get_tmp_params(self, s):
        beta_temp, phi_temp, theta_temp, d_temp = torch.remainder(self.beta - s * self.beta.grad, 2 * pi), \
                                                  torch.remainder(self.phi - s * self.phi.grad, 2 * pi), \
                                                  self.theta - s * self.theta.grad, \
                                                  self.d - s * self.d.grad
        self.tmp_params = {'beta': beta_temp, 'phi': phi_temp, 'theta': theta_temp, 'd': d_temp}

    def _update_params(self, s):
        # self.beta = torch.remainder(self.beta - s * self.beta.grad, 2 * pi)
        # self.phi = torch.remainder(self.phi - s * self.phi.grad, 2 * pi)
        with torch.no_grad():
            self.beta -= s * self.beta.grad
            self.phi -= s * self.phi.grad
            self.theta -= s * self.theta.grad
            self.d -= s * self.d.grad
            self._set_grad_to_zero()

    def _set_grad_to_zero(self):
        self.beta.grad.zero_()
        self.phi.grad.zero_()
        self.theta.grad.zero_()
        self.d.grad.zero_()
