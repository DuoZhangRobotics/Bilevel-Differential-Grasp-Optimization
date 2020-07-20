import torch
from convex_hulls import ConvexHulls
import numpy as np
from numpy.linalg import inv

# define pi in torch
pi = torch.acos(torch.zeros(1)).item() * 2
data_type = torch.double


class BilevelOptimizer(object):
    def __init__(self, hull0: ConvexHulls, target: ConvexHulls):
        self.theta = torch.tensor([0, 0, 0], dtype=data_type).requires_grad_(True)
        self.hull0 = hull0
        self.hull1 = target
        self.distance, self.closest_pos0, self.closest_pos1 = hull0.distance_between_convex_hulls(target)
        self.d = torch.tensor([self.distance / 8], dtype=data_type).requires_grad_(True)
        self.tmp_params = None
        self.centroid0 = torch.tensor([torch.mean(torch.tensor(self.hull0.points[self.hull0.vertices, 0])),
                                       torch.mean(torch.tensor(self.hull0.points[self.hull0.vertices, 1])),
                                       torch.mean(torch.tensor(self.hull0.points[self.hull0.vertices, 2]))
                                       ], dtype=data_type)
        self.centroid1 = torch.tensor([torch.mean(torch.tensor(self.hull1.points[self.hull1.vertices, 0])),
                                       torch.mean(torch.tensor(self.hull1.points[self.hull1.vertices, 1])),
                                       torch.mean(torch.tensor(self.hull1.points[self.hull1.vertices, 2]))
                                       ], dtype=data_type)
        self._initialize_phi_beta()
        self.beta.requires_grad_(True)
        self.phi.requires_grad_(True)
        self.n = self._get_n()

    def reset_params(self):
        self.__init__(self.hull0, self.hull1)

    # calculate the value of objective function
    def obj(self, mode="Normal"):
        # using the parameters from the optimizer
        if mode == "Normal":
            points0 = torch.tensor(self.hull0.points, dtype=data_type)
            points1 = points0 + self.theta
            points2 = torch.tensor(self.hull1.points, dtype=data_type)
            # vertices of two convex hulls
            v1 = points1[self.hull0.vertices, :]
            v2 = points2[self.hull1.vertices, :]
            # value of objective function
            objective = torch.sum(torch.log(torch.matmul(v2, self.n) - self.d)) + torch.sum(
                torch.log(self.d - torch.matmul(v1, self.n))) + torch.log(self.d)
            objective *= -1
            objective += torch.norm(self.centroid0 + self.theta - self.centroid1)
            return objective
        # using the tmp variables during line search
        if mode == "Temp":
            beta: torch.tensor = self.tmp_params['beta']
            phi: torch.tensor = self.tmp_params['phi']
            theta: torch.tensor = self.tmp_params['theta']
            d: torch.tensor = self.tmp_params['d']
            beta.requires_grad_(True)
            phi.requires_grad_(True)
            theta.requires_grad_(True)
            d.requires_grad_(True)
            n = torch.stack(
                [torch.cos(beta) * torch.cos(phi), torch.sin(beta) * torch.cos(phi), torch.sin(phi)]).reshape(
                (3, 1))
            points0 = torch.tensor(self.hull0.points, dtype=data_type)
            #     ch1 = ConvexHull(points0 + theta)
            #     points1 = torch.tensor(ch1.points)
            points1 = points0 + theta
            points2 = torch.tensor(self.hull1.points, dtype=data_type)
            # vertices of two convex hulls
            #     v1 = points1[ch1.vertices, :].T
            v1 = points1[self.hull0.vertices, :]
            v2 = points2[self.hull1.vertices, :]
            # value of objective function
            objective = torch.sum(torch.log(torch.matmul(v2, n) - d)) + torch.sum(
                torch.log(d - torch.matmul(v1, n))) + torch.log(d)

            objective *= -1
            objective += torch.norm(self.centroid0 + theta - self.centroid1)
            objective.backward()
            return objective, beta, phi, theta, d

    def line_search(self, niters: int = 500, tol: float = 1e-10, scale=0.7, c1=1e-4, c2=0.8):
        result = self.obj()
        result.backward()

        for i in range(niters):
            print(f'Round: {i + 1}')
            old_obj = result

            print('last obj', old_obj)
            s = 1
            self._get_tmp_params(s)
            result_temp, _, _, _, _ = self.obj(mode='Temp')
            print('finding feasible step size...')
            # Armijo Condition
            while result_temp > old_obj - c1 * s * torch.sum(
                    torch.tensor(
                        [self.beta.grad, self.phi.grad, self.d.grad, torch.sum(self.theta.grad)], dtype=data_type)) or torch.isnan(
                result_temp):
                s *= scale
                self._get_tmp_params(s)
                result_temp, beta_temp, phi_temp, theta_temp, d_temp = self.obj(mode="Temp")

                if s <= tol:
                    break
            if s <= tol:
                print("step size is too small, line search terminated at armijo condition.")
                break
            print('s2 = ', s)
            # Curvature condition
            # while torch.sum(torch.tensor([beta_temp.grad, phi_temp.grad, d_temp.grad,
            #                               torch.sum(theta_temp.grad)])) < old_obj + c2 * s * torch.sum(torch.tensor(
            #     [self.beta.grad, self.phi.grad, self.d.grad, torch.sum(self.theta.grad)])) or torch.isnan(
            #     result_temp):
            #     s *= scale
            #     self._get_tmp_params(s)
            #     result_temp, beta_temp, phi_temp, theta_temp, d_temp = self.obj(mode="Temp")
            #     if s <= tol:
            #         break
            # if s <= tol:
            #     print("step size is too small, line search terminated at curvature condition.")
            #     break
            # print("s3 = ", s)
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

        rot_y = torch.tensor([[torch.cos(rot_y_angle), 0, torch.sin(rot_y_angle)],
                              [0, 1, 0],
                              [-torch.sin(rot_y_angle), 0, torch.cos(rot_y_angle)]
                              ], dtype=data_type).requires_grad_(True)
        rot_z = torch.tensor([[torch.cos(rot_z_angle), -torch.sin(rot_z_angle), 0],
                              [torch.sin(rot_z_angle), torch.cos(rot_z_angle), 0],
                              [0, 0, 1]
                              ], dtype=data_type).requires_grad_(True)
        self.beta.data.zero_()
        self.phi.data.zero_()
        self.n = rot_z @ rot_y @ self._get_n()

    def _get_n(self) -> torch.tensor:
        return torch.stack([torch.cos(self.beta) * torch.cos(self.phi), torch.sin(self.beta) * torch.cos(self.phi),
                            torch.sin(self.phi)]).reshape(3, 1).requires_grad_(True)

    def _get_tmp_params(self, s):
        beta_temp, phi_temp, theta_temp, d_temp = torch.remainder(self.beta - s * self.beta.grad,
                                                                  2 * pi).detach().clone(), \
                                                  torch.remainder(self.phi - s * self.phi.grad,
                                                                  2 * pi).detach().clone(), \
                                                  (self.theta - s * self.theta.grad).detach().clone(), \
                                                  (self.d - s * self.d.grad).detach().clone()
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

    def _initialize_phi_beta(self):
        closest_vec = torch.tensor(self.closest_pos1 - self.closest_pos0, dtype=data_type)
        closest_vec /= torch.norm(closest_vec)
        sin_phi = closest_vec[2]
        phi = torch.asin(sin_phi)
        beta = torch.atan2(closest_vec[1], closest_vec[0])
        self.phi = phi.detach().clone()
        self.beta = beta.detach().clone()

    def print_grad(self, n: int):
        print(f'++++++++++++++++++++++++++++++++++++++++++++ The {n}th iteration +++++++++++++++++++++++++++++++++++++')
        print(f'beta.grad = {self.beta.grad}')
        print(f'phi.grad = {self.phi.grad}')
        print(f'theta.grad = {self.theta.grad}')
        print(f'd.grad = {self.d.grad}')
        print('+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++')

    def get_params(self):
        points0 = torch.tensor(self.hull0.points, dtype=data_type)
        points1 = points0 + self.theta
        points2 = torch.tensor(self.hull1.points, dtype=data_type)
        v1 = points1[self.hull0.vertices, :]
        v2 = points2[self.hull1.vertices, :]
        return self.beta, self.phi, self.theta, self.d, v1, v2, self.centroid0, self.centroid1

    def _obj_grad_check(self, beta, phi, theta, d, v1, v2, centroid0, centroid1):
        n = torch.stack([torch.cos(beta) * torch.cos(phi), torch.sin(beta) * torch.cos(phi),
                         torch.sin(phi)]).reshape(3, 1).requires_grad_(True)
        # using the parameters from the optimizer
        objective = torch.sum(torch.log(torch.matmul(v2, n) - d)) + torch.sum(
            torch.log(d - torch.matmul(v1, n))) + torch.log(d)
        objective *= -1
        objective += torch.norm(centroid0 + theta - centroid1)
        return objective

    def grad_check(self, eps=1e-10):
        # pytorch numerical gradient
        result = self.obj()
        result.backward()
        self.print_grad(1)

        # pertubation of parameters
        _, _, _, _, v1, v2, centroid0, centroid1 = self.get_params()
        beta = self.beta.detach().clone() + eps
        phi = self.phi.detach().clone() + eps
        theta = self.theta.detach().clone() + eps
        d = self.d.detach().clone() + eps

        beta.requires_grad_(True)
        phi.requires_grad_(True)
        theta.requires_grad_(True)
        d.requires_grad_(True)
        result = self._obj_grad_check(beta, phi, theta, d, v1, v2, centroid0, centroid1)
        result.backward()
        print(beta.grad)
        print(phi.grad)
        print(theta.grad)
        print(d.grad)

