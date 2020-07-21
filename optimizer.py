import torch
from convex_hulls import ConvexHulls
import numpy as np
import matplotlib.pyplot as plt

# define pi in torch
pi = torch.acos(torch.zeros(1)).item() * 2
data_type = torch.double


class BilevelOptimizer(object):
    def __init__(self, hull0: ConvexHulls, target: ConvexHulls):
        self.theta = torch.tensor([0, 0, 0], dtype=data_type).requires_grad_(True)
        self.hull0 = hull0
        self.hull1 = target
        self.distance, self.closest_pos0, self.closest_pos1 = hull0.distance_between_convex_hulls(target)
        # TODO: Change the initialization method of d
        # self.d = torch.tensor([self.distance / 8], dtype=data_type).requires_grad_(True)
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
        self.d = self._initialize_d().requires_grad_(True)
        self.objectives = []

    def reset_params(self):
        self.__init__(self.hull0, self.hull1)

    # TODO: Simplify obj function
    # calculate the value of objective function
    def obj(self, mode="Normal"):
        points0 = torch.tensor(self.hull0.points, dtype=data_type)
        v0 = points0[self.hull0.vertices, :]
        # using the parameters from the optimizer
        if mode == "Normal":
            points2 = torch.tensor(self.hull1.points, dtype=data_type)
            # vertices of two convex hulls
            v1 = v0 + self.theta
            v2 = points2[self.hull1.vertices, :]
            # value of objective function
            objective = torch.sum(torch.log(torch.matmul(v2, self.n) - self.d)) + torch.sum(
                torch.log(self.d - torch.matmul(v1, self.n)))
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
            #     ch1 = ConvexHull(points0 + theta)
            #     points1 = torch.tensor(ch1.points)
            points1 = points0 + theta
            points2 = torch.tensor(self.hull1.points, dtype=data_type)
            # vertices of two convex hulls
            #     v1 = points1[ch1.vertices, :].T
            v1 = v0 + theta
            v2 = points2[self.hull1.vertices, :]
            # value of objective function
            objective = torch.sum(torch.log(torch.matmul(v2, n) - d)) + torch.sum(
                torch.log(d - torch.matmul(v1, n)))

            objective *= -1
            objective += torch.norm(self.centroid0 + theta - self.centroid1)
            objective.backward()
            return objective, beta, phi, theta, d

    def line_search(self, niters: int = 500, tol: float = 1e-10, scale=0.9, c1=1e-4, c2=0.8):
        result = self.obj()
        result.backward()

        for i in range(niters):
            print(f'+++++++++++++++++++++++++++++++++++ The {i + 1}th iteration ++++++++++++++++++++++++++++++++++++++')
            old_obj = result
            self.objectives.append(old_obj)
            print('last obj', old_obj)
            self.print_grad()
            s = 1
            self._get_tmp_params(s)
            result_temp, _, _, _, _ = self.obj(mode='Temp')
            print('finding feasible step size...')
            # Armijo Condition
            # Pre-calculation for Armijo condition, namely, the product of first order derivative and line searching
            # direction
            partial_objective = torch.tensor(
                [self.beta.grad,
                 self.phi.grad,
                 self.d.grad,
                 self.theta[0],
                 self.theta[1],
                 self.theta[2]],
                dtype=data_type).reshape((1, -1))
            line_searching_direction = partial_objective.T
            # TODO: change Armijo condition
            while result_temp > old_obj - c1 * s * partial_objective @ line_searching_direction or \
                    torch.isnan(result_temp):
                s *= scale
                self._get_tmp_params(s)
                result_temp, beta_temp, phi_temp, theta_temp, d_temp = self.obj(mode="Temp")

                if s <= tol:
                    break
            if s <= tol:
                print("step size is too small, line search terminated at armijo condition.")
                break
            print('step size = ', s)
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
            print('+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++')
            if torch.abs(result - old_obj) < tol:
                print("Converged")
                break

    def find_jacobian(self):
        beta, phi, theta, d, v0, v2, centroid0, centroid1 = self.get_params()
        x = torch.tensor([beta.detach().clone(),
                          phi.detach().clone(),
                          theta[0].detach().clone(),
                          theta[1].detach().clone(),
                          theta[2].detach().clone(),
                          d.detach().clone()], dtype=data_type).requires_grad_(True)
        n = torch.stack([torch.cos(x[0]) * torch.cos(x[1]), torch.sin(x[0]) * torch.cos(x[1]),
                         torch.sin(x[1])]).reshape(3, 1).requires_grad_(True)
        v1 = v0 + torch.stack([x[2], x[3], x[4]])
        objective = torch.sum(torch.log(torch.matmul(v2, n) - x[-1])) + torch.sum(
            torch.log(x[-1] - torch.matmul(v1, n)))
        objective *= -1
        objective += torch.norm(centroid0 + torch.stack([x[2], x[3], x[4]]) - centroid1)
        J = torch.autograd.grad(objective, x, retain_graph=True)
        H = torch.autograd.grad(J, x)
        print(f"J = {J}")
        print(f'H = {H}')

    # reset beta and phi and n
    def _reset_n(self):
        rot_y_angle = -self.phi.detach()
        rot_z_angle = self.beta.detach()

        rot_y = torch.tensor([[torch.cos(rot_y_angle), 0, torch.sin(rot_y_angle)],
                              [0, 1, 0],
                              [-torch.sin(rot_y_angle), 0, torch.cos(rot_y_angle)]
                              ], dtype=data_type)
        rot_z = torch.tensor([[torch.cos(rot_z_angle), -torch.sin(rot_z_angle), 0],
                              [torch.sin(rot_z_angle), torch.cos(rot_z_angle), 0],
                              [0, 0, 1]
                              ], dtype=data_type)
        self.beta.data.zero_()
        self.phi.data.zero_()
        self.n = rot_z @ rot_y @ self._get_n()

    def _get_n(self) -> torch.tensor:
        return torch.stack([torch.cos(self.beta) * torch.cos(self.phi), torch.sin(self.beta) * torch.cos(self.phi),
                            torch.sin(self.phi)]).reshape(3, 1).requires_grad_(True)

    # TODO: Use unified function
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

    def _get_vertices(self):
        points0 = torch.tensor(self.hull0.points, dtype=data_type)
        points2 = torch.tensor(self.hull1.points, dtype=data_type)
        v0 = points0[self.hull0.vertices, :]
        v2 = points2[self.hull1.vertices, :]
        return v0, v2

    def _initialize_d(self):
        v0, v2 = self._get_vertices()
        lower_bound = torch.min(v2.mm(self.n))
        upper_bound = torch.max(v0.mm(self.n))
        d = torch.mean(torch.tensor([lower_bound, upper_bound]))
        return d.detach().clone()

    def print_grad(self):
        print(f'beta.grad = {self.beta.grad}')
        print(f'phi.grad = {self.phi.grad}')
        print(f'theta.grad = {self.theta.grad}')
        print(f'd.grad = {self.d.grad}')

    def get_params(self):
        v0, v2 = self._get_vertices()
        return self.beta, self.phi, self.theta, self.d, v0, v2, self.centroid0, self.centroid1

    def plot_objective(self):
        fig = plt.figure()
        plt.plot(self.objectives, marker='d')
        fig.show()

    @staticmethod
    def _obj_grad_check(beta, phi, theta, d, v0, v2, centroid0, centroid1):
        n = torch.stack([torch.cos(beta) * torch.cos(phi), torch.sin(beta) * torch.cos(phi),
                         torch.sin(phi)]).reshape(3, 1).requires_grad_(True)
        v1 = v0 + theta
        # using the parameters from the optimizer
        objective = torch.sum(torch.log(torch.matmul(v2, n) - d)) + torch.sum(
            torch.log(d - torch.matmul(v1, n)))
        objective *= -1
        objective += torch.norm(centroid0 + theta - centroid1)
        return objective

    # TODO: Remove this function
    def grad_check(self, eps=1e-4):
        # pytorch numerical gradient
        result = self.obj()
        result.backward()
        self.print_grad(1)

        # perturbation of parameters
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
        print()
