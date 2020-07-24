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
        self.hull2 = target
        self.distance, self.closest_pos0, self.closest_pos1 = hull0.distance_between_convex_hulls(target)
        self.tmp_params = None
        self.centroid0 = torch.tensor([torch.mean(torch.tensor(self.hull0.points[self.hull0.vertices, 0])),
                                       torch.mean(torch.tensor(self.hull0.points[self.hull0.vertices, 1])),
                                       torch.mean(torch.tensor(self.hull0.points[self.hull0.vertices, 2]))
                                       ], dtype=data_type)
        self.centroid2 = torch.tensor([torch.mean(torch.tensor(self.hull2.points[self.hull2.vertices, 0])),
                                       torch.mean(torch.tensor(self.hull2.points[self.hull2.vertices, 1])),
                                       torch.mean(torch.tensor(self.hull2.points[self.hull2.vertices, 2]))
                                       ], dtype=data_type)
        self._initialize_phi_beta()
        self.beta.requires_grad_(True)
        self.phi.requires_grad_(True)
        self.rotation_matrix = torch.eye(3, dtype=data_type)
        self._reset_rotation_matrix()
        self.d = self._initialize_d().requires_grad_(True)
        self.objectives = []
        self.s = 1

    def reset_params(self):
        self.__init__(self.hull0, self.hull2)

    def obj(self, mode="Normal"):
        # calculate the value of objective function
        v0, v2 = self._get_vertices()
        # using the parameters from the optimizer
        if mode == "Normal":
            # vertices of two convex hulls
            objective = self.obj_fun(self.beta, self.phi, self.theta, self.d, self.rotation_matrix, v0, v2,
                                     self.centroid0, self.centroid2)
            return objective
        # using the tmp variables during line search
        if mode == "Temp":
            beta = self.tmp_params['beta']
            phi = self.tmp_params['phi']
            theta = self.tmp_params['theta']
            d = self.tmp_params['d']
            beta.requires_grad_(True)
            phi.requires_grad_(True)
            theta.requires_grad_(True)
            d.requires_grad_(True)
            objective = self.obj_fun(beta, phi, theta, d, self.rotation_matrix, v0, v2, self.centroid0, self.centroid2)
            objective.backward()
            return objective

    def obj_fun(self, beta, phi, theta, d, rotation_matrix, v0, v2, centroid0, centroid1, gamma=0.01):
        n = self.get_n(rotation_matrix, beta, phi)
        v1 = v0 + theta
        # using the parameters from the optimizer
        objective = -gamma * torch.sum(torch.log(torch.matmul(v2, n) - d)) - gamma * torch.sum(
            torch.log(d - torch.matmul(v1, n)))
        objective += torch.norm(centroid0 + theta - centroid1)
        return objective

    def optimize(self, niters: int = 100000, tol: float = 1e-20, tolg: float = 1e-5, scale=0.9, invscale=2., c1=1e-4):
        result = self.obj()
        result.backward()

        self.s = 1.
        for i in range(niters):
            old_obj = result
            self.objectives.append(old_obj)
            info = 'iter=%s obj0=%s ' % (str(i + 1), str(old_obj.detach().numpy()))
            # info += self.print_grad()
            self._get_tmp_params(self.s)
            result_temp = self.obj(mode='Temp')

            # line search using Armijo Condition
            # Pre-calculation for Armijo condition, namely, the product of first order derivative and line searching
            grad = torch.tensor(
                [self.beta.grad,
                 self.phi.grad,
                 self.d.grad,
                 self.theta.grad[0],
                 self.theta.grad[1],
                 self.theta.grad[2]],
                dtype=data_type).reshape((1, -1))
            hessian = self.hessian()
            print((hessian[1][0]))
            break
            line_searching_direction = grad.T
            while result_temp > old_obj - c1 * self.s * grad @ line_searching_direction or torch.isnan(result_temp):
                self.s *= scale
                self._get_tmp_params(self.s)
                result_temp = self.obj(mode="Temp")
                if self.s <= tol:
                    break
            if self.s <= tol:
                print("Step size is too small, line search terminated at armijo condition.")
                break
            info += ' s=%s' % str(self.s)

            # adopt line search
            self._update_params(self.s)
            info += self.print_params()
            info += self._reset_rotation_matrix()
            result = self.obj()
            result.backward()
            if i % 1000 == 0:
                print(info)
            if self._grad_norm() < tolg:
                print("Converged")
                break
            self.s *= invscale

    @staticmethod
    def obj_fun_hessian(beta, phi, theta, d, rotation_matrix, v0, v2, centroid0, centroid1,
                        gamma=torch.tensor(0.01, dtype=data_type)):
        n = rotation_matrix @ (torch.stack([torch.cos(beta) * torch.cos(phi), torch.sin(beta) * torch.cos(phi),
                                            torch.sin(phi)]).reshape(3, 1).requires_grad_(True))
        v1 = v0 + theta
        # using the parameters from the optimizer
        objective = -gamma * torch.sum(torch.log(torch.matmul(v2, n) - d)) - gamma * torch.sum(
            torch.log(d - torch.matmul(v1, n)))
        objective += torch.norm(centroid0 + theta - centroid1)
        return objective

    def hessian(self):
        beta, phi, theta, d, rotation_matrix, v0, v2, centroid0, centroid1 = self.get_params()
        gamma = torch.tensor(0.01, dtype=data_type)
        inputs = (beta, phi, theta, d, rotation_matrix, v0, v2, centroid0, centroid1, gamma)
        # v0, v2 = self._get_vertices()
        # params = torch.cat([self.beta, self.phi, self.theta, self.d]).requires_grad_(True)
        hessian_matrix = torch.autograd.functional.hessian(self.obj_fun_hessian, inputs)
        temp_hessian = []
        for i in range(4):
            two = torch.stack(hessian_matrix[i][:2]).reshape((2, 1))
            three = torch.cat((two, hessian_matrix[i][2].reshape((3, 1))))
            print(torch.cat((three, hessian_matrix[i][3].reshape((1, 1)))))

        # objective = self.obj_fun(self.beta, self.phi, self.theta, self.d, self.rotation_matrix, v0, v2, self.centroid0,
        #                          self.centroid2)
        # jacobian_matrix = torch.autograd.grad(objective, params, create_graph=True, retain_graph=True)
        # print(jacobian_matrix)
        # hessian_matrix = torch.autograd.grad(jacobian_matrix, params)

        return hessian_matrix[:4]

    def _reset_rotation_matrix(self, show_info=False):
        # reset beta and phi and n
        info = ''
        n = self.get_n(self.rotation_matrix, self.beta, self.phi)
        if show_info:
            info += (' n0=%s' % str(n.detach().numpy().T))
            info += (' beta=%s, phi=%s' % (str(self.beta.detach().numpy()), str(self.phi.detach().numpy())))
        rot_y_angle = -self.phi.detach()
        rot_z_angle = self.beta.detach()
        # print('rot_y_angle = %s, rot_z_angle = %s'%(rot_y_angle.detach().numpy(),rot_z_angle.detach().numpy()))
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
        self.rotation_matrix @= rot_z @ rot_y
        n = self.get_n(self.rotation_matrix, self.beta, self.phi)
        if show_info:
            info += (' beta=%s, phi=%s' % (str(self.beta.detach().numpy()), str(self.phi.detach().numpy())))
            info += (' n1=%s' % str(n.detach().numpy().T))
        return info

    @staticmethod
    def get_n(rotation_matrix, beta, phi) -> torch.tensor:
        return rotation_matrix @ (torch.stack([torch.cos(beta) * torch.cos(phi), torch.sin(beta) * torch.cos(phi),
                                               torch.sin(phi)]).reshape(3, 1).requires_grad_(True))

    def _get_tmp_params(self, s):
        beta_temp, phi_temp, theta_temp, d_temp = (self.beta - s * self.beta.grad).detach().clone(), \
                                                  (self.phi - s * self.phi.grad).detach().clone(), \
                                                  (self.theta - s * self.theta.grad).detach().clone(), \
                                                  (self.d - s * self.d.grad).detach().clone()
        self.tmp_params = {'beta': beta_temp, 'phi': phi_temp, 'theta': theta_temp, 'd': d_temp}

    def _update_params(self, s):
        with torch.no_grad():
            self.beta -= s * self.beta.grad
            # self.beta = torch.remainder(self.beta, 2 * pi).detach().clone().requires_grad_(True)
            self.phi -= s * self.phi.grad
            # self.phi = torch.remainder(self.phi, 2 * pi).detach().clone().requires_grad_(True)
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
        points2 = torch.tensor(self.hull2.points, dtype=data_type)
        v0 = points0[self.hull0.vertices, :]
        v2 = points2[self.hull2.vertices, :]
        return v0, v2

    def _initialize_d(self):
        v0, v2 = self._get_vertices()
        n = self.get_n(self.rotation_matrix, self.beta, self.phi)
        lower_bound = torch.min(v2.mm(n))
        upper_bound = torch.max(v0.mm(n))
        d = torch.mean(torch.tensor([lower_bound, upper_bound]))
        return d.detach().clone()

    def _grad_norm(self):
        ret = 0.
        ret += np.linalg.norm(self.beta.grad.detach().numpy())
        ret += np.linalg.norm(self.phi.grad.detach().numpy())
        ret += np.linalg.norm(self.theta.grad.detach().numpy())
        ret += np.linalg.norm(self.d.grad.detach().numpy())
        return ret

    def print_grad(self):
        info = (' beta.grad=%s' % str(self.beta.grad.detach().numpy()))
        info += (' phi.grad=%s' % str(self.phi.grad.detach().numpy()))
        info += (' theta.grad=%s' % str(self.theta.grad.detach().numpy()))
        info += (' d.grad=%s' % str(self.d.grad.detach().numpy()))
        return info

    def print_params(self):
        info = (' beta=%s' % str(self.beta.detach().numpy()))
        info += (' phi=%s' % str(self.phi.detach().numpy()))
        info += (' theta=%s' % str(self.theta.detach().numpy()))
        info += (' d=%s' % str(self.d.detach().numpy()))
        return info

    def get_params(self):
        v0, v2 = self._get_vertices()
        return self.beta, self.phi, self.theta, self.d, self.rotation_matrix, v0, v2, self.centroid0, self.centroid2

    def plot_objective(self):
        fig = plt.figure()
        plt.plot(self.objectives, marker='d')
        plt.xlabel("iterations")
        plt.ylabel("value of objective function")
        fig.show()
