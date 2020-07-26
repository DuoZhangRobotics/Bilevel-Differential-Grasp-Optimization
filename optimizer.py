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
        self.centroid0, self.centroid2 = self._initialize_centroids()
        self.beta, self.phi = self._initialize_phi_beta()
        self.rotation_matrix = torch.eye(3, dtype=data_type)
        self._reset_rotation_matrix()
        self.d = self._initialize_d().requires_grad_(True)
        self.objectives = []
        self.s = 1

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

    def optimize(self, niters: int = 100000, tol: float = 1e-10, tolg: float = 1e-5, scale=0.9, invscale=2., c1=1e-4):
        result = self.obj()
        result.backward()

        self.s = 1.
        for i in range(niters):
            old_obj = result
            self.objectives.append(old_obj)
            info = 'iter=%s obj0=%s ' % (str(i + 1), str(old_obj.detach().numpy()))
            info += self.print_grad()
            direction = torch.ones((6, 1))
            self._get_tmp_params(self.s)
            result_temp = self.obj(mode='Temp')
            # line search using Armijo Condition
            # Pre-calculation for Armijo condition, namely, the product of first order derivative and line searching
            grad = torch.tensor(
                [self.beta.grad,
                 self.phi.grad,,
                 self.theta.grad[0],
                 self.theta.grad[1],
                 self.theta.grad[2],
                 self.d.grad],
                dtype=data_type).reshape((1, -1))
            hessian = self.hessian()
            # print(np.linalg.cholesky(hessian.detach().numpy()))
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
            info += self.print_grad()
            self._update_params(self.s, line_searching_direction)
            info += self.print_params()
            info += self._reset_rotation_matrix()
            result = self.obj()
            result.backward()
            if i % 1 == 0:
                print(info)
            if self._grad_norm() < tolg:
                print("Converged")
                break
            self.s *= invscale

    # Hessian using torch.autograd.functional.hessian
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

    # Hessian using torch.autograd.grad twice
    @staticmethod
    def obj_hessian(params, rotation_matrix, v0, v2, centroid0, centroid1, gamma):
        beta = params[0, 0]
        phi = params[0, 1]
        theta = params[0, 2:5]
        d = params[0, 5]
        n = rotation_matrix @ (torch.stack([torch.cos(beta) * torch.cos(phi), torch.sin(beta) * torch.cos(phi),
                                            torch.sin(phi)]).reshape(3, 1).requires_grad_(True))
        v1 = v0 + theta
        # using the parameters from the optimizer
        objective = -gamma * torch.sum(torch.log(torch.matmul(v2, n) - d)) - gamma * torch.sum(
            torch.log(d - torch.matmul(v1, n)))
        objective += torch.norm(centroid0 + theta - centroid1)
        return objective

    def hessian(self):
        # beta, phi, theta, d, rotation_matrix, v0, v2, centroid0, centroid1 = self.get_params()
        # beta = beta.detach().clone().requires_grad_(True)
        # phi = phi.detach().clone().requires_grad_(True)
        # theta = theta.detach().clone().requires_grad_(True)
        # d = d.detach().clone().requires_grad_(True)
        # gamma = torch.tensor(0.01, dtype=data_type)
        # inputs = (beta, phi, theta, d, rotation_matrix, v0, v2, centroid0, centroid1, gamma)
        # hessian_matrix = torch.autograd.functional.hessian(self.obj_fun_hessian, inputs)
        # temp_hessian = []
        # for i in range(4):
        #     temp = []
        #     for j in range(4):
        #         if hessian_matrix[i][j].shape == ():
        #             temp.append(hessian_matrix[i][j].reshape((1, 1)))
        #         if hessian_matrix[i][j].shape == torch.Size([3]):
        #             if j == 2:
        #                 temp.append(hessian_matrix[i][j].reshape((1, 3)))
        #             if i == 2:
        #                 temp.append(hessian_matrix[i][j].reshape((3, 1)))
        #         if hessian_matrix[i][j].shape == torch.Size([3, 3]):
        #             temp.append(hessian_matrix[i][j])
        #     temp_hessian.append(torch.cat(temp, 1))
        # hessian_matrix = torch.cat(temp_hessian)

        """ split line """
        beta, phi, theta, d, rotation_matrix, v0, v2, centroid0, centroid1 = self.get_params()
        beta = beta.detach().clone().requires_grad_(True)
        phi = phi.detach().clone().requires_grad_(True)
        theta = theta.detach().clone().requires_grad_(True)
        d = d.detach().clone().requires_grad_(True)
        gamma = torch.tensor(0.01, dtype=data_type)
        params = torch.cat([beta.reshape((1, 1)), phi.reshape(1, 1), theta.reshape(1, 3), d.reshape((1, 1))],
                           1).requires_grad_(True)
        objective = self.obj_hessian(params, rotation_matrix, v0, v2, centroid0, centroid1, gamma)
        # print(torch.autograd.gradcheck(self.obj_hessian, (params, rotation_matrix, v0, v2, centroid0, centroid1, gamma)))
        jacobian = torch.autograd.grad(objective, params, retain_graph=True, create_graph=True)[0]
        print(f'jacobian = {jacobian}')
        length = jacobian.size(1)
        hessian_matrix = torch.zeros((length, length), dtype=data_type)
        for i in range(length):
            hessian_matrix[i, :] = torch.autograd.grad(jacobian[0, i], params, retain_graph=True)[0]
        return hessian_matrix

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

    def _get_tmp_params(self, s):
        beta_temp, phi_temp, theta_temp, d_temp = (self.beta - s * self.beta.grad).detach().clone(), \
                                                  (self.phi - s * self.phi.grad).detach().clone(), \
                                                  (self.theta - s * self.theta.grad).detach().clone(), \
                                                  (self.d - s * self.d.grad).detach().clone()
        self.tmp_params = {'beta': beta_temp, 'phi': phi_temp, 'theta': theta_temp, 'd': d_temp}

    def _update_params(self, s, direction):
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
        return phi.detach().clone().requires_grad_(True), beta.detach().clone().requires_grad_(True)

    def _initialize_centroids(self):
        centroid0 = self._get_centroid(self.hull0)
        centroid2 = self._get_centroid(self.hull2)
        return centroid0, centroid2

    @staticmethod
    def _get_centroid(hull):
        return torch.tensor([torch.mean(torch.tensor(hull.points[hull.vertices, 0])),
                             torch.mean(torch.tensor(hull.points[hull.vertices, 1])),
                             torch.mean(torch.tensor(hull.points[hull.vertices, 2]))
                             ], dtype=data_type)

    def _get_vertices(self):
        points0 = torch.tensor(self.hull0.points, dtype=data_type)
        points2 = torch.tensor(self.hull2.points, dtype=data_type)
        v0 = points0[self.hull0.vertices, :]
        v2 = points2[self.hull2.vertices, :]
        return v0, v2,,

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
