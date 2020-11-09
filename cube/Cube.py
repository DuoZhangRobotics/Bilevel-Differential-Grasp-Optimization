from ConvexHulls import ConvexHull
import numpy as np
import torch
from SQPInterface import SQP
import cvxpy as cp


class Cube:
    def __init__(self, length):
        self.length = length
        self.cube = torch.tensor([[-self.length / 2, -self.length / 2, -self.length / 2],
                                  [-self.length / 2, self.length / 2, -self.length / 2],
                                  [self.length / 2, -self.length / 2, -self.length / 2],
                                  [self.length / 2, self.length / 2, -self.length / 2],
                                  [-self.length / 2, self.length / 2, self.length / 2],
                                  [self.length / 2, -self.length / 2, self.length / 2],
                                  [-self.length / 2, -self.length / 2, self.length / 2],
                                  [self.length / 2, self.length / 2, self.length / 2]], dtype=torch.double)

        self.current_pos = self.cube

    def forward(self, x):
        self.current_pos = self.cube + x
        return self.current_pos

    def convex_hull(self):
        return ConvexHull(self.current_pos.detach().numpy())


class CubeTarget:
    def __init__(self, cube: Cube, target):
        self.cube = cube
        self.target = target
        self.rotation_matrix = torch.eye(3, dtype=torch.double)
        self.params = torch.zeros((1, 3 + 2 + 1), dtype=torch.double).requires_grad_(True)
        self.initialize_params()
        self.chart_reset()

    def initialize_params(self):
        hull = self.cube.convex_hull()
        d, c0, c1 = hull.distance_to(self.target)
        beta, phi = self._get_beta_phi(c0, c1)
        n = self.get_n(self.rotation_matrix, beta, phi)
        lower_bound0 = torch.min(torch.tensor(hull.surface_vertices()) @ n)
        upper_bound0 = torch.max(torch.tensor(hull.surface_vertices()) @ n)
        upper_bound1 = torch.max(torch.tensor(self.target.surface_vertices()) @ n)
        lower_bound1 = torch.min(torch.tensor(self.target.surface_vertices()) @ n)
        d = None
        if lower_bound0 >= upper_bound1:
            d = torch.mean(torch.tensor([lower_bound0, upper_bound1]))
        if lower_bound1 >= upper_bound0:
            d = torch.mean(torch.tensor([lower_bound1, upper_bound0]))
        self.params[:, 3] = beta
        self.params[:, 4] = phi
        self.params[:, 5] = d

    def chart_reset(self):
        beta = self.params[0, 3]
        phi = self.params[0, 4]
        rot_y_angle = -phi.detach()
        rot_z_angle = beta.detach()
        rot_y = torch.tensor([[torch.cos(rot_y_angle), 0, torch.sin(rot_y_angle)],
                              [0, 1, 0],
                              [-torch.sin(rot_y_angle), 0, torch.cos(rot_y_angle)]], dtype=torch.double)
        rot_z = torch.tensor([[torch.cos(rot_z_angle), -torch.sin(rot_z_angle), 0],
                              [torch.sin(rot_z_angle), torch.cos(rot_z_angle), 0],
                              [0, 0, 1]], dtype=torch.double)
        self.params[:, 3].data.zero_()
        self.params[:, 4].data.zero_()
        self.rotation_matrix = self.rotation_matrix @ rot_z @ rot_y

    @staticmethod
    def _get_beta_phi(c0, c1):
        closest_vec = torch.tensor(c1 - c0, dtype=torch.double)
        closest_vec /= torch.norm(closest_vec)
        sin_phi = closest_vec[2]
        phi = torch.asin(sin_phi)
        beta = torch.atan2(closest_vec[1], closest_vec[0])
        return beta.detach().clone().requires_grad_(True), phi.detach().clone().requires_grad_(True)

    @staticmethod
    def get_n(rotation_matrix, beta, phi):
        return rotation_matrix @ (
            torch.stack([torch.cos(beta) * torch.cos(phi),
                         torch.sin(beta) * torch.cos(phi),
                         torch.sin(phi)]).reshape(3, 1).requires_grad_(True))

    def get_log_barrier(self, params, p):
        v0 = p
        v1 = torch.tensor(self.target.surface_vertices(), dtype=torch.double)

        beta = params[0, 3]
        phi = params[0, 4]
        d = params[0, 5]

        n = self.get_n(self.rotation_matrix, beta, phi)
        upper_bound0 = torch.max(v0 @ n)
        lower_bound0 = torch.min(v0 @ n)
        lower_bound1 = torch.min(v1 @ n)
        upper_bound1 = torch.max(v1 @ n)
        if lower_bound0 > upper_bound1:
            objective = -torch.sum(torch.log(torch.atan(v0 @ n - d))) - torch.sum(torch.log(torch.atan(d - v1 @ n)))
        elif lower_bound1 > upper_bound0:
            objective = -torch.sum(torch.log(torch.atan(v1 @ n - d))) - torch.sum(torch.log(torch.atan(d - v0 @ n)))
        else:
            objective = -torch.sum(torch.log(torch.atan(v1 @ n - d))) - torch.sum(torch.log(torch.atan(d - v0 @ n)))

        return objective

    def alg2_objective(self, params, gamma):
        p = self.cube.forward(params[:, :3])
        objective = self.get_log_barrier(params, p)
        objective = objective * gamma
        return objective

    def friction_cone_constraint(self, params, f, gamma, mu, scale_closeness=1.):
        n = self.get_n(self.rotation_matrix, self.params[0, 3], self.params[0, 4]).T
        closeness = self.closeness(params)
        closeness = closeness / scale_closeness
        lamb = gamma / closeness
        # constraints = (Q - torch.min(torch.sum(self.sampled_directions @ f.T, dim=1))).reshape((-1, 1))
        constraints = (n @ f.T - lamb).reshape((1, 1))
        friction = (torch.norm(torch.norm(f - n @ f.T * n)) - mu * n @ f.T).reshape((1, 1))
        constraints = torch.cat((constraints, friction), dim=0)
        return constraints

    def closeness(self, params):
        v0 = self.cube.forward(params[:, :3])
        t = torch.tensor(self.target.surface_vertices(), dtype=torch.double)
        n = self.get_n(self.rotation_matrix, self.params[0, 3], self.params[0, 4])
        d = self.params[0, 5]
        upper_bound0 = torch.max(v0 @ n)
        lower_bound0 = torch.min(v0 @ n)
        lower_bound1 = torch.min(t @ n)
        upper_bound1 = torch.max(t @ n)
        if lower_bound0 > upper_bound1:
            closeness = torch.sum(v0 @ n - d) + torch.sum(d - t @ n)
        elif lower_bound1 > upper_bound0:
            closeness = torch.sum(t @ n - d) + torch.sum(d - v0 @ n)
        else:
            closeness = torch.sum(t @ n - d) + torch.sum(d - v0 @ n)
        return closeness


class Solve(object):
    def __init__(self, cube_target: CubeTarget, sampled_directions, gamma1=torch.tensor(0.001, dtype=torch.double), mu=0.1):
        self.cube_target = cube_target
        self.sampled_directions = sampled_directions
        self.gamma1 = gamma1
        self.mu = mu

    def qf_solver(self):
        Q = cp.Variable(1)
        n = self.cube_target.get_n(self.cube_target.rotation_matrix, self.cube_target.params[0, 3], self.cube_target.params[0, 4])
        n = n.detach().numpy()
        n = n.T
        f = cp.Variable((1, 3))
        constraints = [Q <= cp.min(cp.sum(self.sampled_directions @ f.T, axis=1))]
        lamb = self.gamma1 / self.cube_target.closeness(self.cube_target.params)
        lamb = lamb.detach().numpy()
        constraints.append(n @ f.T <= lamb)
        constraints.append(cp.norm(f - n @ f.T @ n) <= self.mu * n @ f.T)
        prob = cp.Problem(cp.Maximize(Q), constraints)
        prob.solve()
        print(f"Q.value = {Q.value}, f.value = {f.value}")
        return torch.tensor(Q.value, dtype=torch.double), torch.tensor(f.value, dtype=torch.double)

    def obj_func(self, params):
        return self.cube_target.alg2_objective(params=params, gamma=self.gamma1)
        # return self.hand_target.hand_target_objective(params=params, gamma=self.gamma)

    def constraints_func(self, params):
        return self.cube_target.friction_cone_constraint(params=params, f=self.F, gamma=self.gamma1, mu=self.mu)

    def solve(self, x0, niters=100000, plot_interval=30):
        x = x0
        p = self.cube_target.cube.forward(x[:, :3])
        self.Q, self.F = self.qf_solver()
        self.old_Q = self.Q
        for i in range(niters):
            sqp_solver = SQP(self.obj_func, self.constraints_func)
            j = 1
            while self.gamma1 > 0.000000001:
                x = sqp_solver.solve(x, plot_interval=plot_interval)
                self.gamma1 *= 0.9
                print(f"Gamma shrinkage{j}: gamma={self.gamma1} x={x.detach().numpy()} dist={torch.norm(x[:, :3])}")
                j += 1
            # if x is None:
            #     print("SQP failed")
            #     break
            # cubetarget.reset_parameters(x, True)
            # p, _ = hand_target.hand.forward(x[:, :hand_target.front])
            # self.old_Q = self.Q
            # self.Q, self.F = self.qf_solver()
            # print(f'Alg2 Iter{i}: OBJ={self.obj_func(x)}')
            # if self.converged():
            #     print('Alg2 Converged!')
            #     self.x_optimal = x
            #     sqp_solver.plot_meshes()
            #     # self.u_optimal = u
            #     break
        return x  # , u


if __name__ == "__main__":
    cube = Cube(0.3)
    target = ConvexHull(np.array([[-0.3, -0.3, -0.3],
                                  [-0.3, 0.3, -0.3],
                                  [0.3, -0.3, -0.3],
                                  [0.3, 0.3, -0.3],
                                  [-0.3, 0.3, 0.3],
                                  [0.3, -0.3, 0.3],
                                  [-0.3, -0.3, 0.3],
                                  [0.3, 0.3, 0.3]]) + np.array([0., 0., 0.5]))
    cube_target = CubeTarget(cube, target)
    gamma1 = torch.tensor(0.001, dtype=torch.double)
    sampled_directions = np.array([[0, 0, 1]])
    solver = Solve(cube_target, sampled_directions, gamma1)
    solver.solve(cube_target.params)

