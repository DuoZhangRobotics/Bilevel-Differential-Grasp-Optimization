import torch
from convex_hulls import ConvexHulls
import numpy as np
from numpy.linalg import inv

pi = torch.acos(torch.zeros(1)).item() * 2


class BilevelOptimizer(object):
    def __init__(self, hull0: ConvexHulls, target: ConvexHulls):
        self.beta = torch.Tensor([0.0]).requires_grad_(True)
        self.phi = torch.Tensor([0.0]).requires_grad_(True)
        self.theta = torch.Tensor([0.0, 0.0, 0.0]).requires_grad_(True)
        self.hull0 = hull0
        self.hull1 = target
        self.distance, self.closest_pos0, self.closest_pos1 = hull0.distance_between_convex_hulls(target)
        self.d = torch.Tensor([self.distance / 10]).requires_grad_(True)
        self.n = torch.stack([torch.cos(self.beta) * torch.cos(self.phi), torch.sin(self.beta) * torch.cos(self.phi),
                              torch.sin(self.phi)]).reshape(3, 1).requires_grad_(True)

    def obj(self):
        points0 = torch.Tensor(self.hull0.points)
        #     ch1 = ConvexHull(points0 + theta)
        #     points1 = torch.Tensor(ch1.points)
        points1 = points0 + self.theta
        points2 = torch.Tensor(self.hull1.points)
        # vertices of two convex hulls
        #     v1 = points1[ch1.vertices, :].T
        v1 = points1[self.hull0.vertices, :]
        v2 = points2[self.hull1.vertices, :]
        # value of objective function
        objective = torch.sum(torch.log(torch.matmul(v2, self.n) - self.d)) + torch.sum(
            torch.log(self.d - torch.matmul(v1, self.n)))
        objective *= -1
        return objective

    def line_search(self, niter=200, tol=1e-4):
        pass

    def _reset_n(self):
        x_axis = np.array([1, 0, 0]).reshape((3, 1))
        rot_y_angle = np.pi / 4
        rot_z_angle = np.pi / 4
        beta = np.pi / 4
        phi =  np.pi / 3
        n = np.array([np.cos(beta) * np.cos(phi), np.sin(beta) * np.cos(phi), np.sin(phi)]).reshape((3, 1))
        rot_y = np.array([[np.cos(rot_y_angle), 0, -np.sin(rot_y_angle)],
                          [0, 1, 0],
                          [np.sin(rot_y_angle), 0, np.cos(rot_y_angle)]])
        rot_z = np.array([[np.cos(rot_z_angle), -np.sin(rot_z_angle), 0],
                          [np.sin(rot_z_angle), np.cos(rot_z_angle), 0],
                          [0, 0, 1]])
        print(n)
        rotation_matrix = inv(rot_z @ rot_y)
        print(rotation_matrix @ n)
