from Hand import Hand, vtk_render, vtk_add_from_hand, Link
from ConvexHulls import ConvexHull
import vtk,torch,trimesh
import numpy as np
data_type = torch.double

class HandTarget(object):
    def __init__(self, hand: Hand, target: ConvexHull):
        self.hand = hand
        self.target = target
        if self.hand.use_eigen:
            self.params = torch.zeros((1, self.hand.extrinsic_size + self.hand.eg_num + 3 * self.hand.link_num), dtype=data_type)
            self.front = self.hand.extrinsic_size + self.hand.eg_num
        else:
            self.params = torch.zeros((1, self.hand.extrinsic_size + self.hand.nr_dof() + 3 * self.hand.link_num), dtype=data_type)
            self.front = self.hand.extrinsic_size + self.hand.nr_dof()
        self.p, _ = self.hand.forward(torch.zeros((1, self.front)))
        self.rotation_matrix = torch.eye(3, dtype=data_type).repeat((self.hand.link_num, 1, 1))
        self._initialize_params(self.hand.palm, self.target, 0, 0)
        self.chart_reset(self.hand.palm, 0)
        self.params = self.params.detach().clone().requires_grad_(True)

    def _initialize_params(self, root: Link, target: ConvexHull, start, idx):
        hull = ConvexHull(self.p[:, :, start: start + len(root.mesh.vertices)].view(3, -1).T.numpy())
        d, c0, c1=hull.distance_to(target)
        
        beta, phi = self._get_beta_phi(c0, c1)
        n = self._get_n(self.rotation_matrix[idx, :, :], beta, phi)
        lower_bound0 = torch.min(torch.tensor(hull.surface_vertices()) @ n)
        upper_bound0 = torch.max(torch.tensor(hull.surface_vertices()) @ n)
        upper_bound1 = torch.max(torch.tensor(target.surface_vertices()) @ n)
        lower_bound1 = torch.min(torch.tensor(target.surface_vertices()) @ n)
        d = None
        if lower_bound0 >= upper_bound1:
            d = torch.mean(torch.tensor([lower_bound0, upper_bound1]))
        if lower_bound1 >= upper_bound0:
            d = torch.mean(torch.tensor([lower_bound1, upper_bound0]))

        self.params[0, self.front + 3 * idx + 0] = beta.detach().clone()
        self.params[0, self.front + 3 * idx + 1] = phi.detach().clone()
        self.params[0, self.front + 3 * idx + 2] = d.detach().clone()
        idx += 1
        start += len(root.mesh.vertices)
        for child in root.children:
            start, idx = self._initialize_params(child, target, start, idx)
        return start, idx

    def chart_reset(self, root: Link, idx):
        beta = self.params[0, self.front + 3 * idx + 0]
        phi  = self.params[0, self.front + 3 * idx + 1]
        rot_y_angle = -phi.detach()
        rot_z_angle = beta.detach()
        rot_y = torch.tensor([[torch.cos(rot_y_angle), 0, torch.sin(rot_y_angle)],
                              [0, 1, 0],
                              [-torch.sin(rot_y_angle), 0, torch.cos(rot_y_angle)]], dtype=data_type)
        rot_z = torch.tensor([[torch.cos(rot_z_angle), -torch.sin(rot_z_angle), 0],
                              [torch.sin(rot_z_angle), torch.cos(rot_z_angle), 0],
                              [0, 0, 1]], dtype=data_type)
        self.params[:, self.front + 3 * idx + 0].data.zero_()
        self.params[:, self.front + 3 * idx + 1].data.zero_()
        self.rotation_matrix[idx, :, :] = self.rotation_matrix[idx, :, :] @ rot_z @ rot_y
        idx += 1
        for child in root.children:
            idx = self.chart_reset(child, idx)
        return idx

    @staticmethod
    def _get_beta_phi(c0, c1):
        closest_vec = torch.tensor(c0-c1, dtype=data_type)
        closest_vec /= torch.norm(closest_vec)
        sin_phi = closest_vec[2]
        phi = torch.asin(sin_phi)
        beta = torch.atan2(closest_vec[1], closest_vec[0])
        return beta.detach().clone().requires_grad_(True), phi.detach().clone().requires_grad_(True)

    @staticmethod
    def _get_n(rotation_matrix, beta, phi):
        return rotation_matrix @ (
            torch.stack([torch.cos(beta) * torch.cos(phi),
                         torch.sin(beta) * torch.cos(phi),
                         torch.sin(phi)]).reshape(3, 1).requires_grad_(True))

    def reset_parameters(self, params, chart_reset=True):
        self.params = params
        if chart_reset:
            self.chart_reset(self.hand.palm, 0)
        self.params = self.params.detach().clone().requires_grad_(True)

    def get_log_barrier(self, root: Link, target: ConvexHull, params, p, start, idx):
        v0 = p[:, :, start: start + len(root.mesh.vertices)].view(3, -1).T
        v1 = torch.tensor(target.surface_vertices(), dtype=torch.double)
        
        beta = params[0, self.front + 3 * idx + 0]
        phi = params[0, self.front + 3 * idx + 1]
        d = params[0, self.front + 3 * idx + 2]
        
        n = self._get_n(self.rotation_matrix[idx, :, :], beta, phi)
        upper_bound0 = torch.max(v0 @ n)
        lower_bound0 = torch.min(v0 @ n)
        lower_bound1 = torch.min(v1 @ n)
        upper_bound1 = torch.max(v1 @ n)
        if lower_bound0 > upper_bound1:
            objective = -torch.sum(torch.log(v0 @ n - d)) - torch.sum(torch.log(d - v1 @ n))
        elif lower_bound1 > upper_bound0:
            objective = -torch.sum(torch.log(v1 @ n - d)) - torch.sum(torch.log(d - v0 @ n))
        else: objective = -torch.sum(torch.log(v1 @ n - d)) - torch.sum(torch.log(d - v0 @ n))
        start += len(root.mesh.vertices)
        idx += 1
        for child in root.children:
            tmp_objective, start, idx = self.get_log_barrier(child, target, params, p, start, idx)
            objective = objective + tmp_objective
        return objective, start, idx
    
    def get_norm(self, root: Link, target: ConvexHull, p, start):
        centroid1 = torch.tensor(target.centroid, dtype=torch.double)
        centroid0 = torch.mean(p[:, :, start: start + len(root.mesh.vertices)], dim=2)
        norm = torch.norm(centroid0 - centroid1)
        start += len(root.mesh.vertices)
        if root.children:
            for child in root.children:
                tmp_norm, start = self.get_norm(child, target, p, start)
                norm = norm + tmp_norm
        return norm, start
    
    def hand_target_objective(self, params, gamma):
        p, t = self.hand.forward(params[:,:self.front])
        norm, _ = self.get_norm(self.hand.palm, self.target, p, 0)
        objective, _, _ = self.get_log_barrier(self.hand.palm, self.target, params, p, 0, 0)
        objective = objective * gamma
        objective = objective + norm
        return objective