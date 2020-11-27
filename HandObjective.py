from ConvexHulls import ConvexHull
from Metric import Metric
from Hand import Hand, Link
import numpy as np
import torch

np.set_printoptions(threshold=np.inf)
data_type = torch.double


class HandObjective(object):
    def __init__(self, hand: Hand, target: list, metric: Metric=None):
        self.hand = hand
        self.target = target
        self.metric = metric
        self.metric.setup_distance(self.hand)
        if self.hand.use_eigen:
            self.param_size = self.hand.extrinsic_size + self.hand.eg_num + 3 * self.hand.link_num * len(self.target)
            self.front = self.hand.extrinsic_size + self.hand.eg_num
        else:
            self.param_size = self.hand.extrinsic_size + self.hand.nr_dof() + 3 * self.hand.link_num * len(self.target)
            self.front = self.hand.extrinsic_size + self.hand.nr_dof()
        self.params = torch.zeros((1, self.param_size), dtype=data_type)
        self.p, _ = self.hand.forward(torch.zeros((1, self.front)))
        self.rotation_matrix = torch.eye(3, dtype=data_type).repeat((self.hand.link_num * len(self.target), 1, 1))
        self._initialize_params(self.hand.palm, self.target, 0, 0)
        self.chart_reset(self.hand.palm, 0)
        self.params = self.params.detach().clone().requires_grad_(True)

    def _initialize_params(self, root: Link, target: list, start, idx):
        hull = ConvexHull(self.p[:, :, start: start + len(root.mesh.vertices)].view(3, -1).T.numpy())
        for i, t in enumerate(target):
            d, c0, c1 = hull.distance_to(t)

            beta, phi = self._get_beta_phi(c0, c1)
            n = self.get_n(self.rotation_matrix[idx * len(target) + i, :, :], beta, phi)
            lower_bound0 = torch.min(torch.tensor(hull.surface_vertices()) @ n)
            upper_bound0 = torch.max(torch.tensor(hull.surface_vertices()) @ n)
            upper_bound1 = torch.max(torch.tensor(t.surface_vertices()) @ n)
            lower_bound1 = torch.min(torch.tensor(t.surface_vertices()) @ n)
            d = None
            if lower_bound0 >= upper_bound1:
                d = torch.mean(torch.tensor([lower_bound0, upper_bound1]))
            if lower_bound1 >= upper_bound0:
                d = torch.mean(torch.tensor([lower_bound1, upper_bound0]))

            self.params[0, self.front + 3 * (idx * len(target) + i) + 0] = beta.detach().clone()
            self.params[0, self.front + 3 * (idx * len(target) + i) + 1] = phi.detach().clone()
            self.params[0, self.front + 3 * (idx * len(target) + i) + 2] = d.detach().clone()
        idx += 1
        start += len(root.mesh.vertices)
        for child in root.children:
            start, idx = self._initialize_params(child, target, start, idx)
        return start, idx

    def chart_reset(self, root: Link, idx):
        for i, t in enumerate(self.target):
            beta = self.params[0, self.front + 3 * (idx * len(self.target) + i) + 0]
            phi = self.params[0, self.front + 3 * (idx * len(self.target) + i) + 1]
            rot_y_angle = -phi.detach()
            rot_z_angle = beta.detach()
            rot_y = torch.tensor([[torch.cos(rot_y_angle), 0, torch.sin(rot_y_angle)],
                                  [0, 1, 0],
                                  [-torch.sin(rot_y_angle), 0, torch.cos(rot_y_angle)]], dtype=data_type)
            rot_z = torch.tensor([[torch.cos(rot_z_angle), -torch.sin(rot_z_angle), 0],
                                  [torch.sin(rot_z_angle), torch.cos(rot_z_angle), 0],
                                  [0, 0, 1]], dtype=data_type)
            self.params[:, self.front + 3 * (idx * len(self.target) + i) + 0].data.zero_()
            self.params[:, self.front + 3 * (idx * len(self.target) + i) + 1].data.zero_()
            self.rotation_matrix[idx * len(self.target) + i, :, :] = self.rotation_matrix[idx * len(self.target) + i, :,
                                                                     :] @ rot_z @ rot_y
        idx += 1
        for child in root.children:
            idx = self.chart_reset(child, idx)
        return idx

    @staticmethod
    def _get_beta_phi(c0, c1):
        closest_vec = torch.tensor(c1 - c0, dtype=data_type)
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

    def reset_parameters(self, params, chart_reset=True):
        self.params = params
        if chart_reset:
            self.chart_reset(self.hand.palm, 0)
        self.params = self.params.detach().clone().requires_grad_(True)

    def get_log_barrier(self, root: Link, target: list, params, p, start, idx):
        v0 = p[:, :, start: start + len(root.mesh.vertices)].view(3, -1).T
        objective = 0
        for i, t in enumerate(target):
            v1 = torch.tensor(t.surface_vertices(), dtype=torch.double)

            beta = params[0, self.front + 3 * (idx * len(target) + i) + 0]
            phi = params[0, self.front + 3 * (idx * len(target) + i) + 1]
            d = params[0, self.front + 3 * (idx * len(target) + i) + 2]

            n = self.get_n(self.rotation_matrix[idx * len(target) + i, :, :], beta, phi)
            upper_bound0 = torch.max(v0 @ n)
            lower_bound0 = torch.min(v0 @ n)
            lower_bound1 = torch.min(v1 @ n)
            upper_bound1 = torch.max(v1 @ n)
            if lower_bound0 > upper_bound1:
                objective += -torch.sum(torch.log(v0 @ n - d)) - torch.sum(torch.log(d - v1 @ n))
            elif lower_bound1 > upper_bound0:
                objective += -torch.sum(torch.log(v1 @ n - d)) - torch.sum(torch.log(d - v0 @ n))
            else:
                objective += -torch.sum(torch.log(v1 @ n - d)) - torch.sum(torch.log(d - v0 @ n))

        start += len(root.mesh.vertices)
        idx += 1
        for child in root.children:
            tmp_objective, start, idx = self.get_log_barrier(child, target, params, p, start, idx)
            objective = objective + tmp_objective
        return objective, start, idx

    def get_norm(self, root: Link, target: list, p, start):
        centroid0 = torch.mean(p[:, :, start: start + len(root.mesh.vertices)], dim=2)
        norm = 0
        for t in target:
            centroid1 = torch.tensor(t.centroid, dtype=torch.double)
            norm = norm + torch.norm(centroid0 - centroid1)
        start += len(root.mesh.vertices)
        if root.children:
            for child in root.children:
                tmp_norm, start = self.get_norm(child, target, p, start)
                norm = norm + tmp_norm
        return norm, start

    def hand_objective(self, params, gamma, normCoef=1.):
        p, t = self.hand.forward(params[:, :self.front])
        
        #collision objective
        if gamma>0.:
            objective, _, _ = self.get_log_barrier(self.hand.palm, self.target, params, p, 0, 0)
            objective = objective * gamma
        else: 
            objective = 0.
            
        #debug objective, distance between hand and objective
        if normCoef>0.:
            norm, _ = self.get_norm(self.hand.palm, self.target, p, 0)
            objective = objective + norm * normCoef
        
        return objective

    def Q_metric_objective(self, params, gamma, alpha):
        p, t = self.hand.forward(params[:, :self.front])

        #collision objective
        if gamma>0.:
            objective, _, _ = self.get_log_barrier(self.hand.palm, self.target, params, p, 0, 0)
            objective = objective * gamma
        else: 
            objective = 0.

        # Q metric objective
        Q_metric = self.metric.compute_metric_torch(self.hand)
        
        objective = objective + Q_metric
        
        return objective
    
