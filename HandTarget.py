import numpy as np
import vtk
from Hand import Hand, vtk_render, vtk_add_from_hand, Link
import torch
import trimesh
from ConvexhullSettings import ConvexHullSettings
# define pi in torch
pi = torch.acos(torch.zeros(1)).item() * 2
data_type = torch.double


class HandTarget(object):
    def __init__(self, hand: Hand, target: trimesh.Trimesh):
        self.hand = hand
        self.target = target
        self.centroid1 = self.target.centroid
        if self.hand.use_eigen:
            self.params = torch.zeros((1, self.hand.extrinsic_size + self.hand.eg_num + 3 * self.hand.link_num),
                                      dtype=data_type)
            self.front = self.hand.extrinsic_size + self.hand.eg_num
        else:
            self.params = torch.zeros((1, self.hand.extrinsic_size + self.hand.nr_dof() + 3 * self.hand.link_num),
                                      dtype=data_type)
            self.front = self.hand.extrinsic_size + self.hand.nr_dof()

        self.p, _ = self.hand.forward(torch.zeros((1, self.front)))
        self.rotation_matrix = torch.eye(3, dtype=data_type).repeat((self.hand.link_num, 1, 1))
        self._initialize_params(self.hand.palm, self.target, 0, 0)
        self.chart_reset(self.hand.palm, 0)
        self._initialize_d(self.hand.palm, self.target, 0, 0)
        self.params = self.params.detach().clone().requires_grad_(True)

    def _initialize_params(self, root: Link, target: trimesh.Trimesh, start, idx):
        vertices = self.p[:, :, start: start + len(root.mesh.vertices)].view(3, -1).T
        chsettings = ConvexHullSettings(np.array(vertices), np.array(target.vertices), chart_reset=False)
        # closest_pos, closest_pos1 = chsettings.closest_pos0, chsettings.closest_pos1
        # closest_vec = torch.from_numpy(closest_pos1 - closest_pos)
        # closest_vec /= torch.norm(closest_vec)
        # sin_phi = closest_vec[2]
        #
        # phi = torch.asin(sin_phi)
        # beta = torch.atan2(closest_vec[1], closest_vec[0])

        # d = self._initialize_d(root, target, start)
        beta = chsettings.beta
        phi = chsettings.phi

        self.params[0, self.front + 3 * idx] = beta.detach().clone()
        self.params[0, self.front + 3 * idx + 1] = phi.detach().clone()
        # self.params[0, self.front + 3 * start + 2] = d.detach().clone()
        idx += 1
        start += len(root.mesh.vertices)
        for child in root.children:
            start, idx = self._initialize_params(child, target, start, idx)
        return start, idx

    def chart_reset(self, root: Link, start):
        beta, phi = self.get_beta_phi(start)
        rotation_matrix = self.get_rotation_matrix(start)
        # print(f'before beta = {beta}, phi = {phi}, n = {self.get_n(self.get_rotation_matrix(start), beta, phi).T}')
        rot_y_angle = -phi.detach()
        rot_z_angle = beta.detach()
        rot_y = torch.tensor([[torch.cos(rot_y_angle), 0, torch.sin(rot_y_angle)],
                              [0, 1, 0],
                              [-torch.sin(rot_y_angle), 0, torch.cos(rot_y_angle)]
                              ], dtype=data_type)
        rot_z = torch.tensor([[torch.cos(rot_z_angle), -torch.sin(rot_z_angle), 0],
                              [torch.sin(rot_z_angle), torch.cos(rot_z_angle), 0],
                              [0, 0, 1]
                              ], dtype=data_type)
        self.params[:, self.front + 3 * start].data.zero_()
        self.params[:, self.front + 3 * start + 1].data.zero_()
        # beta, phi = self.get_beta_phi(start)
        self.rotation_matrix[start, :, :] = rotation_matrix @ rot_z @ rot_y
        # print(f'after beta = {beta}, phi = {phi}, n = {self.get_n(self.get_rotation_matrix(start), beta, phi).T}')
        # print(f'after rotation matrix = \n{rotation_matrix}')
        start += 1
        for child in root.children:
            start = self.chart_reset(child, start)
        return start

    def _initialize_d(self, root: Link, target: trimesh.Trimesh, start, idx):
        v0 = self.p[:, :, start: start + len(root.mesh.vertices)].view(3, -1).T
        v2 = torch.tensor(target.vertices, dtype=data_type)
        beta, phi = self.get_beta_phi(idx)
        n = self.get_n(self.get_rotation_matrix(idx), beta, phi)
        lower_bound2 = torch.min(v2 @ n)
        upper_bound2 = torch.max(v2 @ n)
        upper_bound = torch.max(v0 @ n)
        lower_bound = torch.min(v0 @ n)
        # print(lower_bound, upper_bound)
        # print(lower_bound2, upper_bound2)
        d = None
        if lower_bound >= upper_bound2:
            d = torch.mean(torch.tensor([lower_bound, upper_bound2]))
        if lower_bound2 >= upper_bound:
            d = torch.mean(torch.tensor([lower_bound2, upper_bound]))
        self.params[0, self.front + 3 * idx + 2] = d.detach().clone()
        idx += 1
        start += len(root.mesh.vertices)
        for child in root.children:
            start, idx = self._initialize_d(child, target, start, idx)
        return start, idx

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

    def get_beta_phi(self, idx):
        beta = self.params[0, self.front + 3 * idx]
        phi = self.params[0, self.front + 3 * idx + 1]
        return beta, phi

    def get_d(self, idx):
        return self.params[0, self.front + 3 * idx + 2]

    def get_rotation_matrix(self, idx):
        return self.rotation_matrix[idx, :, :]


