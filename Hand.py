import numpy as np
import vtk
from hand import Hand, vtk_render, vtk_add_from_hand, Link
import torch
import trimesh

# define pi in torch
pi = torch.acos(torch.zeros(1)).item() * 2
data_type = torch.double


class HandTarget(object):
    def __init__(self, hand: Hand, target: trimesh.Trimesh):
        self.hand = hand
        self.target = target
        self.centroid1 = self.target.centroid
        self.root = self.hand.palm
        self._initialize_params(self.root, self.target)
        self.chart_reset(self.root)
        if self.hand.use_eigen:
            self.params = torch.zeros((1, self.hand.extrinsic_size + self.hand.eg_num), 
                                      dtype=data_type).requires_grad_(True)
        else:
            self.params = torch.zeros((1, self.hand.extrinsic_size + self.hand.nr_dof()),
                                      dtype=data_type).requires_grad_(True)

    def _initialize_params(self, root: Link, target: trimesh.Trimesh):
        centroid0 = root.centroid
        centroid1 = torch.tensor(target.centroid, dtype=torch.double)
        closest_vec = torch.tensor(centroid1 - centroid0, dtype=torch.double)
        closest_vec /= torch.norm(closest_vec)
        sin_phi = closest_vec[2]
        root.phi = torch.asin(sin_phi)
        root.beta = torch.atan2(closest_vec[1], closest_vec[0])
        root.d = self._get_d(root, target)
        for child in root.children:
            self._initialize_params(child, target)

    def chart_reset(self, root: Link):
        rot_y_angle = -root.phi.detach()
        rot_z_angle = root.beta.detach()
        rot_y = torch.tensor([[torch.cos(rot_y_angle), 0, torch.sin(rot_y_angle)],
                              [0, 1, 0],
                              [-torch.sin(rot_y_angle), 0, torch.cos(rot_y_angle)]
                              ], dtype=data_type)
        rot_z = torch.tensor([[torch.cos(rot_z_angle), -torch.sin(rot_z_angle), 0],
                              [torch.sin(rot_z_angle), torch.cos(rot_z_angle), 0],
                              [0, 0, 1]
                              ], dtype=data_type)
        root.beta.data.zero_()
        root.phi.data.zero_()
        root.rotation_matrix @= rot_z @ rot_y
        for child in root.children:
            self.chart_reset(child)

    def _get_d(self, root: Link, target: trimesh.Trimesh):
        v0, v2 = torch.tensor(root.mesh.vertices, dtype=data_type), torch.tensor(target.vertices, dtype=data_type)
        n = self.get_n(root.rotation_matrix, root.beta, root.phi)
        lower_bound = torch.min(v2 @ n)
        upper_bound = torch.max(v0 @ n)
        d = torch.mean(torch.tensor([lower_bound, upper_bound]))
        return d.detach().clone()

    @staticmethod
    def get_n(rotation_matrix, beta, phi):
        return rotation_matrix @ (
            torch.stack([torch.cos(beta) * torch.cos(phi),
                         torch.sin(beta) * torch.cos(phi),
                         torch.sin(phi)]).reshape(3, 1).requires_grad_(True))

    def reset_parameters(self, params):
        self.params = params
        self.chart_reset(self.root)


if __name__ == '__main__':
    path = 'hand/ShadowHand/'
    scale = 0.01
    use_eigen = True
    hand = Hand(path, scale, use_joint_limit=True, use_quat=True, use_eigen=use_eigen)
    if hand.use_eigen:
        dofs = np.zeros(hand.eg_num)
    else:
        dofs = np.zeros(hand.nr_dof())
    hand.forward_kinematics(np.zeros(hand.extrinsic_size), dofs)

    hand.value_check(10)
    hand.grad_check(2)
    # hand.write_limits()
    renderer = vtk.vtkRenderer()
    vtk_add_from_hand(hand, renderer, scale)
    vtk_render(renderer, axes=False)
