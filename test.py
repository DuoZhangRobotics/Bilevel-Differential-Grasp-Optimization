import numpy as np
import vtk
from hand import Hand, vtk_render, vtk_add_from_hand, Link
import torch
import trimesh



path = r'./cube/'
path = 'hand/ShadowHand/'

hand = Hand(path, 1.0, use_joint_limit=False, use_quat=False, use_eigen=True, use_contacts=False)
if hand.use_eigen:
    dofs = np.zeros(hand.eg_num)
    params = 100 * torch.zeros((1, hand.extrinsic_size + hand.eg_num))
else:
    dofs = np.zeros(hand.nr_dof())
    params = 100 * torch.zeros((1, hand.extrinsic_size + hand.nr_dof()))
p, t = hand.forward(params)
print(hand.palm.mesh.vertices)
print(p[:, :, :len(hand.palm.mesh.vertices)])
print(p.shape)
print('====')
print(t.shape)
# hand.forward_kinematics(np.zeros(hand.extrinsic_size), dofs)
# hand.write_limits()
renderer = vtk.vtkRenderer()
vtk_add_from_hand(hand, renderer, 1.0, use_torch=True)
vtk_render(renderer, axes=False)
