import torch
from convex_hulls import ConvexHulls
import numpy as np
from ConvexhullSettings import ConvexHullSettings
from LineSearcher import LineSearcher
from Objectives import obj_fun, hand_obj_fun
from Optimizer import Optimizer
from Plot import Plot
import trimesh
from Hand import Hand, Link, vtk_render, vtk_add_from_hand
from HandTarget import HandTarget
import vtk
# define pi in torch
pi = torch.acos(torch.zeros(1)).item() * 2
data_type = torch.double


if __name__ == "__main__":
    # path = r'./cube/'
    path = 'hand/ShadowHand/'
    scale = 0.01
    use_eigen = True
    hand = Hand(path, scale, use_joint_limit=True, use_quat=False, use_eigen=use_eigen)
    if hand.use_eigen:
        dofs = np.zeros(hand.eg_num)
        params = 100 * torch.zeros((1, hand.extrinsic_size + hand.eg_num))
    else:
        dofs = np.zeros(hand.nr_dof())
        params = 100 * torch.zeros((1, hand.extrinsic_size + hand.nr_dof()))
    hand.forward(params)
    cube = np.array([[-0.5, -0.5, -0.5],
                     [-0.5, 0.5, -0.5],
                     [0.5, -0.5, -0.5],
                     [0.5, 0.5, -0.5],
                     [-0.5, 0.5, 0.5],
                     [0.5, -0.5, 0.5],
                     [-0.5, -0.5, 0.5],
                     [0.5, 0.5, 0.5]])

    points2 = np.random.rand(30, 3) + 3.
    hull2 = ConvexHulls(points2)
    hull2 = ConvexHulls(hull2.points[hull2.vertices, :])
    target = trimesh.Trimesh(vertices=hull2.points[hull2.vertices, :], faces=hull2.simplices)

    gamma = torch.tensor(0.01, dtype=data_type)
    hand_target = HandTarget(hand, target)
    print("hand target rotation matrix: \n", hand_target.rotation_matrix)
    optimizer = Optimizer(hand_obj_fun, params=[hand_target.params, hand_target, gamma], mode='Armijo', method='GD')
    optimizer.optimize(niters=10000)
    params = optimizer.params[0][:, :hand_target.front].detach()
    hand.forward(params)
    renderer = vtk.vtkRenderer()
    vtk_add_from_hand(hand, renderer, 1.0, use_torch=True)
    vtk_render(renderer, axes=False)
    # plotter = Plot(bioptimizer)
    # plotter.plot_convex_hulls()
    # plotter.plot_obj()
    # plotter.plot_high_quality()
