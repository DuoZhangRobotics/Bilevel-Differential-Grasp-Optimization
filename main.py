import torch
from convex_hulls import ConvexHulls
import numpy as np
from ConvexhullSettings import ConvexHullSettings
from LineSearcher import LineSearcher
from Objectives import obj_fun, hand_obj_fun
from Optimizer import Optimizer
from Plot import Plot
import trimesh
from hand import Hand, Link
from HandTarget import HandTarget
# define pi in torch
pi = torch.acos(torch.zeros(1)).item() * 2
data_type = torch.double

if __name__ == "__main__":
    path = r'./cube/'
    path = 'hand/ShadowHand/'
    scale = 0.01
    use_eigen = True
    hand = Hand(path, scale, use_joint_limit=True, use_quat=True, use_eigen=use_eigen)
    if hand.use_eigen:
        dofs = np.zeros(hand.eg_num)
        params = 100 * torch.ones((1, hand.extrinsic_size + hand.eg_num))
    else:
        dofs = np.zeros(hand.nr_dof())
        params = 100 * torch.ones((1, hand.extrinsic_size + hand.nr_dof()))
    hand.forward(params)
    cube = np.array([[-0.5, -0.5, -0.5],
                     [-0.5, 0.5, -0.5],
                     [0.5, -0.5, -0.5],
                     [0.5, 0.5, -0.5],
                     [-0.5, 0.5, 0.5],
                     [0.5, -0.5, 0.5],
                     [-0.5, -0.5, 0.5],
                     [0.5, 0.5, 0.5]])

    points2 = np.random.rand(30, 3) + 10.
    hull2 = ConvexHulls(points2)
    hull2 = ConvexHulls(hull2.points[hull2.vertices, :])
    target = trimesh.Trimesh(vertices=hull2.points[hull2.vertices, :], faces=hull2.simplices)

    gamma = torch.tensor(0.01, dtype=data_type)
    hand_target = HandTarget(hand, target)
    optimizer = Optimizer(hand_target, hand_obj_fun, params=[hand_target.params, hand_target, gamma], mode='Armijo', method='Newton')
    optimizer.optimize(niters=100)
    # plotter = Plot(bioptimizer)
    # plotter.plot_convex_hulls()
    # plotter.plot_obj()
    # plotter.plot_high_quality()
