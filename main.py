import torch
from convex_hulls import ConvexHulls
import numpy as np
from ConvexhullSettings import ConvexHullSettings
from LineSearcher import LineSearcher
from function import obj_fun, hand_obj_fun
from Optimizer import Optimizer
from Plot import Plot
import trimesh
from hand import Hand, Link
from Hand import HandTarget
# define pi in torch
pi = torch.acos(torch.zeros(1)).item() * 2
data_type = torch.double

if __name__ == "__main__":
    path = 'hand/ShadowHand/'
    scale = 0.01
    use_eigen = True
    hand = Hand(path, scale, use_joint_limit=True, use_quat=True, use_eigen=use_eigen)
    if hand.use_eigen:
        dofs = np.zeros(hand.eg_num)
    else:
        dofs = np.zeros(hand.nr_dof())
    cube = np.array([[-0.5, -0.5, -0.5],
                     [-0.5, 0.5, -0.5],
                     [0.5, -0.5, -0.5],
                     [0.5, 0.5, -0.5],
                     [-0.5, 0.5, 0.5],
                     [0.5, -0.5, 0.5],
                     [-0.5, -0.5, 0.5],
                     [0.5, 0.5, 0.5]])

    mesh = trimesh.load_mesh(r'cube/off/cube.off').apply_scale(1.0)
    print(mesh.is_watertight)
    # mesh.show()

    path = r'./cube/'
    hand = Hand(path, 1.0, use_joint_limit=False, use_quat=False, use_eigen=False)
    if hand.use_eigen:
        dofs = np.zeros(hand.eg_num)
    else:
        dofs = np.zeros(hand.nr_dof())
    hand.forward_kinematics(np.zeros(hand.extrinsic_size), dofs)

    points2 = np.random.rand(30, 3) + 3
    hull2 = ConvexHulls(points2)
    hull2 = ConvexHulls(hull2.points[hull2.vertices, :])
    target = trimesh.Trimesh(vertices=hull2.points[hull2.vertices, :], faces=hull2.simplices)

    gamma = torch.tensor(0.01, dtype=data_type)
    hand_target = HandTarget(hand, target)
    optimizer = Optimizer(hand_target, hand_obj_fun, params=[hand_target.params, hand, target, gamma], mode='Armijo', method='Newton')
    optimizer.optimize(niters=10000)
    # plotter = Plot(bioptimizer)
    # plotter.plot_convex_hulls()
    # plotter.plot_obj()
    # plotter.plot_high_quality()
