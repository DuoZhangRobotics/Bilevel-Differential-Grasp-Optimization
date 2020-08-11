import torch
import numpy as np
from ConvexhullSettings import ConvexHullSettings, ConvexHulls
from LineSearcher import LineSearcher
from Objectives import obj_fun, hand_obj_fun
from Optimizer import Optimizer
from Plot import Plot
import trimesh
from Hand import Hand, Link, vtk_render, vtk_add_from_hand
from HandTarget import HandTarget
import vtk
import matplotlib.pyplot as plt
import os
# define pi in torch
pi = torch.acos(torch.zeros(1)).item() * 2
data_type = torch.double


if __name__ == "__main__":
    path = 'hand/ShadowHand/'
    scale = 0.01
    use_eigen = True
    hand = Hand(path, scale, use_joint_limit=False, use_quat=False, use_eigen=use_eigen, use_contacts=False)
    if hand.use_eigen:
        dofs = np.zeros(hand.eg_num)
        params = torch.zeros((1, hand.extrinsic_size + hand.eg_num))
    else:
        dofs = np.zeros(hand.nr_dof())
        params = torch.zeros((1, hand.extrinsic_size + hand.nr_dof()))
    p, t = hand.forward(params)
    mesh = hand.draw(scale_factor=1, show_to_screen=False, use_torch=True)
    meshes = [mesh]

    cube = np.array([[-0.5, -0.5, -0.5],
                     [-0.5, 0.5, -0.5],
                     [0.5, -0.5, -0.5],
                     [0.5, 0.5, -0.5],
                     [-0.5, 0.5, 0.5],
                     [0.5, -0.5, 0.5],
                     [-0.5, -0.5, 0.5],
                     [0.5, 0.5, 0.5]])

    points2 = np.random.rand(30, 3) + np.array([0., -0.2, 1.3])
    hull2 = ConvexHulls(points2)
    hull2 = ConvexHulls(hull2.points[hull2.vertices, :])
    target = trimesh.Trimesh(vertices=hull2.points[hull2.vertices, :], faces=hull2.simplices)
    meshes.append(target)

    # renderer = vtk.vtkRenderer()
    # vtk_add_from_hand(meshes, renderer, 1.0, use_torch=True)
    #
    # vtk_render(renderer, axes=True)
    # exit(1)
    gamma = torch.tensor(0.001, dtype=data_type)
    hand_target = HandTarget(hand, target)
    optimizer = Optimizer(hand_obj_fun, params=[hand_target.params, hand_target, gamma], mode='Armijo',
                          method='Newton', adaptive=False, gamma=gamma)
    optimizer.optimize(niters=500, plot_interval=20)
    # params = optimizer.params[0][:, :hand_target.front].detach()
    # hand.forward(params)
    renderer = vtk.vtkRenderer()
    vtk_add_from_hand(optimizer.meshes, renderer, 1.0, use_torch=True)
    vtk_render(renderer, axes=True)
    meshes = [target]
    fig = plt.figure()
    plt.plot(optimizer.objectives, marker='d')
    plt.xlabel("iterations")
    plt.ylabel("value of objective function")
    fig.show()
    fig.savefig('Objective_function_values.png')

    fig = plt.figure()
    plt.plot(optimizer.grad_norms, marker='d')
    plt.xlabel("iterations")
    plt.ylabel("norm of gradients")
    fig.show()
    fig.savefig('grad_norm.png')
    # plotter = Plot(bioptimizer)
    # plotter.plot_convex_hulls()
    # plotter.plot_obj()
    # plotter.plot_high_quality()
