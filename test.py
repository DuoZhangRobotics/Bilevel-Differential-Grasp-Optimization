import numpy as np
import vtk
from Hand import Hand, vtk_render, vtk_add_from_hand, Link
import torch
import trimesh
from ConvexhullSettings import ConvexHullSettings, ConvexHulls
from Objectives import hand_obj_fun
from Optimizer import Optimizer
from Plot import Plot
path = 'hand/ShadowHand/'

hand = Hand(path, 0.01, use_joint_limit=True, use_quat=False, use_eigen=True, use_contacts=False)
if hand.use_eigen:
    dofs = np.zeros(hand.eg_num)
    params = 100 * torch.zeros((1, hand.extrinsic_size + hand.eg_num))
else:
    dofs = np.zeros(hand.nr_dof())
    params = 100 * torch.zeros((1, hand.extrinsic_size + hand.nr_dof()))
p, t = hand.forward(params)

palm = hand.palm.mesh.vertices


points2 = np.random.rand(30, 3) + np.array([1, 1, 1])
hull2 = ConvexHulls(points2)
hull2 = ConvexHulls(hull2.points[hull2.vertices, :])
target = trimesh.Trimesh(vertices=hull2.points[hull2.vertices, :], faces=hull2.simplices)

chsettings = ConvexHullSettings(points2, palm)
print(chsettings.d)
plotter = Plot(chsettings)
plotter.plot_high_quality()
plotter.plot_convex_hulls()



renderer = vtk.vtkRenderer()
mesh = [hand.draw(1, show_to_screen=False, use_torch=True)]
vtk_add_from_hand(mesh, renderer, 1.0, use_torch=True)
vtk_render(renderer, axes=False)
