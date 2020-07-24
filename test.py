import torch
import torch.nn as nn
from pprint import pprint
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from convex_hulls import ConvexHulls
from optimizer import BilevelOptimizer
from torch import cos, sin
import numpy as np
from plotting import plotting
import plotly.graph_objects as go

pi = torch.acos(torch.zeros(1, dtype=torch.double)) * 2

if __name__ == '__main__':
    cube = np.array([[-0.5, -0.5, -0.5],
                     [-0.5, 0.5, -0.5],
                     [0.5, -0.5, -0.5],
                     [0.5, 0.5, -0.5],
                     [-0.5, 0.5, 0.5],
                     [0.5, -0.5, 0.5],
                     [-0.5, -0.5, 0.5],
                     [0.5, 0.5, 0.5]])

    center = np.zeros((8, 3))

    hull0 = ConvexHulls(cube)

    points2 = np.random.rand(30, 3) + 3
    hull2 = ConvexHulls(points2)

    optimizer = BilevelOptimizer(hull0, hull2)
    inputs = optimizer.get_params()

    n = optimizer.get_n(optimizer.rotation_matrix, optimizer.beta, optimizer.phi).detach().numpy()
    optimizer.optimize(niters=1)
    # optimizer.plot_objective()
    #
    # cube += optimizer.theta.detach().numpy()
    # hull = ConvexHulls(cube)
    # fig = plt.figure()
    # ax = fig.add_subplot(111, projection='3d')
    # ax.scatter3D(cube[:, 0], cube[:, 1], cube[:, 2], c='blue')
    # for simplex1 in hull.simplices:
    #     ax.plot3D(cube[simplex1, 0], cube[simplex1, 1], cube[simplex1, 2], 'orange')
    #
    # ax.scatter3D(points2[:, 0], points2[:, 1], points2[:, 2], c='tomato')
    # for simplex2 in hull2.simplices:
    #     ax.plot3D(points2[simplex2, 0], points2[simplex2, 1], points2[simplex2, 2], 'lightblue')
    #
    # direction = np.concatenate((np.array([0, 0, 0]).reshape((3, 1)), n), axis=1)
    # ax.plot3D(direction[0, :], direction[1, :], direction[2, :], lw=5)
    # fig.savefig('plot.png')
    # fig.show()

    # plotting(cube, hull2.points)

