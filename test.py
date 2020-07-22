import plotly.graph_objects as go
import torch
import torch.nn as nn
from pprint import pprint
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from convex_hulls import ConvexHulls
import numpy as np
from optimizer import BilevelOptimizer
from torch import cos, sin, square

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

    # optimizer.find_jacobian()

    n = optimizer.get_n(optimizer.rotation_matrix, optimizer.beta, optimizer.phi).detach().numpy()
    optimizer.line_search(niters=1500)
    optimizer.plot_objective()

    cube += optimizer.theta.detach().numpy()
    hull = ConvexHulls(cube)
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter3D(cube[:, 0], cube[:, 1], cube[:, 2], c='blue')
    for simplex1 in hull.simplices:
        ax.plot3D(cube[simplex1, 0], cube[simplex1, 1], cube[simplex1, 2], 'orange')

    points2 = np.random.rand(30, 3) + 3
    hull2 = ConvexHulls(points2)
    ax.scatter3D(points2[:, 0], points2[:, 1], points2[:, 2], c='tomato')
    for simplex2 in hull2.simplices:
        ax.plot3D(points2[simplex2, 0], points2[simplex2, 1], points2[simplex2, 2], 'lightblue')

    direction = np.concatenate((np.array([0, 0, 0]).reshape((3, 1)), n), axis=1)
    ax.plot3D(direction[0, :], direction[1, :], direction[2, :], lw=5)
    fig.show()


    # x = torch.tensor([1, 0, 0], dtype=torch.double).reshape((3, 1))
    # beta = torch.tensor(0.4659915043137905, dtype=torch.double)
    # phi = torch.tensor(13.030464595970653, dtype=torch.double)
    # nn = torch.stack([torch.cos(beta) * torch.cos(phi), torch.sin(beta) * torch.cos(phi),
    #                   torch.sin(phi)]).reshape(3, 1).requires_grad_(True)
    # theta = torch.acos(x.T @ nn)
    # axis = torch.cross(x, nn)
    # ux, uy, uz = axis[0], axis[1], axis[2]
    # rotation_matrix = torch.tensor([[square(ux) + (1 - square(ux)) * cos(theta),
    #                                  ux * uy * (1 - cos(theta)) - uz * sin(theta),
    #                                  ux * uz * (1 - cos(theta) + uy * sin(theta))],
    #                                 [ux * uy * (1 - cos(theta)) + uz * sin(theta),
    #                                  square(uy) + (1 - square(uy)) * cos(theta),
    #                                  uz * uy * (1-cos(theta)) - ux * sin(theta)],
    #                                 [ux * uz * (1-cos(theta)) - uy * sin(theta),
    #                                  uz * uy * (1-cos(theta)) + ux * sin(theta),
    #                                  square(uz) + (1 - square(uz)) * cos(theta)]
    #                                 ])
    #
    # rot_y = torch.tensor([[cos(-phi), 0, sin(-phi)],
    #                       [0, 1, 0],
    #                       [-sin(-phi), 0, cos(-phi)]
    #                       ], dtype=torch.double)
    # rot_z = torch.tensor([[cos(beta), -sin(beta), 0],
    #                       [sin(beta), cos(beta), 0],
    #                       [0, 0, 1]
    #                       ], dtype=torch.double)
    # # result = axis * (axis.T @ x) + torch.cross(torch.cos(theta) * torch.cross(axis, x), axis) + torch.sin(theta) * torch.cross(axis, x)
    # result = rot_z @ rot_y @ x
    # print(nn)
    # print(result)
