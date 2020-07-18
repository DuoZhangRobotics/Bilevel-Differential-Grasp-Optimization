import plotly.graph_objects as go
import torch
import torch.nn as nn
from pprint import pprint
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from convex_hulls import ConvexHulls
import numpy as np
from optimizer import BilevelOptimizer

pi = torch.acos(torch.zeros(1)).item() * 2

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


    def control_cube(start: np.array, control_dict: {'d': 0, 'alpha': 0, 'beta': 0, 'gamma': 0}) -> np.array:
        end = start + control_dict['d']
        return end


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
    # fig.show()

    optimizer = BilevelOptimizer(hull, hull2)
    # print(optimizer.centroid1)
    # optimizer.line_search(niters=1500)
    # 
    # cube += optimizer.theta.detach().numpy()
    # hull = ConvexHulls(cube)
    # fig = plt.figure()
    # ax = fig.add_subplot(111, projection='3d')
    # ax.scatter3D(cube[:, 0], cube[:, 1], cube[:, 2], c='blue')
    # for simplex1 in hull.simplices:
    #     ax.plot3D(cube[simplex1, 0], cube[simplex1, 1], cube[simplex1, 2], 'orange')
    # 
    # points2 = np.random.rand(30, 3) + 3
    # hull2 = ConvexHulls(points2)
    # ax.scatter3D(points2[:, 0], points2[:, 1], points2[:, 2], c='tomato')
    # for simplex2 in hull2.simplices:
    #     ax.plot3D(points2[simplex2, 0], points2[simplex2, 1], points2[simplex2, 2], 'lightblue')
    # fig.show()

    x1 = np.linspace(0.1, 0.7, 100)
    x2 = np.linspace(1.3, 1.9, 100)
    theta = np.pi / 4
    w1 = np.cos(theta)
    w2 = np.sin(theta)
    mesh1 = np.meshgrid(x1, x1)
    mesh2 = np.meshgrid(x2, x2)

    figure = plt.figure()

