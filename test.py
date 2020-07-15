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
    optimizer.line_search()
