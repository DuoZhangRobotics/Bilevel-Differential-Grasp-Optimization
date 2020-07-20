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


    hull0 = ConvexHulls(cube)
    points2 = np.random.rand(30, 3) + 3
    hull2 = ConvexHulls(points2)

    optimizer = BilevelOptimizer(hull0, hull2)
    optimizer.grad_check()
    optimizer.reset_params()
    inputs = optimizer.get_params()

    def obj(beta, phi, theta, d, v1, v2, centroid0, centroid1):
        n = torch.stack([torch.cos(beta) * torch.cos(phi), torch.sin(beta) * torch.cos(phi),
                         torch.sin(phi)]).reshape(3, 1).requires_grad_(True)
        # using the parameters from the optimizer
        objective = torch.sum(torch.log(torch.matmul(v2, n) - d)) + torch.sum(
            torch.log(d - torch.matmul(v1, n))) + torch.log(d)
        objective *= -1
        objective += torch.norm(centroid0 + theta - centroid1)
        return objective
        # using the tmp variables during line search


    torch.autograd.gradcheck(obj, inputs)

    optimizer.line_search(niters=1500)

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
    fig.show()




