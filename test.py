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

    hull0 = ConvexHulls(cube)
    points2 = np.random.rand(30, 3) + 3
    hull2 = ConvexHulls(points2)

    optimizer = BilevelOptimizer(hull0, hull2)
    # optimizer.grad_check()
    # optimizer.reset_params()
    inputs = optimizer.get_params()

    # optimizer.find_jacobian()

    def obj(beta, phi, theta, d, v0, v2, centroid0, centroid1):
        n = torch.stack([torch.cos(beta) * torch.cos(phi), torch.sin(beta) * torch.cos(phi),
                         torch.sin(phi)]).reshape(3, 1).requires_grad_(True)
        v1 = v0 + theta
        # using the parameters from the optimizer
        objective = torch.sum(torch.log(torch.matmul(v2, n) - d)) + torch.sum(
            torch.log(d - torch.matmul(v1, n))) + torch.log(d)
        objective *= -1
        objective += torch.norm(centroid0 + theta - centroid1)
        return objective
        # using the tmp variables during line search


    check_result = torch.autograd.gradcheck(obj, inputs)
    print(f'check result is {check_result}')
    n = optimizer.n.detach().numpy()
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
