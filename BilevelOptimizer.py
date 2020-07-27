import torch
import numpy as np
from LineSearcher import LineSearcher
from ConvexhullSettings import ConvexHullSettings
from typing import Callable
import matplotlib.pyplot as plt
from convex_hulls import ConvexHulls

# define pi in torch
pi = torch.acos(torch.zeros(1)).item() * 2
data_type = torch.double


class BilevelOptimizer(object):
    def __init__(self, ch_settings: ConvexHullSettings, line_searcher: LineSearcher, obj_func: Callable):
        self.ch_settings = ch_settings
        self.line_searcher = line_searcher
        self.func = obj_func
        self.objectives = []

    def optimize(self, niters: int = 100000, tol: float = 1e-10, tolg: float = 1e-5, scale=0.9, invscale=2.):
        for i in range(niters):
            self.objectives.append(self.line_searcher.output)
            parameters = self.line_searcher.line_search(tol=tol, scale=scale, invscale=invscale)
            self.ch_settings.reset_parameters(parameters, chart_reset=True)
            parameters = self.ch_settings.get_params()
            # parameters.append()
            self.line_searcher.reset_parameters(self.ch_settings.get_params())
            print(f"iter {i + 1}", self.line_searcher.jacobian_matrix)
            if self.line_searcher.grad_norm() < tolg:
                print("Converged!")
                break
            
    def plot_obj(self):
        fig = plt.figure()
        plt.plot(self.objectives, marker='d')
        plt.xlabel("iterations")
        plt.ylabel("value of objective function")
        fig.show()
        fig.savefig('Objective_function_values.')

    def plot_convex_hulls(self):
        hull, hull2 = self.ch_settings.get_hulls()
        theta = self.ch_settings.theta.detach().numpy()
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        ax.scatter3D(hull.points[:, 0], hull.points[:, 1], hull.points[:, 2], c='blue')
        for simplex1 in hull.simplices:
            ax.plot3D(hull.points[simplex1, 0], hull.points[simplex1, 1], hull.points[simplex1, 2], 'orange')

        hull.points += theta
        ax.scatter3D(hull.points[:, 0], hull.points[:, 1], hull.points[:, 2], c='orchid')
        for simplex1 in hull.simplices:
            ax.plot3D(hull.points[simplex1, 0], hull.points[simplex1, 1], hull.points[simplex1, 2], 'dodgerblue')

        ax.scatter3D(hull2.points[:, 0], hull2.points[:, 1], hull2.points[:, 2], c='tomato')
        for simplex2 in hull2.simplices:
            ax.plot3D(hull2.points[simplex2, 0], hull2.points[simplex2, 1], hull2.points[simplex2, 2], 'lightblue')
        fig.show()
        fig.savefig('Convex_hulls.png')
