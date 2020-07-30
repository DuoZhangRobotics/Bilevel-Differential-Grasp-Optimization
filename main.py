import torch
from convex_hulls import ConvexHulls
import numpy as np
from ConvexhullSettings import ConvexHullSettings
from LineSearcher import LineSearcher
from function import obj_fun
from BilevelOptimizer import BilevelOptimizer
from Plot import Plot


# define pi in torch
pi = torch.acos(torch.zeros(1)).item() * 2
data_type = torch.double

if __name__ == "__main__":
    cube = np.array([[-0.5, -0.5, -0.5],
                     [-0.5, 0.5, -0.5],
                     [0.5, -0.5, -0.5],
                     [0.5, 0.5, -0.5],
                     [-0.5, 0.5, 0.5],
                     [0.5, -0.5, 0.5],
                     [-0.5, -0.5, 0.5],
                     [0.5, 0.5, 0.5]])

    hull0 = ConvexHulls(cube)

    points2 = np.random.rand(30, 3) + 3
    hull2 = ConvexHulls(points2)
    gamma = torch.tensor(0.01, dtype=data_type)
    ch_settings = ConvexHullSettings(hull0, hull2)
    parameters = ch_settings.get_params()
    parameters.append(gamma)
    line_searcher = LineSearcher(obj_fun, parameters, method='Newton')

    bioptimizer = BilevelOptimizer(ch_settings, line_searcher, obj_fun)
    bioptimizer.optimize(niters=10000)
    plotter = Plot(bioptimizer)
    plotter.plot_convex_hulls()
    plotter.plot_obj()
    plotter.plot_high_quality()
