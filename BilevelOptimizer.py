import torch
from LineSearcher import LineSearcher
from ConvexhullSettings import ConvexHullSettings
from typing import Callable

# define pi in torch
pi = torch.acos(torch.zeros(1)).item() * 2
data_type = torch.double


class BilevelOptimizer(object):
    def __init__(self, ch_settings: ConvexHullSettings, line_searcher: LineSearcher, obj_func: Callable,
                 gamma=torch.tensor(0.01, dtype=data_type)):
        self.ch_settings = ch_settings
        self.line_searcher = line_searcher
        self.func = obj_func
        self.gamma = gamma
        self.objectives = []

    def optimize(self, niters: int = 100000, tol: float = 1e-10, tolg: float = 1e-5, scale=0.9, invscale=2.):
        for i in range(niters):
            self.objectives.append(self.line_searcher.output)
            parameters = self.line_searcher.line_search(tol=tol, scale=scale, invscale=invscale)
            self.ch_settings.reset_parameters(parameters, chart_reset=True)
            parameters = self.ch_settings.get_params()
            parameters.append(self.gamma)
            self.line_searcher.reset_parameters(self.ch_settings.get_params())
            print(f"iter {i + 1}", self.line_searcher.jacobian_matrix)
            if self.line_searcher.grad_norm() < tolg:
                print("Converged!")
                break
