import torch
import numpy as np
from typing import Callable
data_type = torch.double


class LineSearcher(object):
    def __init__(self, func: Callable, params: list):
        if not isinstance(func, Callable):
            raise TypeError("The input function should be callable.")
        if not isinstance(params, list):
            raise TypeError("Input parameters should be a list.")
        self.func = func
        self.params = params

    def line_search(self, grad, direction, obj, c1=1e-4, c2=0.9, s=1, tol=1e-20, scale=0.9):
        g_dot_d = (grad.T @ direction)[0, 0]
        param = self.params[0].clone()
        
        new_param = [param - s * torch.tensor(direction.T)] + self.params[1:3]
        new_obj = self.func(*new_param)
        # Armijo condition
        while new_obj > obj - c1 * s * g_dot_d or torch.isnan(new_obj):
            s *= scale
            with torch.no_grad():
                new_param = [param - s * torch.tensor(direction.T)] + self.params[1:3]
            new_param[0].requires_grad_(True)
            new_obj = self.func(*new_param)
            if s <= tol:
                return None, None, None

        return s, new_param, new_obj
