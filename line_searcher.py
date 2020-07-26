import torch
import numpy as np


class LineSearch(object):
    def __init__(self, func: callable, mode='Armijo', c1=1e-4, c2=0.9):
        """

        Parameters
        ----------
        func: the function going to be optimized
        mode: including 'Armijo' or 'Wolfe', if mode is set to Armijo, then only Wolfe first condition is applied.
              if the mode is 'Wolfe', then both the armijo condition and the curvature condition are applied.
        c1:   default is 1e-4, as the parameter for armijo condition
        c2:   optional. default value is 0.9, as the parameter for wolfe condition
        """
        self.func = func
        self.mode = mode
        self.c1, self.c2 = c1, c2
   