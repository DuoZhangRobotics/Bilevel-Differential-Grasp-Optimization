import torch
import numpy as np
from typing import Callable, Union, Sequence
from numpy.linalg import LinAlgError
from function import obj_fun
from optimizer import BilevelOptimizer
from ConvexhullSettings import ConvexHullSettings

# define pi in torch
pi = torch.acos(torch.zeros(1)).item() * 2
data_type = torch.double


def _as_list(inputs, arg_name, fn_name):
    # Ensures that inp is a tuple of Tensors
    # Returns whether or not the original inp was a tuple and the tupled version of the input
    is_input_list = True
    if not isinstance(inputs, list):
        inputs = (inputs,)
        is_input_list = False

    for i, el in enumerate(inputs):
        if not torch.is_tensor(el):
            if is_input_list:
                raise TypeError("The {} given to {} must be either a Tensor or a tuple of Tensors but the"
                                " value at index {} has type {}.".format(arg_name, fn_name, i, type(el)))
            else:
                raise TypeError("The {} given to {} must be either a Tensor or a tuple of Tensors but the"
                                " given {} has type {}.".format(arg_name, fn_name, arg_name, type(el)))

    return inputs, is_input_list


_TensorOrTensors = Union[torch.tensor, Sequence[torch.tensor]]


class LineSearcher(object):
    def __init__(self, func: Callable, params: _TensorOrTensors, mode='Armijo', c1=1e-4, c2=0.9, s=1):
        """

        Parameters
        ----------
        func: the function going to be optimized
        params: the inputs of func, parmas[0] should be the arg tensor shaped in 1xn that needs to be optimized
        mode: including 'Armijo' or 'Wolfe', if mode is set to Armijo, then only Wolfe first condition is applied.
              if the mode is 'Wolfe', then both the armijo condition and the curvature condition are applied.
        c1:   default is 1e-4, as the parameter for armijo condition
        c2:   optional. default value is 0.9, as the parameter for wolfe condition
        """
        if not isinstance(func, callable):
            raise TypeError("The input function should be callable.")

        if not isinstance(params, tuple):
            raise TypeError("Input parameters shouble be a tuple.")

        self.func = func
        self.params, _ = _as_list(params, 'inputs', func.__name__)
        self.params_need_to_be_optimized: torch.tensor = self.params[0]
        self.mode = mode
        self.c1, self.c2 = c1, c2
        self.s = s
        self.direction = torch.ones((1, self.params_need_to_be_optimized.size(1)))
        self.objective_values = []
        self.output = self.func(*params)
        self.jacobian_matrix: torch.tensor = torch.autograd.grad(self.output,
                                                                 self.params_need_to_be_optimized,
                                                                 retain_graph=True,
                                                                 create_graph=True)[0]
        self.hessian_matrix = self.hessian()

    def line_search(self, tol: float = 1e-10, scale=0.9, method='SGD'):
        """

        Parameters
        ----------

        tol:             tolerance for step size
        scale:           every time if the step size doesn't satisfy the wolfe condition, scale it
        method:          method to optimize the objective function, including "SGD" and "Newton"

        Returns         None
        -------

        """

        # line search using Armijo Condition
        # Pre-calculation for Armijo condition, namely, the product of first order derivative and line searching
        grad = self.jacobian_matrix.reshape((1, -1))
        line_searching_direction = grad.T
        if method == 'Newton':
            hessian = self.hessian()
            # check if the hessian matrix is positive definite
            try:
                np.linalg.cholesky(hessian.detach().numpy())
            except LinAlgError:
                # if the hessian matrix is not positive definite, then make it positive definite.
                print("The Hessian matrix is not positive definite. Now trying to correct it...")
                hessian = self.make_it_positive_definite(hessian)
            line_searching_direction = torch.inverse(hessian) @ line_searching_direction

        tmp_output, _ = self.update_output(line_searching_direction)
        while tmp_output > self.output - self.c1 * self.s * grad @ line_searching_direction or torch.isnan(tmp_output):
            self.s *= scale
            tmp_output, _ = self.update_output(line_searching_direction)
            if self.s <= tol:
                print("Step size is too small. Line search terminated when testing the armijo condition")
                break
        if self.mode == 'Wolfe':
            pass

        _, param_need_to_be_updated = self.update_output(line_searching_direction)
        return param_need_to_be_updated

    def hessian(self):
        length = self.jacobian_matrix.size(1)
        hessian_matrix = torch.zeros((length, length), dtype=data_type)
        for i in range(length):
            hessian_matrix[i, :] = torch.autograd.grad(self.jacobian_matrix[0, i], self.params, retain_graph=True)[0]
        return hessian_matrix

    def update_params(self, param0):
        params = self.params
        params[0] = param0
        return params

    def update_output(self, line_search_direction):
        param0 = self.params_need_to_be_optimized - self.s * line_search_direction.T
        params = self.update_params(param0)
        return self.func(*params), param0

    # TODO: fix this function later
    @staticmethod
    def make_it_positive_definite(matrix: torch.tensor):
        return matrix
