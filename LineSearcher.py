import torch
import numpy as np
from typing import Callable, Union, Sequence
from numpy.linalg import LinAlgError

# define pi in torch
pi = torch.acos(torch.zeros(1)).item() * 2
data_type = torch.double


# def _as_list(inputs, arg_name, fn_name):
#     # Ensures that inp is a list of Tensors
#     # Returns whether or not the original inp was a tuple and the listed version of the input
#     is_input_list = True
#     if not isinstance(inputs, list):
#         inputs = (inputs,)
#         is_input_list = False
#
#     for i, el in enumerate(inputs):
#         if not torch.is_tensor(el):
#             if is_input_list:
#                 raise TypeError("The {} given to {} must be either a Tensor or a tuple of Tensors but the"
#                                 " value at index {} has type {}.".format(arg_name, fn_name, i, type(el)))
#             else:
#                 raise TypeError("The {} given to {} must be either a Tensor or a tuple of Tensors but the"
#                                 " given {} has type {}.".format(arg_name, fn_name, arg_name, type(el)))
#
#     return inputs
#
#
# _TensorOrTensors = Union[torch.tensor, Sequence[torch.tensor]]


class LineSearcher(object):
    def __init__(self, func: Callable, params: list):
        """

        Parameters
        ----------
        func: the function going to be optimized
        params: the inputs of func, parmas[0] should be the arg tensor shaped in 1xn that needs to be optimized
        """
        if not isinstance(func, Callable):
            raise TypeError("The input function should be callable.")

        if not isinstance(params, list):
            raise TypeError("Input parameters shoudle be a list.")

        self.func = func
        # self.params = _as_list(params, 'inputs', func.__name__)
        self.params = params
        self.params_need_to_be_optimized: torch.tensor = self.params[0]
        self.direction = torch.ones((1, self.params_need_to_be_optimized.size(1)))

    def line_search(self, grad, hessian_matrix, output, mode='Armijo', method='SGD', c1=1e-4, c2=0.9, s=1,
                    tol=1e-10, scale=0.9):
        """

        Parameters
        ----------

        s
        c2
        c1
        method
        mode
        output
        grad
        hessian_matrix
        tol:             tolerance for step size
        scale:           every time if the step size doesn't satisfy the wolfe condition, scale it

        Returns         None
        -------

        """

        # line search using Armijo Condition
        # Pre-calculation for Armijo condition, namely, the product of first order derivative and line searching
        line_searching_direction = grad.T
        if method == 'Newton':
            line_searching_direction = torch.inverse(hessian_matrix) @ line_searching_direction

        tmp_output, _ = self.update_output(s, line_searching_direction)
        while tmp_output > output - c1 * s * grad @ line_searching_direction or torch.isnan(tmp_output):
            s *= scale
            tmp_output, _ = self.update_output(s, line_searching_direction)
            if s <= tol:
                print("Step size is too small. Line search terminated when testing the armijo condition")
                break
        if mode == 'Wolfe':
            pass
        return s

    def update_output(self, s, line_search_direction):
        param0 = self.params_need_to_be_optimized - s * line_search_direction.T
        params = self.update_params(param0)
        return self.func(*params), param0

    def update_params(self, param0):
        params = self.params
        params[0] = param0
        return params
