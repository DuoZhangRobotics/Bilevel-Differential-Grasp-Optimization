import torch
import numpy as np
from typing import Callable, Union, Sequence
from function import obj_fun
from optimizer import BilevelOptimizer

# define pi in torch
pi = torch.acos(torch.zeros(1)).item() * 2
data_type = torch.double


def _as_tuple(inputs, arg_name, fn_name):
    # Ensures that inp is a tuple of Tensors
    # Returns whether or not the original inp was a tuple and the tupled version of the input
    is_input_tuple = True
    if not isinstance(inputs, tuple):
        inputs = (inputs,)
        is_input_tuple = False

    for i, el in enumerate(inputs):
        if not torch.is_tensor(el):
            if is_input_tuple:
                raise TypeError("The {} given to {} must be either a Tensor or a tuple of Tensors but the"
                                " value at index {} has type {}.".format(arg_name, fn_name, i, type(el)))
            else:
                raise TypeError("The {} given to {} must be either a Tensor or a tuple of Tensors but the"
                                " given {} has type {}.".format(arg_name, fn_name, arg_name, type(el)))

    return inputs, is_input_tuple


_TensorOrTensors = Union[torch.tensor, Sequence[torch.tensor]]


class LineSearcher(object):
    def __init__(self, bioptimizer: BilevelOptimizer, func: Callable, params: _TensorOrTensors, mode='Armijo', c1=1e-4,
                 c2=0.9):
        """

        Parameters
        ----------
        bioptimizer: the convex hull settings
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

        self.bioptimizer = bioptimizer
        self.func = func
        self.params, _ = _as_tuple(params, 'inputs', func.__name__)
        self.params_need_to_be_optimized: torch.tensor = self.params[0]
        self.mode = mode
        self.c1, self.c2 = c1, c2
        self.s = 1
        self.direction = torch.ones((1, self.params_need_to_be_optimized.size(1)))
        self.objective_values = []
        self.output = self.func(*params)
        self.jacobian_matrix: torch.tensor = torch.autograd.grad(self.output,
                                                                 self.params_need_to_be_optimized,
                                                                 retain_graph=True,
                                                                 create_graph=True)[0]
        self.hessian_matrix = self.hessian()

    def optimize(self, niters=10000, tol: float = 1e-10, tolg: float = 1e-5, scale=0.9, invscale=2., verbose=False):
        result: torch.tensor = self.func(*self.params)
        result.backward()

        for i in range(niters):
            old_obj = result
            self.objective_values.append(old_obj)
            result_temp = self.obj(mode='Temp')
            # line search using Armijo Condition
            # Pre-calculation for Armijo condition, namely, the product of first order derivative and line searching
            grad = self.params_need_to_be_optimized.grad.reshape((1, -1))
            hessian = self.hessian()
            # print(np.linalg.cholesky(hessian.detach().numpy()))
            line_searching_direction = grad.T
            while result_temp > old_obj - self.c1 * self.s * grad @ line_searching_direction or torch.isnan(
                    result_temp):
                self.s *= scale
                self._get_tmp_params(self.s)
                result_temp = self.obj(mode="Temp")
                if self.s <= tol:
                    break
            if self.s <= tol:
                print("Step size is too small, line search terminated at armijo condition.")
                break

            # adopt line search
            self._update_params(self.s, line_searching_direction)
            info = self._reset_rotation_matrix()
            result = self.obj()
            result.backward()
            if i % 1 == 0:
                print(info)
            if self._grad_norm() < tolg:
                print("Converged")
                break
            self.s *= invscale

    def hessian(self):
        length = self.jacobian_matrix.size(1)
        hessian_matrix = torch.zeros((length, length), dtype=data_type)
        for i in range(length):
            hessian_matrix[i, :] = torch.autograd.grad(self.jacobian_matrix[0, i], self.params, retain_graph=True)[0]
        return hessian_matrix
