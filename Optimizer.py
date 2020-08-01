import torch
from LineSearcher import LineSearcher
from typing import Callable
import numpy as np
from numpy.linalg import LinAlgError
from Hand import HandTarget
# define pi in torch
pi = torch.acos(torch.zeros(1)).item() * 2
data_type = torch.double


class Optimizer(object):
    def __init__(self, hand_target: HandTarget, obj_func: Callable, params, method='SGD', mode='Armijo',
                 c1=1e-4, c2=0.9, s=1, gamma=torch.tensor(0.01, dtype=data_type)):
        self.hand_target = hand_target
        self.func = obj_func
        self.output = self.func(*params)
        self.params = params
        self.line_searcher = LineSearcher(self.func, self.params)
        self.params_need_to_be_optimized = self.params[0]
        print(self.params_need_to_be_optimized)
        self.objectives = []
        self.method = method
        self.mode = mode
        self.jacobian_matrix: torch.tensor = torch.autograd.grad(self.output,
                                                                 self.params_need_to_be_optimized,
                                                                 retain_graph=True,
                                                                 create_graph=True)[0]
        self.hessian_matrix = self.hessian()
        self.direction = self.find_direction()
        self.c1 = c1
        self.c2 = c2
        self.s = s
        self.gamma = gamma

    def optimize(self, niters: int = 100000, tol: float = 1e-10, tolg: float = 1e-5, scale=0.9, invscale=2.):
        for i in range(niters):
            self.objectives.append(self.output)
            grad = self.jacobian_matrix.reshape((1, -1))
            self.s = self.line_searcher.line_search(grad=grad,
                                                    hessian_matrix=self.hessian_matrix,
                                                    output=self.output,
                                                    mode=self.mode,
                                                    method=self.method,
                                                    tol=tol,
                                                    scale=scale,
                                                    c1=self.c1,
                                                    c2=self.c2,
                                                    s=self.s)

            self.reset_parameters()
            self.s *= invscale
            print(f"iter {i + 1}", self.jacobian_matrix)
            if self.grad_norm() < tolg:
                print("Converged!")
                break

    def hessian(self):
        length = self.jacobian_matrix.size(1)
        hessian_matrix = torch.zeros((length, length), dtype=data_type)
        for i in range(length):
            hessian_matrix[i, :] = torch.autograd.grad(self.jacobian_matrix[0, i], self.params[0], retain_graph=True)[0]
        # check if the hessian matrix is positive definite
        try:
            np.linalg.cholesky(hessian_matrix.detach().numpy())
        except LinAlgError:
            # if the hessian matrix is not positive definite, then make it positive definite.
            print("The Hessian matrix is not positive definite. Now trying to correct it...")
            hessian_matrix = torch.tensor(self.make_it_positive_definite(), dtype=data_type)
        return hessian_matrix

    def make_it_positive_definite(self, scale=0.01):
        eigenvalues, eigenvectors = np.linalg.eig(self.hessian_matrix)
        l = np.linalg.norm(eigenvalues, ord=1)
        eigenvalues += scale * l
        self.hessian_matrix = eigenvectors @ np.diag(eigenvalues) @ np.linalg.inv(eigenvectors)
        try:
            np.linalg.cholesky(self.hessian_matrix)
        except LinAlgError:
            # if the hessian matrix is not positive definite, then make it positive definite.
            print("The Hessian matrix is not positive definite. Now trying to correct it...")
            self.hessian_matrix = self.make_it_positive_definite()
        return self.hessian_matrix

    def grad_norm(self):
        return torch.norm(self.jacobian_matrix)

    def find_direction(self):
        grad = self.jacobian_matrix.reshape((1, -1))
        direction = grad.T
        if self.method == 'Newton':
            direction = torch.inverse(self.hessian_matrix) @ direction

        return direction

    def reset_parameters(self):
        param0 = self.params_need_to_be_optimized - self.s * self.direction.T
        self.params[0] = param0
        self.hand_target.reset_parameters(self.params)
        self.params = self.hand_target.params

    def get_params(self):
        beta =  self.params_need_to_be_optimized[0, 0]
        phi =   self.params_need_to_be_optimized[0, 1]
        theta = self.params_need_to_be_optimized[0, 2: 5]
        d =     self.params_need_to_be_optimized[0, 5]
        


