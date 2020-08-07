import torch
from LineSearcher import LineSearcher
from typing import Callable
import numpy as np
from numpy.linalg import LinAlgError
from HandTarget import HandTarget
import copy

# define pi in torch
pi = torch.acos(torch.zeros(1)).item() * 2
data_type = torch.double


class Optimizer(object):
    def __init__(self, obj_func: Callable, params, method='SGD', mode='Armijo',
                 c1=1e-4, c2=0.9, s=1, gamma=torch.tensor(0.01, dtype=data_type)):
        self.func = obj_func
        self.params = params
        self.output = self.func(*self.params)
        # self.hand_target_copy = copy.deepcopy(self.hand_target)
        print(f"OUTPUT = {self.output}")
        self.line_searcher = LineSearcher(self.func, self.params)
        self.objectives = []
        self.method = method
        self.mode = mode
        self.jacobian_matrix: torch.tensor = torch.autograd.grad(self.output,
                                                                 self.params[0],
                                                                 retain_graph=True,
                                                                 create_graph=True)[0]
        print("self.jacobian matrix = \n", self.jacobian_matrix)
        # print(self.jacobian_matrix.size())
        # self.grad_check()
        if self.method == "Newton":
            self.hessian_matrix = self.hessian()
            print(self.hessian_matrix)
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
                                                    direction=self.direction,
                                                    output=self.output,
                                                    mode=self.mode,
                                                    tol=tol,
                                                    scale=scale,
                                                    c1=self.c1,
                                                    c2=self.c2,
                                                    s=self.s)
            if self.s is not None:
                self.reset_parameters()
                self.s *= invscale
                print(f"iter {i + 1}", self.jacobian_matrix)
                if self.grad_norm() < tolg:
                    print("Converged!")
                    break
            else:
                break

    def hessian(self):
        length = self.jacobian_matrix.size(1)
        hessian_matrix = torch.zeros((length, length), dtype=data_type)
        for i in range(length):
            # with torch.autograd.detect_anomaly():
            # print(self.jacobian_matrix[:, i])
            hessian_matrix[i, :] = torch.autograd.grad(self.jacobian_matrix[:, i],
                                                       self.params[0],
                                                       retain_graph=True)[0]
        # check if the hessian matrix is positive definite
        try:
            np.linalg.cholesky(hessian_matrix.detach().numpy())
        except LinAlgError:
            # if the hessian matrix is not positive definite, ¬then make it positive definite.
            print("The Hessian matrix is not positive definite. Now trying to correct it...")
            hessian_matrix = torch.tensor(self.make_it_positive_definite(hessian_matrix), dtype=data_type)
        return hessian_matrix

    def make_it_positive_definite(self, hessian_matrix, scale=0.01):
        eigenvalues, eigenvectors = np.linalg.eig(hessian_matrix)
        l = np.linalg.norm(eigenvalues, ord=1)
        eigenvalues += scale * l
        hessian_matrix = eigenvectors @ np.diag(eigenvalues) @ np.linalg.inv(eigenvectors)
        try:
            np.linalg.cholesky(hessian_matrix)
        except LinAlgError:
            # if the hessian matrix is not positive definite, then make it positive definite.
            print("The Hessian matrix is not positive definite. Now trying to correct it...")
            hessian_matrix = self.make_it_positive_definite(hessian_matrix)
        return hessian_matrix

    def grad_norm(self):
        return torch.norm(self.jacobian_matrix)

    def find_direction(self):
        grad = self.jacobian_matrix.reshape((1, -1))
        direction = grad.T
        if self.method == 'Newton':
            direction = torch.inverse(self.hessian_matrix) @ direction

        return direction

    def reset_parameters(self):
        # print(f"Original Param0 = {self.params[0]}")
        param0 = self.params[0] - self.s * self.direction.T
        # print(f"New Param0 = {param0}")
        # print(f'hand params = {self.hand_target.params}')
        print("================reset=================")
        self.params[1].reset_parameters(param0, chart_reset=True)
        print('==================reset end ================')
        # print(f'hand params after reset = {self.hand_target.params}')
        self.params[0] = self.params[1].params
        # self.params[1] = self.hand_target
        # print(f"Changed Param0 = {self.params[0]}")
        self._reinit()

    def _reinit(self):
        print('reinit')
        self.output = self.func(*self.params)
        print(f"OUTPUT = {self.output}")
        self.line_searcher = LineSearcher(self.func, self.params)
        self.objectives = []
        self.jacobian_matrix: torch.tensor = torch.autograd.grad(self.output,
                                                                 self.params[0],
                                                                 retain_graph=True,
                                                                 create_graph=True)[0]
        if self.method == "Newton":
            self.hessian_matrix = self.hessian()
        self.direction = self.find_direction()

    def grad_check(self):
        print(torch.autograd.gradcheck(self.func, self.params))
