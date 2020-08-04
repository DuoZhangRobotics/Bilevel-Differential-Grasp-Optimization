import torch
from LineSearcher import LineSearcher
from typing import Callable
import numpy as np
from numpy.linalg import LinAlgError
from HandTarget import HandTarget
from torchviz import make_dot, make_dot_from_trace

# define pi in torch
pi = torch.acos(torch.zeros(1)).item() * 2
data_type = torch.double
torch.autograd.set_detect_anomaly(True)


class Optimizer(object):
    def __init__(self, hand_target: HandTarget, obj_func: Callable, params, method='SGD', mode='Armijo',
                 c1=1e-4, c2=0.9, s=1, gamma=torch.tensor(0.01, dtype=data_type)):
        self.hand_target = hand_target
        self.func = obj_func
        self.output = self.func(*params)
        print(f"OUTPUT = {self.output}")
        self.params = params
        self.line_searcher = LineSearcher(self.func, self.params)
        self.objectives = []
        self.method = method
        self.mode = mode
        self.jacobian_matrix: torch.tensor = torch.autograd.grad(self.output,
                                                                 self.params[0],
                                                                 retain_graph=True,
                                                                 create_graph=True)[0]
        # self.grad_check()
        if self.method == "Newton":
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
            hessian_matrix[i, :] = torch.autograd.grad(self.jacobian_matrix[0, i],
                                                       self.params[0],
                                                       retain_graph=True)[0]
            print(hessian_matrix[i, :])
        # check if the hessian matrix is positive definite
        try:
            np.linalg.cholesky(hessian_matrix.detach().numpy())
        except LinAlgError:
            # if the hessian matrix is not positive definite, ¬then make it positive definite.
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
        param0 = self.params[0] - self.s * self.direction.T
        self.hand_target.reset_parameters(param0)
        self.params[0] = self.hand_target.params
        self._reinit()

    def _reinit(self):
        self.__init__(self.hand_target, self.func, self.params, self.method, self.mode, self.c1, self.c2,
                      self.s, self.gamma)

    def grad_check(self):
        print(torch.autograd.gradcheck(self.func, self.params))
