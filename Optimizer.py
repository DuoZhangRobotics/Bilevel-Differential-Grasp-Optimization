from LineSearcher import LineSearcher
from typing import Callable
import numpy as np
from numpy.linalg import LinAlgError
from HandTarget import HandTarget
from Hand import vtk_add_from_hand, vtk_render
import time,copy,vtk,torch

# define pi in torch
pi = torch.acos(torch.zeros(1)).item() * 2
data_type = torch.double

class Optimizer(object):
    def __init__(self, obj_func: Callable, params, method='SGD', mode='Armijo',
                 c1=1e-4, c2=0.9, s=1., adaptive=True,
                 gamma=torch.tensor(0.1, dtype=data_type)):
        self.func = obj_func
        self.params = params
        self.adaptive = adaptive
        self.gamma = gamma
        if self.adaptive:
            self.params[-1] = self.gamma
        self.output = self.func(*self.params)
        self.last_output = torch.tensor(1e6, dtype=data_type)
        self.line_searcher = LineSearcher(self.func, self.params)
        self.objectives = []
        self.grad_norms = []
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
        self.meshes = [self.params[1].target, self.params[1].hand.draw(scale_factor=1, show_to_screen=False, use_torch=True)]

    def optimize(self, niters: int = 100000, tol: float = 1e-10, tolg: float = 1e-5, scale=0.9, invscale=2., plot_interval=50):
        for i in range(niters):
            self.objectives.append(self.output)
            grad = self.jacobian_matrix.reshape((1, -1))
            last_s = self.s
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
                if self.s==last_s:
                    self.s *= invscale
                if i % plot_interval == 0:
                    self.meshes.append(self.params[1].hand.draw(scale_factor=1, show_to_screen=False, use_torch=True))
                if self.grad_norm() < tolg:
                    print("Converged!")
                    break
            else:
                break
            print(f"Iter{i:3d}: obj={self.output:3.6f} grad={self.grad_norm():3.6f} s={self.s:3.6f}")

    def hessian(self):
        length = self.jacobian_matrix.size(1)
        hessian_matrix = torch.zeros((length, length), dtype=data_type)
        for i in range(length):
            # with torch.autograd.detect_anomaly():
            hessian_matrix[i, :] = torch.autograd.grad(self.jacobian_matrix[:, i], self.params[0], retain_graph=True)[0]
        # check if the hessian matrix is positive definite
        try:
            np.linalg.cholesky(hessian_matrix.detach().numpy())
        except LinAlgError:
            # if the hessian matrix is not positive definite, then make it positive definite.
            hessian_matrix = torch.tensor(self.make_positive_definite(hessian_matrix), dtype=data_type)
        return hessian_matrix

    def make_positive_definite(self, hessian_matrix, scale=1., min_cond=0.00001):
        eigenvalues, eigenvectors = np.linalg.eig(hessian_matrix)
        l = np.max(np.abs(eigenvalues))
        for i in range(len(eigenvalues)):
            eigenvalues[i] = max(eigenvalues[i], l * min_cond)
        hessian_matrix = eigenvectors @ np.diag(eigenvalues) @ np.linalg.inv(eigenvectors)
        try:
            np.linalg.cholesky(hessian_matrix)
        except LinAlgError:
            # if the hessian matrix is not positive definite, then make it positive definite.
            hessian_matrix = self.make_positive_definite(hessian_matrix)
        return hessian_matrix

    def grad_norm(self):
        #return torch.norm(self.jacobian_matrix)
        return torch.max(torch.abs(self.jacobian_matrix))

    def find_direction(self):
        grad = self.jacobian_matrix.reshape((1, -1))
        direction = grad.T
        if self.method == 'Newton':
            direction = torch.inverse(self.hessian_matrix) @ direction
        return direction

    def reset_parameters(self):
        with torch.no_grad():
            param0 = self.params[0] - self.s * self.direction.T
            self.params[1].reset_parameters(param0, chart_reset=True)
        self.params[0] = self.params[1].params.requires_grad_(True)
        if self.adaptive and torch.abs(self.output - self.last_output) < 10 * self.gamma and\
                self.gamma > torch.tensor(1e-4, dtype=data_type):
            self.gamma *= 0.1
            self.params[-1] = self.gamma
        self._reinit()

    def _reinit(self):
        self.last_output = self.output
        self.output = self.func(*self.params)
        self.grad_norms.append(self.grad_norm())
        self.line_searcher = LineSearcher(self.func, self.params)
        self.jacobian_matrix: torch.tensor = torch.autograd.grad(self.output,
                                                                 self.params[0],
                                                                 retain_graph=True,
                                                                 create_graph=True)[0]
        if self.method == "Newton":
            self.hessian_matrix = self.hessian()
        self.direction = self.find_direction()

    def grad_check(self):
        print(torch.autograd.gradcheck(self.func, self.params))
