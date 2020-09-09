import torch
import scipy
import numpy as np
from MeritFunction import MeritFunction

data_type = torch.double

class QP(object):
    def __init__(self, function, constraints):
        self.function = function
        self.constraints = constraints

    def get_derivative(self, x):
        return torch.autograd.grad(self.function, x, retain_graph=True, create_graph=True)[0]

    def get_hessian(self, x, jacobian):
        length = jacobian.size(1)
        hessian = torch.zeros((length, length), dtype=data_type)
        for i in range(length):
            # with torch.autograd.detect_anomaly():
            hessian[i, :] = torch.autograd.grad(jacobian[:, i], x, retain_graph=True)[0]
        try:
            # check if the hessian matrix is positive definite
            hessian = hessian.detach().numpy()
            return hessian, scipy.linalg.cho_factor(hessian)
        except np.linalg.LinAlgError:
            # if the hessian matrix is not positive definite, then make it positive definite.
            hessian = torch.tensor(self.make_positive_definite(hessian), dtype=data_type)
            hessian = hessian.detach().numpy()
            return hessian, scipy.linalg.cho_factor(hessian)

    def get_qp_objective(self, x):
        derivative = self.get_derivative(x)
        hessian = self.get_hessian(x, derivative)

    def make_positive_definite(self, hessian,  min_cond=0.00001):
        eigenvalues, eigenvectors = np.linalg.eig(hessian)
        l = np.max(np.abs(eigenvalues))
        for i in range(len(eigenvalues)):
            eigenvalues[i] = max(eigenvalues[i], l * min_cond)
        return eigenvectors @ np.diag(eigenvalues) @ np.linalg.inv(eigenvectors)

class SQP(object):
    def __init__(self, functions, constraints, merit_function: MeritFunction):
        self.functions = functions
        self.constraints = constraints
        self.merit_function = merit_function

    def get_qp(self, x):


    def get_jacobian(self, x):
        pass

    def get_hessian(self, x):
        pass

    def make_hessian_positive_definite(self, x):
        pass

    def line_search(self, x, merit_func):
        pass

    def solve(self, x0, u0):
        pass
