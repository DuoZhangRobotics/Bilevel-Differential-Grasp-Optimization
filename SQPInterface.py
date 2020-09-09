import torch
import scipy
import numpy as np
from MeritFunction import MeritFunction

data_type = torch.double


class QP(object):
    def __init__(self, function, constraints, u):
        self.function = function
        self.constraints = constraints
        self.u = u

    def lagrangian(self, x):
        return self.function(x) + self.u * self.constraints(x)

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

    def get_qp_objective(self, x, xk):
        d = x - xk
        derivative = torch.autograd.grad(self.function(xk), xk)[0]
        dLagrangian = torch.autograd.grad(self.lagrangian(xk), xk)[0]
        hessian = self.get_hessian(xk, dLagrangian)
        return derivative.T * d + 0.5 * d.T * hessian * d

    def get_qp_constraints(self, x, xk):
        d = x - xk
        jacobian_constraints = torch.autograd.grad(self.constraints(xk), xk)[0]
        return jacobian_constraints.T * d + self.constraints(xk)

    def get_derivatives(self, xk):
        derivative = torch.autograd.grad(self.function(xk), xk)[0]
        dLagrangian = torch.autograd.grad(self.lagrangian(xk), xk)[0]
        bk = self.get_hessian(xk, dLagrangian)
        dh = torch.autograd.grad(self.constraints(xk), xk)[0]
        return derivative, bk, dh

    @staticmethod
    def make_positive_definite(hessian,  min_cond=0.00001):
        eigenvalues, eigenvectors = np.linalg.eig(hessian)
        l = np.max(np.abs(eigenvalues))
        for i in range(len(eigenvalues)):
            eigenvalues[i] = max(eigenvalues[i], l * min_cond)
        return eigenvectors @ np.diag(eigenvalues) @ np.linalg.inv(eigenvectors)


class SQP(object):
    def __init__(self, function, constraints, mf: MeritFunction, u):
        self.function = function
        self.constraints = constraints
        self.merit_function = mf.merit_function
        self.u = u
        self.qp = QP(self.function, self.constraints, self.u)

    def lagrangian(self, x):
        return self.function(x) + self.u * self.constraints(x)

    def line_search(self, x, grad, direction, obj, c1=1e-4, s=1, tol=1e-20, scale=0.9):
        g_dot_d = (grad.T @ direction)[0, 0]

        new_x = x - s * torch.tensor(direction.T)
        new_obj = self.merit_function(new_x)
        # Armijo condition
        while new_obj > obj - c1 * s * g_dot_d or torch.isnan(new_obj) :
            s *= scale
            with torch.no_grad():
                new_x = x - s * torch.tensor(direction.T)
            new_x.requires_grad_(True)
            new_obj = self.merit_function(new_x)
            if s <= tol:
                return None, None, None
        return s, new_x, new_obj

    def solve(self, x0, u0, niters: int = 100000, tol: float = 1e-10, tolg: float = 1e-5, plot_interval=50):
        df, bk, dh = self.qp.get_derivatives(x0)
        u = u0
        x = x0
        for i in range(niters):


            if torch.norm

