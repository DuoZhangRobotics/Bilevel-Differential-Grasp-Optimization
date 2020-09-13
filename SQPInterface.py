import torch
import scipy
import numpy as np
from MeritFunction import MeritFunction
from LineSearcher import LineSearcher

data_type = torch.double


class QP(object):
    def __init__(self, function, constraints, eta):
        self.function = function
        self.constraints = constraints
        self.eta = eta

    def lagrangian(self, x):
        # TODO: should be active sets or using slack variables
        return self.function(x) + self.eta * torch.sum(self.constraints(x))

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
            return torch.tensor(hessian, dtype=data_type), scipy.linalg.cho_factor(hessian)

    def get_qp_objective(self, x, xk):
        d = x - xk
        df = torch.autograd.grad(self.function(xk), xk)[0]
        dLagrangian = torch.autograd.grad(self.lagrangian(xk), xk)[0]
        hessian, _ = self.get_hessian(xk, dLagrangian)
        return df.T * d + 0.5 * d.T * hessian * d

    def get_qp_constraints(self, x, xk):
        d = x - xk
        dh = torch.autograd.grad(self.constraints(xk), xk)[0]
        return dh.T * d + self.constraints(xk)

    def get_derivatives(self, xk):
        df: torch.tensor = torch.autograd.grad(self.function(xk), xk, retain_graph=True, create_graph=True)[0]
        dLagrangian = torch.autograd.grad(self.lagrangian(xk), xk, retain_graph=True, create_graph=True)[0]
        bk, _ = self.get_hessian(xk, dLagrangian)
        constraints = self.constraints(xk).view((-1, 1))
        dh: torch.tensor = torch.zeros((constraints.size(0), xk.size(1)))
        for i in range(constraints.size(0)):
            dh[i, :] = torch.autograd.grad(constraints[i], xk, retain_graph=True)[0]

        return df, bk, dh

    def solve(self, xk):
        df, hl, dh = self.get_derivatives(xk)
        A = torch.cat((torch.cat((hl, -dh.T), dim=1),
                       torch.cat((dh, torch.zeros((dh.size(0), dh.size(0)))), dim=1)), dim=0)
        B = torch.cat((df.T, self.constraints(xk).view(-1, 1)))
        d, _ = torch.solve(-B, A)
        dx = d[:dh.size(1), :]
        du = d[dh.size(1):, :]
        return dx, du, df

    @staticmethod
    def make_positive_definite(hessian, min_cond=0.00001):
        eigenvalues, eigenvectors = np.linalg.eig(hessian)
        l = np.max(np.abs(eigenvalues))
        for i in range(len(eigenvalues)):
            eigenvalues[i] = max(eigenvalues[i], l * min_cond)
        return eigenvectors @ np.diag(eigenvalues) @ np.linalg.inv(eigenvectors)


class SQP(object):
    def __init__(self, function, constraints, mf: MeritFunction, u):
        self.function = function
        self.constraints = constraints
        self.mf = mf
        self.u = u
        self.qp = QP(self.function, self.constraints, self.u)

    def lagrangian(self, x):
        return self.function(x) + self.u * self.constraints(x)

    def solve(self, x0, niters: int = 100000, tol: float = 1e-10, tolg: float = 1e-5, plot_interval=50):
        s = 1.
        invscale = 2.
        scale = 0.9
        self.mf_values = []
        self.grad_norms = []
        self.objectives = []
        dx, du, _ = self.qp.solve(x0)
        x: torch.tensor = x0
        u: torch.tensor = self.u.detach().clone()
        mf_val = self.mf.merit_function(x)
        line_searcher = LineSearcher(self.mf.merit_function, [x])
        for i in range(niters):
            last_s = s
            s, new_x, mf_val = line_searcher.line_search(obj=mf_val,
                                                         directional_derivative=self.mf.get_directional_derivative(x,
                                                                                                                   dx),
                                                         direction=dx,
                                                         s=last_s, scale=scale, tol=tol,
                                                         use_directional_derivative=True)
            if s is None:
                print("Line-Search failed!")
                break
            with torch.no_grad():
                x = new_x[0]
                u = u - s * du
            x.requires_grad_(True)
            u.requires_grad_(True)

            if s == last_s:
                s *= invscale

            self.qp.u = u
            dx, du, df = self.qp.solve(x)
            self.grad_norms.append(torch.max(torch.abs(df)).detach().numpy())
            self.mf_values.append(mf_val.detach().numpy())
            self.objectives.append(self.lagrangian(x).detach().numpy())
            print(f"Iter{i:3d}: obj={self.objectives[-1]} grad={self.grad_norms[-1]} s={s:3.6f}")
            if self.grad_norms[-1] < tolg:
                print("Converged!")
        return x, u
