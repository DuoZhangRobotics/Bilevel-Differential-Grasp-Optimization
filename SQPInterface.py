import torch
import scipy
import numpy as np
from MeritFunction import MeritFunction
from LineSearcher import LineSearcher
import cvxpy as cp

data_type = torch.double


class QP(object):
    def __init__(self, function, constraints):
        self.function = function
        self.constraints = constraints
        # self.u = u

    def lagrangian(self, x):
        # TODO: should be active sets or using slack variables
        return self.function(x) # + self.u * torch.sum(self.constraints(x))

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

    def get_derivatives(self, xk):
        # df: torch.tensor = torch.autograd.grad(self.function(xk), xk, retain_graph=True, create_graph=True)[0]
        dLagrangian = torch.autograd.grad(self.lagrangian(xk), xk, retain_graph=True, create_graph=True)[0]
        bk, _ = self.get_hessian(xk, dLagrangian)
        constraints = self.constraints(xk)
        dh: torch.tensor = torch.zeros((constraints.size(0), xk.size(1)))
        for i in range(constraints.size(0)):
            dh[i, :] = torch.autograd.grad(constraints[i], xk, retain_graph=True)[0]
        return dLagrangian, bk, dh

    def solve(self, xk):
        dl, hl, dh = self.get_derivatives(xk)
        dl = dl.detach().numpy()
        hl = hl.detach().numpy()
        dh = dh.detach().numpy()
        dx = cp.Variable(xk.detach().numpy().T.shape)
        c = self.constraints(xk).detach().numpy()
        print("c = ", c)
        constraints = [dh @ dx + c <= 0]
        prob = cp.Problem(cp.Minimize(dl @ dx + 1 / 2 * cp.quad_form(dx, hl)), constraints=constraints)
        prob.solve()
        du = constraints[0].dual_value
        return torch.tensor(dx.value, dtype=data_type), torch.tensor(du, dtype=data_type)

    @staticmethod
    def make_positive_definite(hessian, min_cond=0.00001):
        eigenvalues, eigenvectors = np.linalg.eig(hessian)
        l = np.max(np.abs(eigenvalues))
        for i in range(len(eigenvalues)):
            eigenvalues[i] = max(eigenvalues[i], l * min_cond)
        return eigenvectors @ np.diag(eigenvalues) @ np.linalg.inv(eigenvectors)


class SQP(object):
    def __init__(self, function, constraints):
        self.function = function
        self.constraints = constraints
        # self.u = u
        self.qp = QP(self.function, self.constraints)

    def lagrangian(self, x):
        return self.function(x)  # + self.u * self.constraints(x)

    def solve(self, x0, niters: int = 100000, tol: float = 1e-10, tolg: float = 1e-5):
        s = 1.
        invscale = 2.
        scale = 0.9
        self.mf_values = []
        self.grad_norms = []
        self.objectives = []
        dx, du = self.qp.solve(x0)
        print("dx = ", dx)
        x: torch.tensor = x0
        # u: torch.tensor = self.u.detach().clone()
        self.mf = MeritFunction(self.function, self.constraints, x, dx, tol=tolg)
        mf_val = self.mf.merit_function(x)
        print("mf_val = ", mf_val)
        print("directional derivative = ", self.mf.directional_derivative)
        line_searcher = LineSearcher(self.mf.merit_function, [x])
        for i in range(niters):
            last_s = s
            s, new_x, mf_val = line_searcher.line_search(obj=mf_val,
                                                         directional_derivative=self.mf.directional_derivative,
                                                         direction=dx,
                                                         s=last_s, scale=scale, tol=tol,
                                                         use_directional_derivative=True)
            if s is None:
                print("Line-Search failed!")
                break
            with torch.no_grad():
                x = new_x[0]
                # u = u - s * du
            x.requires_grad_(True)
            # u.requires_grad_(True)

            if s == last_s:
                s *= invscale

            # self.qp.u = u
            dx, du = self.qp.solve(x)
            # self.grad_norms.append(torch.max(torch.abs(df)).detach().numpy())
            self.mf_values.append(mf_val.detach().numpy())
            self.objectives.append(self.lagrangian(x).detach().numpy())
            self.grad_norms.append(self.mf.directional_derivative)
            self.mf = MeritFunction(self.function, self.constraints, x, dx, tol=tolg)
            print(f"Iter{i:3d}: obj={self.objectives[-1]} s={s:3.6f}")
            if self.mf.converged:
                print("Converged!")
        return x
