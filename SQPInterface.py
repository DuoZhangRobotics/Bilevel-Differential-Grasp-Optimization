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
        return self.function(x)  # + self.u.T @ self.constraints(x)

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
        df: torch.tensor = torch.autograd.grad(self.function(xk), xk, retain_graph=True, create_graph=True)[0]
        dLagrangian = torch.autograd.grad(self.lagrangian(xk), xk, retain_graph=True, create_graph=True)[0]
        hl, _ = self.get_hessian(xk, dLagrangian)
        constraints = self.constraints(xk)
        dh: torch.tensor = torch.zeros((constraints.size(0), xk.size(1)))
        for i in range(constraints.size(0)):
            dh[i, :] = torch.autograd.grad(constraints[i], xk, retain_graph=True)[0]
        return df, hl, dh

    def solve(self, xk):
        df, hessianL, dh = self.get_derivatives(xk)
        df = df.detach().numpy()
        hessianL = hessianL.detach().numpy()
        dh = dh.detach().numpy()
        h = self.constraints(xk).detach().numpy()
        xk = xk.detach().numpy()
        dx = cp.Variable(xk.T.shape)
        # print(f'df={np.max(np.abs(df))} dh={np.max(np.abs(dh))} hl={np.max(np.abs(hl))} xk={np.max(np.abs(xk))}')
        # print("c = ", c)
        constraints = [dh @ dx + h <= 0]
        prob = cp.Problem(cp.Minimize(df @ dx + 0.5 * cp.quad_form(dx, hessianL)), constraints=constraints)
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
        # self.u = torch.ones((36, 1), dtype=torch.double)
        self.qp = QP(self.function, self.constraints)

    def lagrangian(self, x):
        return self.function(x)  # + self.u.T @ self.constraints(x)

    def solve(self, x0, niters: int = 100000, tol: float = 1e-10, tolg: float = 1e-5):
        s = 1.
        invscale = 2.
        scale = 0.9
        self.mf_values = []
        self.grad_norms = []
        self.objectives = []
        x: torch.tensor = x0
        # u: torch.tensor = self.u.detach().clone()
        dx, du = self.qp.solve(x)
        print("Start")
        self.mf = MeritFunction(self.function, self.constraints, x, dx, tol=tolg)
        mf_val = self.mf.merit_function(x)
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
                return None  # , None
            with torch.no_grad():
                x = new_x[0]
                # u = u - s * du
            x.requires_grad_(True)
            # u.requires_grad_(True)

            if s == last_s:
                s *= invscale

            # self.qp.u = u
            self.mf_values.append(mf_val.detach().numpy())
            self.objectives.append(self.lagrangian(x).detach().numpy())
            self.grad_norms.append(np.abs(self.mf.directional_derivative.detach().numpy()))
            dx_norm = np.max(np.abs(dx.detach().numpy()))
            dx, du = self.qp.solve(x)
            print("End")
            self.mf = MeritFunction(self.function, self.constraints, x, dx, tol=tolg)
            print(f"Iter{i:3d}: obj={self.objectives[-1]} grad={self.grad_norms[-1]} mf_val={mf_val} dx_norm={dx_norm} eta={self.mf.eta} s={s:3.6f}")
            # TODO: Convergence condition should not be mf.converged but KKT condition for primal problem. Must Introduce u to the structure.
            if self.mf.converged:
                print("SQP Converged!")
                break
        return x  # , u
