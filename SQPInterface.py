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
        return self.function(x) #+ self.u.T @ self.constraints(x)

    def get_hessian(self, x, jacobian):
        length = jacobian.size(1)
        hessian = torch.zeros((length, length), dtype=data_type)
        for i in range(length):
            # with torch.autograd.detect_anomaly():
            hessian[i, :] = torch.autograd.grad(jacobian[:, i], x, retain_graph=True)[0]

        try:
            # check if the hessian matrix is positive definite
            hessian = hessian.detach().numpy()
            # print("check hessian symmetric in try: ", self.check_symmetric(hessian))
            return torch.tensor(hessian, dtype=data_type), scipy.linalg.cho_factor(hessian)
        except np.linalg.LinAlgError:
            # if the hessian matrix is not positive definite, then make it positive definite.
            hessian = torch.tensor(self.make_positive_definite(hessian), dtype=data_type)
            hessian = hessian.detach().numpy()
            # print("check hessian symmetric in except: ", self.check_symmetric(hessian))
            return torch.tensor(hessian, dtype=data_type), scipy.linalg.cho_factor(hessian)

    def get_derivatives(self, xk):
        dLagrangian = torch.autograd.grad(self.lagrangian(xk), xk, retain_graph=True, create_graph=True)[0]
        hl, _ = self.get_hessian(xk, dLagrangian)
        constraints = self.constraints(xk)
        dh: torch.tensor = torch.zeros((constraints.size(0), xk.size(1)))
        for i in range(constraints.size(0)):
            dh[i, :] = torch.autograd.grad(constraints[i], xk, retain_graph=True)[0]
        return dLagrangian, hl, dh

    def solve(self, xk):
        dL, hessianL, dh = self.get_derivatives(xk)
        dL = dL.detach().numpy()
        hessianL = hessianL.detach().numpy()
        dh = dh.detach().numpy()
        h = self.constraints(xk).detach().numpy()
        xk = xk.detach().numpy()
        dx = cp.Variable(xk.T.shape)
        # print(f'df={np.max(np.abs(df))} dh={np.max(np.abs(dh))} hl={np.max(np.abs(hl))} xk={np.max(np.abs(xk))}')
        constraints = [dh @ dx + h <= 0]
        prob = cp.Problem(cp.Minimize(dL @ dx + 0.5 * cp.quad_form(dx, hessianL)), constraints=constraints)
        prob.solve(solver=cp.MOSEK, verbose=True)
        du = constraints[0].dual_value
        # print(f'x={xk} dL={dL}, hessianL={hessianL}, dh={dh.T}, h={h.T}, dx={dx.value.T}, du={du.T}')
        return torch.tensor(dx.value, dtype=data_type), torch.tensor(du, dtype=data_type)

    @staticmethod
    def make_positive_definite(hessian, min_cond=0.00001):
        eigenvalues, eigenvectors = np.linalg.eig(hessian)
        tau = min_cond - eigenvalues
        tau = np.where(tau < 0, 0, tau)
        return hessian + eigenvectors @ np.diag(tau) @ np.linalg.inv(eigenvectors)
        # for i in range(len(eigenvalues)):
        #     eigenvalues[i] = max(eigenvalues[i], min_cond)
        # return eigenvectors @ np.diag(eigenvalues) @ np.linalg.inv(eigenvectors)

    @staticmethod
    def check_symmetric(a, rtol=1e-05, atol=1e-08):
        return np.allclose(a, a.T, rtol=rtol, atol=atol)


class SQP(object):
    def __init__(self, function, constraints):
        self.function = function
        self.constraints = constraints

    def initialize_u(self, x):
        df = torch.autograd.grad(self.function(x), x, retain_graph=True, create_graph=True)[0]
        constraints = self.constraints(x)
        dh: torch.tensor = torch.zeros((constraints.size(0), x.size(1)))
        for i in range(constraints.size(0)):
            dh[i, :] = torch.autograd.grad(constraints[i], x, retain_graph=True)[0]

        # print("dh = ", dh)
        # print("dh dh.T = ", dh @ dh.T)
        # print("inverse of dh dh.T = ", torch.inverse(dh @ dh.T))
        # print("inverse  dh dh.T = ", torch.inverse(dh @ dh.T) @ dh)
        # print("df = ", df)

        # u0 = -torch.inverse(dh @ dh.T) @ dh @ df
        u0 = torch.ones(constraints.shape, dtype=data_type)
        u0[torch.where(constraints <= 0)] = 0
        return u0

    def solve(self, x0, niters: int = 100000, tol: float = 1e-30, tolg: float = 1e-5, hand_target=None, plot_interval=50):
        s = 1.
        invscale = 2.
        scale = 0.9
        self.mf_values = []
        self.grad_norms = []
        self.objectives = []
        if hand_target:
            self.meshes = [i.mesh() for i in hand_target.target]
            self.meshes.append(hand_target.hand.draw(scale_factor=1, show_to_screen=False, use_torch=True))
        x: torch.tensor = x0
        # u: torch.tensor = self.initialize_u(x)
        # self.u = u
        self.qp = QP(self.function, self.constraints)
        dx, du = self.qp.solve(x)
        self.mf = MeritFunction(self.function, self.constraints, x, dx, tol=tolg)
        mf_val = self.mf.merit_function(x)
        L = self.function(x) + du.T @ self.constraints(x)
        dL = torch.autograd.grad(L, x, retain_graph=True)[0]
        line_searcher = LineSearcher(self.mf.merit_function, [x])
        for i in range(niters):
            # self.grad_norms.append(np.abs(self.mf.directional_derivative.detach().numpy()))
            max_c = np.max(self.constraints(x).detach().numpy())
            print(f"Iter{i}: obj={L.detach().numpy()} dL={np.max(np.abs(dL.detach().numpy()))} x={x[:, :3].detach().numpy()} max_du={np.max(du.detach().numpy())} dfdx={self.mf.dfdx} mf_grad={self.mf.directional_derivative.detach().numpy()}  max_constraint={max_c} eta={self.mf.eta.detach().numpy()} s={s}")
            last_s = s
            s, new_x, mf_val = line_searcher.line_search(obj=mf_val,
                                                         directional_derivative=-self.mf.directional_derivative,
                                                         direction=-dx,
                                                         s=last_s, scale=scale, tol=tol,
                                                         use_directional_derivative=True)
            if s is None:
                if hand_target:
                    self.plot_meshes()
                print("Line-Search failed!")
                return None  # , None
            with torch.no_grad():
                x = new_x[0]
                # print("x = ", x[:, :hand_target.front])
                # u = u + s * du
                # cons = self.constraints(x)
                # u[cons <= 0] = 0
                # u[torch.where(cons > 1)] = torch.log(cons[torch.where(cons > 0)])
                # u[torch.where(0 < cons < 1)] = -torch.log(cons[torch.where(0 < cons < 1)])
                # u[cons > 0] = 10

            x.requires_grad_(True)
            # u.requires_grad_(True)

            # self.u = u
            if hand_target:
                hand_target.hand.forward(x[:, :hand_target.front])
                self.meshes.append(hand_target.hand.draw(scale_factor=1, show_to_screen=False, use_torch=True))
                if (i + 1) % plot_interval == 0:
                    self.plot_meshes()

            converged, L, dL = self.converge_condition(x, du, tolg)
            if converged:
                print(f"Iter{i + 1}: obj={L.detach().numpy()} dL={np.max(np.abs(dL.detach().numpy()))} max_du={np.max(du.detach().numpy())} x={x[:, :3].detach().numpy()} mf_grad={self.mf.directional_derivative.detach().numpy()}  max_constraint={max_c} eta={self.mf.eta.detach().numpy()} s={s}")
                print("SQP converged!")
                break

            self.qp = QP(self.function, self.constraints)
            dx, du = self.qp.solve(x)
            self.mf = MeritFunction(self.function, self.constraints, x, dx, tol=tolg, last_eta=self.mf.eta)
            mf_val = self.mf.merit_function(x)
            line_searcher = LineSearcher(self.mf.merit_function, [x])
            if s == last_s:
                s *= invscale
        if hand_target:
            self.plot_meshes()
        return x  # , u

    def converge_condition(self, x, u, tol=1e-9):
        L = self.function(x) + u.T @ self.constraints(x)
        dL = torch.autograd.grad(L, x, retain_graph=True)[0]
        u = u.detach().numpy()
        converged = True
        if not np.max(np.abs(dL.detach().numpy())) <= tol:
            converged = False

        constraints = self.constraints(x).detach().numpy()
        if not np.all(constraints <= 0):
            tmp_c = constraints[np.where(constraints > 0)]
            if not np.max(np.abs(tmp_c)) < tol:
                converged = False
        if not np.all(u >= 0):
            tmp_u = u[np.where(u < 0)]
            if not np.max(np.abs(tmp_u)) < tol:
                converged = False

        if not np.max(np.abs(u * constraints)) < tol:
            converged = False
        return converged, L, dL

    def plot_meshes(self):
        # show hand
        import vtk
        from Hand import vtk_add_from_hand, vtk_render
        renderer = vtk.vtkRenderer()
        # print(len(self.meshes))
        vtk_add_from_hand(self.meshes, renderer, 1.0, use_torch=True)
        vtk_render(renderer, axes=True)


