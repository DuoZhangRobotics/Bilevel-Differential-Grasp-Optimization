from LineSearcher import LineSearcher
from typing import Callable
import numpy as np
import torch, scipy
import pickle

data_type = torch.double

class Optimizer(object):
    def __init__(self, obj_func: Callable, params, method='Newton'):
        self.func = obj_func
        self.params = params
        self.method = method
        self.line_searcher = LineSearcher(self.func, self.params)

    def optimize(self, niters: int = 100000, tol: float = 1e-10, tolg: float = 1e-5, plot_interval=50):
        scale = 0.9
        invscale = 2.
        s = 1.0
        self.objectives = []
        self.grad_norms = []
        self.meshes = [i.mesh() for i in self.params[1].target]
        self.meshes.append(self.params[1].hand.draw(scale_factor=1, show_to_screen=False, use_torch=True))

        line_searcher = LineSearcher(self.func, self.params)
        for i in range(niters):
            last_s = s

            # find search direction
            obj = self.func(*self.params)
            jacobian: torch.tensor = torch.autograd.grad(obj, self.params[0], retain_graph=True, create_graph=True)[0]
            if self.method == "Newton":
                hessian, L = self.hessian(jacobian)
                jacobian = jacobian.detach().numpy()
                direction = scipy.linalg.cho_solve(L, jacobian.T)
            else:
                jacobian = jacobian.detach().numpy()
                direction = jacobian.T

            # line search
            s, new_params, obj = line_searcher.line_search(grad=jacobian.T, direction=direction, obj=obj, tol=tol,
                                                           scale=scale, s=last_s)
            if s is None:
                print("Line-Search failed!")
                break

            # Riemann optimization: chart reset
            with torch.no_grad():
                self.params[1].reset_parameters(new_params[0], chart_reset=True)
            self.params[0] = self.params[1].params.requires_grad_(True)

            # adaptive scaling of search length
            if s == last_s:
                s *= invscale
            if i % plot_interval == 0:
                self.meshes.append(self.params[1].hand.draw(scale_factor=1, show_to_screen=False, use_torch=True))
                self.plot_meshes()
                import os
                if not os.path.exists('./log'):
                    os.mkdir('./log')
                    
                with open(rf'./log/Optimizer_{i}.pkl', 'wb') as output:
                    pickle.dump(self, output, pickle.HIGHEST_PROTOCOL)

                with open(rf'./log/HandTarget_{i}.pkl', 'wb') as output:
                    pickle.dump(self.params[1], output, pickle.HIGHEST_PROTOCOL)

            # record/print
            self.objectives.append(obj)
            self.grad_norms.append(np.max(np.abs(jacobian)))
            print(f"Iter{i:3d}: obj={self.objectives[-1]:3.6f} grad={self.grad_norms[-1]:3.6f} x={self.params[0][:, :3].detach().numpy()} s={s:3.6f}")

            # convergence check
            if self.grad_norms[-1] < tolg:
                print("Converged!")
                
                import os
                if not os.path.exists('./log'):
                    os.mkdir('./log')
                    
                with open(rf'./log/Optimizer_final.pkl', 'wb') as output:
                    pickle.dump(self, output, pickle.HIGHEST_PROTOCOL)

                with open(rf'./log/HandTarget_final.pkl', 'wb') as output:
                    pickle.dump(self.params[1], output, pickle.HIGHEST_PROTOCOL)
                break

    def hessian(self, jacobian):
        length = jacobian.size(1)
        hessian = torch.zeros((length, length), dtype=data_type)
        for i in range(length):
            # with torch.autograd.detect_anomaly():
            hessian[i, :] = torch.autograd.grad(jacobian[:, i], self.params[0], retain_graph=True)[0]
        try:
            # check if the hessian matrix is positive definite
            hessian = hessian.detach().numpy()
            return hessian, scipy.linalg.cho_factor(hessian)
        except np.linalg.LinAlgError:
            # if the hessian matrix is not positive definite, then make it positive definite.
            hessian = torch.tensor(self.make_positive_definite(hessian), dtype=data_type)
            hessian = hessian.detach().numpy()
            return hessian, scipy.linalg.cho_factor(hessian)

    def make_positive_definite(self, hessian, min_cond=0.00001):
        eigenvalues, eigenvectors = np.linalg.eig(hessian)
        l = np.max(np.abs(eigenvalues))
        for i in range(len(eigenvalues)):
            eigenvalues[i] = max(eigenvalues[i], l * min_cond)
        return eigenvectors @ np.diag(eigenvalues) @ np.transpose(eigenvectors)

    def grad_check(self):
        print(torch.autograd.gradcheck(self.func, self.params))

    def plot_meshes(self):
        # show hand
        import vtk
        from Hand import vtk_add_from_hand, vtk_render
        renderer = vtk.vtkRenderer()
        vtk_add_from_hand(self.meshes, renderer, 1.0, use_torch=True)
        vtk_render(renderer, axes=True)

    def plot_history(self):
        # show convergence
        import matplotlib.pyplot as plt
        fig = plt.figure()
        plt.plot(self.objectives, marker='d', label="Objective")
        plt.plot(self.grad_norms, marker='d', label="Gradient")
        plt.title("Convergence History")
        plt.legend(loc="upper right")
        plt.xlabel("#Iterations")
        plt.ylabel("Value")
        fig.show()
        return fig
