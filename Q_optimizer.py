import _pickle as pickle
import torch
import cvxpy as cp
import numpy as np
from HandTarget import HandTarget
from Optimizer import Optimizer
from Directions import Directions


data_type = torch.double


def obj_func(param, hand_target, gamma):
    return hand_target.hand_target_objective(param, gamma)


def load_optimizer(index: int = None):
    if index is None:
        with open(f'./log/HandTarget_final.pkl', 'rb') as input:
            hand_target: HandTarget = pickle.load(input)

        with open(f'./log/Optimizer_final.pkl', 'rb') as input:
            optimizer: Optimizer = pickle.load(input)

    else:
        with open(f'./log/HandTarget_{index}.pkl', 'rb') as input:
            hand_target: HandTarget = pickle.load(input)

        with open(f'./log/Optimizer_{index}.pkl', 'rb') as input:
            optimizer: Optimizer = pickle.load(input)

    return hand_target, optimizer


class QOptimizer(object):
    def __init__(self, hand_target, sampled_directions, gamma=0.001, mu=0.1):
        self.hand_target = hand_target
        self.sampled_directions = sampled_directions
        self.gamma = gamma
        self.mu = mu
        self.n, self.d = self.get_n_d()
        self.lamb = self.get_lambda()

    def optimize(self):
        Q = cp.Variable(1)
        f = cp.Variable((self.n.shape[0], 3))
        constraints = [Q <= cp.min(cp.sum(self.sampled_directions @ f.T, axis=1))]
        for i in range(self.n.shape[0]):
            n = self.n[i, :]
            constraints.append(n @ f[i, :].T <= self.lamb)
            constraints.append(cp.norm(f[i, :] - n @ f[i, :].T * n) <= self.mu * n @ f[i, :].T)

        prob = cp.Problem(cp.Maximize(Q), constraints)
        prob.solve()
        print(f"Q.value = {Q.value}, f.value = {f.value}")
        return Q.value, f.value

    def get_lambda(self):
        p, _ = self.hand_target.hand.forward(self.hand_target.params[:, :self.hand_target.front])
        closeness, _, _ = self.hand_target.closeness(self.hand_target.hand.palm, self.hand_target.target,
                                                     self.hand_target.params, p, 0, 0)
        closeness = closeness.detach().numpy()
        return self.gamma / closeness

    def get_n_d(self):
        n, d, _ = self.hand_target.get_n_d(self.hand_target.params, self.hand_target.hand.palm, 0)
        n = n.detach().numpy()
        d = d.detach().numpy().reshape((-1, 1))
        return n, d


if __name__ == '__main__':
    hand_target, optimizer = load_optimizer()
    # optimizer.plot_meshes()
    sampled_directions = np.array(Directions(res=2, dim=3).dirs)
    # sampled_directions = np.array([[1., 1, 1]])
    qoptimizer = QOptimizer(hand_target, sampled_directions)
    qoptimizer.optimize()
