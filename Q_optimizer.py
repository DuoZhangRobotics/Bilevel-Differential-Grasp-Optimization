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
    def __init__(self, hand_target, non_linear_optimizer, sampled_directions, gamma=0.001, mu=0.1):
        self.hand_target = hand_target
        self.non_linear_optimizer = non_linear_optimizer
        self.sampled_directions = sampled_directions
        self.gamma = gamma
        self.mu = mu

    def optimize(self, tol=1e-5):
        v_num = sum([len(t.vertices) for t in self.hand_target.target])
        v = np.array([np.mean(np.array(t.points[t.vertices]), axis=0) for t in self.hand_target.target]).reshape((-1, 3))
        Q = cp.Variable(1)
        f = cp.Variable((1, 3))
        n, d, _ = self._get_n_d(self.hand_target.hand.palm, 0)
        n = np.array(n)
        d = np.array(d).reshape((-1, 1))

        constraints = [Q <= cp.sum(self.sampled_directions @ f.T, axis=1),
                       n @ f.T <= self.gamma * np.abs(d - n @ v.T),
                       # cp.sum_squares(n @ f.T) * self.mu >=
                       ]
        prob = cp.Problem(cp.Maximize(Q), constraints)
        prob.solve()
        print(f'The optimal value is: {prob.value}')
        print(f'The optimal value of Q is : {Q.value}')
        print(f'The optimal value of f is: {f.value}')

    def _get_n_d(self, root, idx):
        n = []
        d = []
        for i, t in enumerate(self.hand_target.target):
            beta = self.hand_target.params[0, self.hand_target.front + 3 * (idx * len(self.hand_target.target) + i) + 0]
            phi = self.hand_target.params[0, self.hand_target.front + 3 * (idx * len(self.hand_target.target) + i) + 1]
            rotation_matrix = self.hand_target.rotation_matrix[idx * len(self.hand_target.target) + i, :, :]
            n.append(self.hand_target.get_n(rotation_matrix, beta, phi).detach().numpy().reshape((-1)))
            d.append(self.hand_target.params[0, self.hand_target.front + 3 * (idx * len(self.hand_target.target) + i) + 2])
        idx += 1
        for child in root.children:
            tmp_n, tmp_d, idx = self._get_n_d(child, idx)
            n.extend(tmp_n)
            d.extend(tmp_d)
        return n, d, idx


if __name__ == '__main__':
    hand_target, optimizer = load_optimizer()
    sampled_directions = np.array(Directions(res=2, dim=3).dirs)
    qoptimizer = QOptimizer(hand_target, optimizer, sampled_directions)
    qoptimizer.optimize()