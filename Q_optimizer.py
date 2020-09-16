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

    def optimize(self):
        v_num = sum([len(t.vertices) for t in self.hand_target.target])
        v = np.array([np.array(t.points[t.vertices]) for t in self.hand_target.target]).reshape((-1, 3))
        Q = cp.Variable(1)
        f = cp.Variable((v_num, 3))
        n, d, _ = self.hand_target.get_n_d(self.hand_target.params, self.hand_target.hand.palm, 0)
        n = n.detach().numpy()
        d = d.detach().numpy().reshape((-1, 1))

        constraints = [Q <= cp.sum(self.sampled_directions @ f.T, axis=1),
                       n @ f.T <= self.gamma * np.abs(d - n @ v.T),
                       # cp.sum_squares(n @ f.T) * self.mu >=
                       ]
        prob = cp.Problem(cp.Maximize(Q), constraints)
        prob.solve()

        return Q.value, f.value


if __name__ == '__main__':
    hand_target, optimizer = load_optimizer()
    optimizer.plot_meshes()
    sampled_directions = np.array(Directions(res=2, dim=3).dirs)
    qoptimizer = QOptimizer(hand_target, sampled_directions)
    qoptimizer.optimize()
