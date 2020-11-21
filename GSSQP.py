import torch, cvxpy
import numpy as np
from HandTarget import HandTarget


class DirectionSolver:
    def __init__(self, hand_target:HandTarget, epsilon=torch.tensor(0.1, dtype=torch.double), alpha=torch.tensor(1, dtype=torch.double)):
        self.hand_target = hand_target
        self.epsilon = epsilon
        self.alpha = alpha

    def solve(self):
        pass
