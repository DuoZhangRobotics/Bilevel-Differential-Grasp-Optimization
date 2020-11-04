from ConvexHulls import ConvexHull
from HandTarget import HandTarget
from Q_optimizer import QOptimizer, load_optimizer, obj_func
from MeritFunction import MeritFunction
from SQPInterface import SQP
from Directions import Directions
from Hand import Hand
import torch
import numpy as np


def obj(x: torch.tensor):
    # return 0.5 * (torch.exp(x) + torch.exp(1/x))
    # return 1/3 * torch.pow(x, 3) - x
    # return torch.square(x)

    # beta = x[:, 3]
    # phi = x[:, 4]
    return torch.log(x)


def constraints(x: torch.tensor):
    c = torch.cat(((x-1) * (x-5), x - 4))
    return c


def derivative(x: torch.tensor):
    # return 0.5 * (torch.exp(x) - (1/torch.square(x)) * torch.exp(1/x))
    return torch.square(x) - 1


if __name__ == "__main__":
    x = torch.tensor(100, dtype=torch.double).reshape((1, 1)).requires_grad_(True)
    gamma = 1.
    obj = lambda x: gamma * torch.log(x)
    while gamma > 1e-5:
        print("Gamma = ", gamma)
        sqp_solver = SQP(obj, constraints)
        x = sqp_solver.solve(x)
        gamma *= 0.1
    print(x)

