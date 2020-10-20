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
    return 1/3 * torch.pow(x, 3) - x
    # return torch.square(x)


def constraints(x: torch.tensor):
    c = torch.cat((3 - x, x - 4))
    return c


def derivative(x: torch.tensor):
    # return 0.5 * (torch.exp(x) - (1/torch.square(x)) * torch.exp(1/x))
    return torch.square(x) - 1


if __name__ == "__main__":
    sqp_solver = SQP(obj, constraints)
    x0 = torch.tensor(100, dtype=torch.double).reshape((1, 1)).requires_grad_(True)
    x_optimal = sqp_solver.solve(x0)
    print(x_optimal)
    #
    # x0 = torch.tensor(3.69916941, dtype=torch.double).reshape((1, 1)).requires_grad_(True)
    # x_optimal = sqp_solver.solve(x0)
    # print(x_optimal)
    # x0 = torch.tensor(3.97302029, dtype=torch.double).reshape((1, 1)).requires_grad_(True)
    # x_optimal = sqp_solver.solve(x0)
    # print(x_optimal)
