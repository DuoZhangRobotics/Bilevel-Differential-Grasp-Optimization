import torch
import numpy as np


def get_variables(params: torch.tensor):
    beta = params[0, 0]
    phi = params[0, 1]
    theta = params[0, 2: 5]
    d = params[0, 5]
    return beta, phi, theta, d


def get_n(rotation_matrix, beta, phi) -> torch.tensor:
    return rotation_matrix @ (torch.stack([torch.cos(beta) * torch.cos(phi), torch.sin(beta) * torch.cos(phi),
                                           torch.sin(phi)]).reshape(3, 1).requires_grad_(True))


def obj_fun(input_tensor_tuple):
    if isinstance(input_tensor_tuple, tuple):
        params, rotation_matrix, v0, v2, centroid0, centroid1, gamma = input_tensor_tuple
    else:
        raise TypeError("The input is not a tuple of tensor.")
    beta, phi, theta, d = get_variables(params)
    n = get_n(rotation_matrix, beta, phi)
    v1 = v0 + theta
    # using the parameters from the optimizer
    objective = -gamma * torch.sum(torch.log(torch.matmul(v2, n) - d)) - gamma * torch.sum(
        torch.log(d - torch.matmul(v1, n)))
    objective += torch.norm(centroid0 + theta - centroid1)
    return objective


def exp_reducer(x):
    return x.exp().sum()


input = torch.tensor([0, 0, 0, 1, 1, 1, 2, 2, 2], dtype=torch.double).requires_grad_(True)
out = exp_reducer(input)
j = torch.autograd.grad(out, input, retain_graph=True)
print(j)
print(input.grad)
out.backward()
print(input.grad)
