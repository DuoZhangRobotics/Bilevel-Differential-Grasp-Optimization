import torch
import numpy as np
from Hand import Link, Hand
import trimesh
from ConvexhullSettings import ConvexHullSettings
from HandTarget import HandTarget


def get_variables(params: torch.tensor):
    beta = params[0, 0]
    phi = params[0, 1]
    theta = params[0, 2: 5]
    d = params[0, 5]
    return beta, phi, theta, d


def get_n(rotation_matrix, beta, phi) -> torch.tensor:
    return rotation_matrix @ (torch.stack([torch.cos(beta) * torch.cos(phi), torch.sin(beta) * torch.cos(phi),
                                           torch.sin(phi)]).reshape(3, 1).requires_grad_(True))


def obj_fun(params, ch_settings: ConvexHullSettings, gamma=torch.tensor(0.001, dtype=torch.double)):
    beta, phi, theta, d = get_variables(params)
    v0, v2 = ch_settings.get_vertices()
    centroid0, centroid1 = ch_settings.get_centroids()
    n = get_n(ch_settings.rotation_matrix, beta, phi)
    v1 = v0 + theta
    # using the parameters from the optimizer
    objective = -gamma * torch.sum(torch.log(torch.matmul(v2, n) - d)) - gamma * torch.sum(
        torch.log(d - torch.matmul(v1, n)))
    objective += torch.norm(centroid0 + theta - centroid1)
    return objective


def hand_obj_fun(params, hand_target: HandTarget, gamma=torch.tensor(0.01, dtype=torch.double)):
    motion_params = params[:, :hand_target.front]
    # objective = torch.square(params).sum()
    p, t = hand_target.hand.forward(motion_params)
    norm, _ = get_norm(hand_target.hand.palm, hand_target.target, p, 0)
    print(f'Norm = {norm}')
    objective, _, _ = get_constrain(hand_target.hand.palm, hand_target.target, hand_target, params, p, 0, 0)
    print(f'Constraint = {objective}')
    objective = objective * gamma
    objective = objective + norm
    return objective


def get_norm(root: Link, target_link: trimesh.Trimesh, p, start=0):
    centroid1 = torch.tensor(target_link.centroid, dtype=torch.double)
    centroid0 = torch.mean(p[:, :, start: start + len(root.mesh.vertices)], dim=2)
    norm = torch.norm(centroid0 - centroid1)
    # if start == 0:
    #     #     norm = norm * 100
    start += len(root.mesh.vertices)
    if root.children:
        for child in root.children:
            tmp_norm, start = get_norm(child, target_link, p, start)
            norm = norm + tmp_norm
    return norm, start


def get_constrain(root: Link, target_link: trimesh.Trimesh, hand_target: HandTarget, params, p, start, idx):
    # print(f"mesh = {root.mesh}")
    v1 = p[:, :, start: start + len(root.mesh.vertices)].view(3, -1).T
    v2 = torch.tensor(target_link.vertices, dtype=torch.double)
    rotation_matrix = hand_target.get_rotation_matrix(idx)
    beta, phi, d = get_beta_phi_d(params, hand_target, idx)
    n = get_n(rotation_matrix, beta, phi)
    lower_bound2 = torch.min(v2 @ n)
    upper_bound2 = torch.max(v2 @ n)
    upper_bound = torch.max(v1 @ n)
    lower_bound = torch.min(v1 @ n)
    # print(f'idx, beta, phi, n, d = {idx} {beta} {phi} {n.T} {d}')
    # print(f'rotation matrix = \n{rotation_matrix}')

    # print(f'v1 @ n - d = {torch.sum(torch.log(v1 @ n - d))}')
    # print(f'v2 @ n - d = {torch.sum(torch.log(v2 @ n - d))}')
    # print(f'd - v1 @ n = {torch.sum(torch.log(d - v1 @ n))}')
    # print(f'd - v2 @ n = {torch.sum(torch.log(d - v2 @ n))}')
    if lower_bound > upper_bound2:
        objective = -torch.sum(torch.log(v1 @ n - d)) - torch.sum(torch.log(d - v2 @ n))
        # objective = -torch.sum(v1 @ n - d) - torch.sum(d - v2 @ n)

    elif lower_bound2 > upper_bound:
        objective = -torch.sum(torch.log(v2 @ n - d)) - torch.sum(torch.log(d - v1 @ n))
        # objective = -torch.sum(v2 @ n - d) - torch.sum(d - v1 @ n)

    else:
        objective = -torch.sum(torch.log(v2 @ n - d)) - torch.sum(torch.log(d - v1 @ n))
        # objective = -torch.sum(v2 @ n - d) - torch.sum(d - v1 @ n)

    start += len(root.mesh.vertices)
    idx += 1
    for child in root.children:
        tmp_objective, start, idx = get_constrain(child, target_link, hand_target, params, p, start, idx)
        objective = objective + tmp_objective
    return objective, start, idx


def get_beta_phi_d(params, hand_target, idx):
    beta = params[0, hand_target.front + 3 * idx]
    phi = params[0, hand_target.front + 3 * idx + 1]
    d = params[0, hand_target.front + 3 * idx + 2]

    return beta, phi, d


