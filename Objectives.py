import torch
import numpy as np
from hand import Link, Hand
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
    p, t = hand_target.hand.forward(motion_params)
    norm = get_norm(hand_target.hand.palm, hand_target.target, p, 0)
    objective = gamma * get_constrain(hand_target.hand.palm, hand_target.target, hand_target, p, 0, 0)
    print(objective)
    print(norm)
    objective += norm
    return norm


def get_norm(root: Link, target_link: trimesh.Trimesh, p, start=0):
    centroid1 = torch.tensor(target_link.centroid, dtype=torch.double)
    centroid0 = torch.mean(p[:, :, start: start + len(root.mesh.vertices)], dim=2)
    print(centroid0, centroid1)
    norm = torch.norm(centroid0 - centroid1)
    start += len(root.mesh.vertices)
    if root.children:
        for child in root.children:
            norm += get_norm(child, target_link, p, start)
    return norm


def get_constrain(root: Link, target_link: trimesh.Trimesh, hand_target: HandTarget, p, start, idx):
    v1 = p[:, :, start: start + len(root.mesh.vertices)].view(3, -1).T
    v2 = torch.tensor(target_link.vertices, dtype=torch.double)
    rotation_matrix = hand_target.get_rotation_matrix(idx)
    beta, phi = hand_target.get_beta_phi(idx)
    d = hand_target.get_d(idx)
    n = get_n(rotation_matrix, beta, phi)
    objective = -torch.sum(torch.log(v1 @ n - d)) - torch.sum(torch.log(d - v2 @ n))

    start += len(root.mesh.vertices)
    idx += 1
    for child in root.children:
        objective += get_constrain(child, target_link, hand_target, p, start, idx)

    return objective
