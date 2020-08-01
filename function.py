import torch
import numpy as np
from hand import Link, Hand
import trimesh
from ConvexhullSettings import ConvexHullSettings


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


def hand_obj_fun(params, hand: Hand, target_mesh: trimesh.Trimesh, gamma=torch.tensor(0.01, dtype=torch.double)):
    # if hand.use_eigen:
    #     extrinsic_mat = params[0, 0:7]
    #     dofs = params[0, 7:]
    # else:
    #     extrinsic_mat = params[0, 0:6]
    #     dofs = params[0, 6:]
    # hand.forward_kinematics(extrinsic_mat, dofs)
    hand.forward(params)
    norm = get_norm(hand.palm, target_mesh)
    objective = gamma * get_constrain(hand.palm, target_mesh)
    objective += norm
    return norm


def get_norm(root: Link, target_link: trimesh.Trimesh):
    centroid1 = target_link.centroid
    norm = torch.norm(root.centroid - centroid1)
    for child in root.children:
        norm += torch.norm(child.center - centroid1)
        if len(child.children) != 0:
            norm += get_norm(child, target_link)
    return norm


def get_constrain(root: Link, target_link: trimesh.Trimesh):
    v1, v2 = torch.tensor(root.mesh.vertices, dtype=torch.double), torch.tensor(target_link.vertices,
                                                                                dtype=torch.double)
    n = get_n(root.rotation_matrix, root.beta, root.phi)
    objective = -torch.sum(torch.log(v2 @ n - root.d)) - torch.sum(torch.log(root.d - v1 @ n))
    for child in root.children:
        objective += get_constrain(child, target_link)

    return objective
