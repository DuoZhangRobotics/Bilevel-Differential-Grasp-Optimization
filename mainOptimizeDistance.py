from HandObjective import HandObjective
from ConvexHulls import ConvexHull
from Optimizer import Optimizer
from Hand import Hand
import numpy as np
import torch

data_type = torch.double

if __name__ == "__main__":
    path = 'hand/BarrettHand/'
    hand = Hand(path, scale=0.01, use_joint_limit=False, use_quat=False, use_eigen=False, use_contacts=False)
    if hand.use_eigen:
        dofs = np.zeros(hand.eg_num)
        params = torch.zeros((1, hand.extrinsic_size + hand.eg_num))
    else:
        dofs = np.zeros(hand.nr_dof())
        params = torch.zeros((1, hand.extrinsic_size + hand.nr_dof()))
    p, t = hand.forward(params)

    # create object
    target = [ConvexHull(np.array([[-1.0,-1.0,-1.0],
                                   [-1.0, 1.0,-1.0],
                                   [ 1.0,-1.0,-1.0],
                                   [ 1.0, 1.0,-1.0],
                                   [-1.0, 1.0, 1.0],
                                   [ 1.0,-1.0, 1.0],
                                   [-1.0,-1.0, 1.0],
                                   [ 1.0, 1.0, 1.0]]) + 1.5),
              ConvexHull(np.array([[-1.0,-1.0,-1.0],
                                   [-1.0, 1.0,-1.0],
                                   [ 1.0,-1.0,-1.0],
                                   [ 1.0, 1.0,-1.0],
                                   [-1.0, 1.0, 1.0],
                                   [ 1.0,-1.0, 1.0],
                                   [-1.0,-1.0, 1.0],
                                   [ 1.0, 1.0, 1.0]]) + 2.0)]
    obj = HandObjective(hand, target)
    gamma = torch.tensor(0.001, dtype=data_type)

    def obj_func(param, hand_target):
        return obj.hand_objective(param, gamma)

    optimizer = Optimizer(obj_func, params=[obj.params, obj], method='Newton')
    optimizer.optimize(niters=1000, plot_interval=20)
    optimizer.plot_meshes()
    optimizer.plot_history().savefig("history.png")
