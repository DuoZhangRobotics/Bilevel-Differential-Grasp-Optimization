from ConvexHulls import ConvexHull
from HandTarget import HandTarget
from Optimizer import Optimizer
from Hand import Hand
import torch,trimesh
import numpy as np

data_type = torch.double

if __name__ == "__main__":
    path = 'hand/ShadowHand/'
    hand = Hand(path, scale=0.01, use_joint_limit=False, use_quat=False, use_eigen=False, use_contacts=False)
    if hand.use_eigen:
        dofs = np.zeros(hand.eg_num)
        params = torch.zeros((1, hand.extrinsic_size + hand.eg_num))
    else:
        dofs = np.zeros(hand.nr_dof())
        params = torch.zeros((1, hand.extrinsic_size + hand.nr_dof()))
    p, t = hand.forward(params)
    
    #create object
    hull = ConvexHull(np.random.rand(30, 3) + np.array([0., -0.2, 1.3]))
    hand_target = HandTarget(hand, hull)
    gamma=torch.tensor(0.001, dtype=data_type)
    def obj_func(param, hand_target, gamma):
        return hand_target.hand_target_objective(param, gamma)
    optimizer = Optimizer(obj_func, params=[hand_target.params, hand_target, gamma], method='Newton', gamma=gamma)
    optimizer.optimize(niters=500, plot_interval=20)
    optimizer.plot_history()
    optimizer.plot_meshes()
