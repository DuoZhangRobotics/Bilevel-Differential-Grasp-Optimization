from ConvexHulls import ConvexHull
from HandTarget import HandTarget
from Q_optimizer import QOptimizer
from SQPInterface import SQP
from Hand import Hand
import torch
import numpy as np


def Alg2(hand_target: HandTarget, sampled_direction, x0, u0, niters=100000):
    # TODO: Fill out the frame
    x = x0
    u = u0
    p, _ = hand_target.hand.forward(x[:, :hand_target.front])
    for i in range(niters):
        qoptimizer = QOptimizer(hand_target, sampled_direction)
        Q, F = qoptimizer.optimize()
        function = hand_target.get_log_barrier(hand_target.hand.palm, hand_target.target, x, p, 0, 0)
        constraints = hand_target.friction_cone_constraint(F, 0.001)
        sqp_solver = SQP(function, constraints, u=u)
        x, u = sqp_solver.solve()
        hand_target.reset_parameters(x, True)
        p, _ = hand_target.hand.forward(x[:, :hand_target.front])


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
    target = [ConvexHull(np.random.rand(4, 3) + np.array([0.0, 0.0, 1.0]))]
    hand_target = HandTarget(hand, target)
    gamma = torch.tensor(0.001, dtype=torch.double)

    qoptimizer = QOptimizer(hand_target, np.array([0.0, 0.0, 1.0]))
    Q, F = qoptimizer.optimize()
