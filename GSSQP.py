import torch
import cvxpy as cp
import numpy as np
from HandTarget import HandTarget
from Directions import Directions
from ConvexHulls import ConvexHull
from Hand import Hand


class DirectionSolver:
    def __init__(self, hand_target: HandTarget, epsilon=torch.tensor(0.001, dtype=torch.double), alpha=torch.tensor(1, dtype=torch.double)):
        self.hand_target = hand_target
        self.epsilon = epsilon
        self.alpha = alpha

    def solve(self):
        params = self.hand_target.params.detach().clone().requires_grad_(True)
        D = cp.Variable(params[:, :self.hand_target.front].shape)
        Phi = cp.Variable(1)
        constraints = []
        p, t = self.hand_target.hand.forward(params[:, :self.hand_target.front])
        log_barrier, _, _ = self.hand_target.get_log_barrier(self.hand_target.hand.palm, self.hand_target.target, params, p, 0, 0)
        g_R = self.hand_target.linear_piece_enumeration(p, self.epsilon, self.alpha)
        for i in range(g_R.shape[1]):
            g_DR = torch.autograd.grad(g_R[0, i], params, retain_graph=True, create_graph=True)[0].detach().numpy()
            constraints.append(Phi >= g_DR[:, :self.hand_target.front] @ D.T)
        dlog_barrier = torch.autograd.grad(log_barrier, params)[0].detach().numpy()
        prob = cp.Problem(cp.Minimize(dlog_barrier[:, :self.hand_target.front] @ D.T + Phi), constraints=constraints)
        prob.solve()
        return torch.tensor(D.value, dtype=torch.double)


class GSSQP:
    def __init__(self):
        pass

    def optimize(self):
        pass


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
    count = 100
    sampled_directions = np.array(Directions(2, 3).dirs)
    # target = [ConvexHull(2 * np.random.rand(4, 3) + 2.)]
    target = [ConvexHull(np.array([[-1.0, -1.0, -1.0],
                                   [-1.0, 1.0, -1.0],
                                   [1.0, -1.0, -1.0],
                                   [1.0, 1.0, -1.0],
                                   [-1.0, 1.0, 1.0],
                                   [1.0, -1.0, 1.0],
                                   [-1.0, -1.0, 1.0],
                                   [1.0, 1.0, 1.0]]) + 2.0)]
    hand_target = HandTarget(hand, target, sampled_directions)
    direction_solver = DirectionSolver(hand_target=hand_target)
    direction_solver.solve()
