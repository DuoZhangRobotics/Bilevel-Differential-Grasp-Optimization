from ConvexHulls import ConvexHull
from Metric import Metric

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
    Metric(target)