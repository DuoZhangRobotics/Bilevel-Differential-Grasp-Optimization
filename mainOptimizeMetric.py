from ConvexHulls import ConvexHull
from Metric import Metric
from Hand import Hand
import numpy as np
import torch

if __name__ == "__main__":
    path = 'hand/BarrettHand/'
    hand = Hand(path, scale=0.01, use_joint_limit=False, use_quat=False, use_eigen=False, use_contacts=False)
    if hand.use_eigen:
        params = torch.rand((1, hand.extrinsic_size + hand.eg_num))
    else:
        params = torch.rand((1, hand.extrinsic_size + hand.nr_dof()))
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
    metric = Metric(target)
    metric.setup_distance(hand)
    metric.draw_samples(.1, 0.5)
    for i in range(10):
        use_numpy = i%2
        if use_numpy:
            print("numpy: ",metric.compute_metric_numpy(hand))
        else:
            print("torch: ",metric.compute_metric_torch(hand))