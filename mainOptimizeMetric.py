from ConvexHulls import ConvexHull
from HandObjective import HandObjective
from Optimizer import Optimizer
from Metric import Metric
from Hand import Hand
import numpy as np
import torch
import argparse


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--alpha", type=float, default=10)
    parser.add_argument("--gamma", type=float, default=0.0001)
    parser.add_argument("--soft_version", type=bool, default=False)
    parser.add_argument("--min_version", type=bool, default=False)
    parser.add_argument("--centroid_version", type=bool, default=False)
    parser.add_argument("--sample_count", type=int, default=100)
    args = parser.parse_args()

    gamma = torch.tensor(args.gamma, dtype=torch.double)
    alpha = torch.tensor(args.alpha, dtype=torch.double)
    soft_version = args.soft_version
    min_version = args.min_version
    centroid_version = args.centroid_version
    sample_count = args.sample_count
    print("GAMMA = ", gamma)
    print("ALPHA = ", alpha)
    print("SOFT VERSION: ", soft_version)
    print("MIN VERSION: ", min_version)
    print("CENTROID VERSION: ", centroid_version)
    print("SAMPLE COUNT: ", sample_count)

    path = 'hand/BarrettHand/'
    hand = Hand(path, scale=0.01, use_joint_limit=False, use_quat=False, use_eigen=False, use_contacts=False)
    if hand.use_eigen:
        params = torch.rand((1, hand.extrinsic_size + hand.eg_num))
    else:
        params = torch.rand((1, hand.extrinsic_size + hand.nr_dof()))
    p, lsp, t = hand.forward(params)

    # create object
    target = [ConvexHull(np.array([[-0.5,-0.5,-0.5],
                                   [-0.5, 0.5,-0.5],
                                   [ 0.5,-0.5,-0.5],
                                   [ 0.5, 0.5,-0.5],
                                   [-0.5, 0.5, 0.5],
                                   [ 0.5,-0.5, 0.5],
                                   [-0.5,-0.5, 0.5],
                                   [ 0.5, 0.5, 0.5]]) + 1.0),
            #   ConvexHull(np.array([[-0.5,-0.5,-0.5],
            #                        [-0.5, 0.5,-0.5],
            #                        [ 0.5,-0.5,-0.5],
            #                        [ 0.5, 0.5,-0.5],
            #                        [-0.5, 0.5, 0.5],
            #                        [ 0.5,-0.5, 0.5],
            #                        [-0.5,-0.5, 0.5],
            #                        [ 0.5, 0.5, 0.5]]) + 2.0)
                                   ]
    metric = Metric(target)
    # metric.setup_distance(hand)
    # metric.draw_samples(.1, 0.5)
    # for i in range(10):
    #     use_numpy = i%2
    #     if use_numpy:
    #         print("numpy: ",metric.compute_metric_numpy(hand))
    #     else:
    #         print("torch: ",metric.compute_metric_torch(hand))
    obj = HandObjective(hand, target, metric)


    def obj_func(param, hand_target):
        return obj.Q_metric_objective(param, gamma=gamma,
                                             alpha=alpha, 
                                             soft_version=soft_version, 
                                             min_version=min_version, 
                                             centroid_version=centroid_version)

    optimizer = Optimizer(obj_func, params=[obj.params, obj], method='Newton')
    optimizer.optimize(niters=100000, plot_interval=10000)
    optimizer.plot_meshes()
    optimizer.plot_history().savefig("history.png")
    print("GAMMA = ", gamma)
    print("ALPHA = ", alpha)
    print("SOFT VERSION: ", soft_version)
    print("MIN VERSION: ", min_version)
    print("CENTROID VERSION: ", centroid_version)
    print("SAMPLE COUNT: ", sample_count)

