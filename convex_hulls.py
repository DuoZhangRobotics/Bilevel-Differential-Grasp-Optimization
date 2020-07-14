import numpy as np
import cvxpy as cp
from scipy.spatial import ConvexHull


class Convex_Hull(object):
    def __init__(self, points_cloud1: np.array):
        self.points1 = points_cloud1
        self.hull = self.__generate_convex_hulls()
        self.A = self.hull.equations[:, :-1]
        self.b = -self.hull.equations[:, -1]
        self.simplices = self.hull.simplices
        self.vertices = self.hull.vertices

    def __generate_convex_hulls(self):
        return ConvexHull(self.points1)

    # distance between two convex hulls
    def distance_between_convex_hulls(self, target_hull):
        A2, b2 = target_hull.A, target_hull.b

        x1 = cp.Variable(3)
        x2 = cp.Variable(3)

        objective = cp.Minimize(cp.sum_squares(x1 - x2))
        constraints = [self.A @ x1 <= self.b, A2 @ x2 <= b2]
        prob = cp.Problem(objective, constraints)
        min_dist = prob.solve()
        return min_dist, x1.value, x2.value
