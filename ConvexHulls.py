import torch, scipy, trimesh
import numpy as np
import cvxpy as cp

data_type = torch.double


class ConvexHull(object):
    def __init__(self, points: np.array):
        self.points = points
        self.reset()
        
    def reset(self):
        self.hull = scipy.spatial.ConvexHull(self.points)
        self.A = self.hull.equations[:, :-1]
        self.b = -self.hull.equations[:, -1]
        self.simplices = self.hull.simplices
        self.vertices = self.hull.vertices
        self.centroid = self.mesh().centroid

    def distance_to(self, target_hull):
        # distance between two convex hulls
        x1 = cp.Variable(3)

        if isinstance(target_hull, ConvexHull):
            x2 = cp.Variable(3)
            A2, b2 = target_hull.A, target_hull.b
            constraints = [self.A @ x1 <= self.b, A2 @ x2 <= b2]
        else:
            x2 = [target_hull[0], target_hull[1], target_hull[2]]
            constraints = [self.A @ x1 <= self.b]

        objective = cp.Minimize(cp.sum_squares(x1 - x2))
        prob = cp.Problem(objective, constraints)
        min_dist = prob.solve(solver=cp.MOSEK)

        if isinstance(target_hull, ConvexHull):
            return min_dist, x1.value, x2.value
        else:
            return min_dist, x1.value, x2

    def surface_vertices(self):
        return self.points[self.vertices, :]

    def surface_indices(self):
        faces = []
        for f in self.simplices:
            if f[0] in self.vertices and f[1] in self.vertices and f[2] in self.vertices:
                faces.append([self.vertices.tolist().index(f[0]), self.vertices.tolist().index(f[1]),
                              self.vertices.tolist().index(f[2])])
        return np.array(faces)

    def mesh(self):
        return trimesh.Trimesh(vertices=self.surface_vertices(), faces=self.surface_indices())

    def mass(self):
        m = self.mesh()
        return m.center_mass,m.mass

    def translate(self, dt):
        self.points += dt
        self.reset()

    def contain(self, x2, thres=1e-3):
        return self.distance_to(x2)[0] < thres
