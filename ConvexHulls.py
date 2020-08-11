import torch,scipy,trimesh
import numpy as np
import cvxpy as cp
data_type = torch.double

class ConvexHull(object):
    def __init__(self, points: np.array):
        self.points = points
        self.hull = scipy.spatial.ConvexHull(self.points)
        self.A = self.hull.equations[:, :-1]
        self.b = -self.hull.equations[:, -1]
        self.simplices = self.hull.simplices
        self.vertices = self.hull.vertices
        self.centroid = self.mesh().centroid

    def distance_to(self, target_hull):
        # distance between two convex hulls
        A2, b2 = target_hull.A, target_hull.b

        x1 = cp.Variable(3)
        x2 = cp.Variable(3)

        objective = cp.Minimize(cp.sum_squares(x1 - x2))
        constraints = [self.A @ x1 <= self.b, A2 @ x2 <= b2]
        prob = cp.Problem(objective, constraints)
        min_dist = prob.solve(solver='SCS')
        return min_dist, x1.value, x2.value
    
    def surface_vertices(self):
        return self.points[self.vertices, :]

    def surface_indices(self):
        faces=[]
        for f in self.simplices:
            if f[0] in self.vertices and f[1] in self.vertices and f[2] in self.vertices:
                faces.append([self.vertices.tolist().index(f[0]),self.vertices.tolist().index(f[1]),self.vertices.tolist().index(f[2])])
        return np.array(faces)
    
    def mesh(self):
        return trimesh.Trimesh(vertices=self.surface_vertices(), faces=self.surface_indices())
