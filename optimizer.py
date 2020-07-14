import torch
from convex_hulls import ConvexHulls

pi = torch.acos(torch.zeros(1)).item() * 2


class BilevelOptimizer(object):
    def __init__(self, hull0: Convex_Hull, target: Convex_Hull):
        self.beta = torch.Tensor([0.0]).requires_grad_(True)
        self.phi = torch.Tensor([0.0]).requires_grad_(True)
        self.theta = torch.Tensor([0.0, 0.0, 0.0]).requires_grad_(True)
        self.hull0 = hull0
        self.hull1 = target
        self.distance, self.closest_pos0, self.closest_pos1 = hull0.distance_between_convex_hulls(target)
        self.d = torch.Tensor([self.distance / 10]).requires_grad_(True)
        self.n = torch.cat([torch.cos(self.beta) * torch.cos(self.phi), torch.sin(self.beta) * torch.cos(self.phi),
                            torch.sin(self.phi)]).reshape(1, 3).requires_grad_(True)

    def objective(self):
        # partial derivative of n with respect to beta
        pn_pbeta = torch.Tensor(
            [-torch.sin(self.beta) * torch.cos(self.phi), torch.cos(self.beta) * torch.cos(self.phi),
             torch.sin(self.phi)]).reshape(3, 1)
        # partial derivative of n with respect to phi
        pn_pphi = torch.Tensor(
            [-torch.cos(self.beta) * torch.sin(self.phi), -torch.sin(self.beta) * torch.sin(self.phi),
             torch.cos(self.phi)]).reshape(3, 1)
        # all the points contained in or on convex hulls
        points0 = torch.Tensor(self.hull0.points)
        #     ch1 = ConvexHull(points0 + theta)
        #     points1 = torch.Tensor(ch1.points)
        points1 = points0 + self.theta
        points2 = torch.Tensor(self.hull1.points)
        # vertices of two convex hulls
        #     v1 = points1[ch1.vertices, :].T
        v1 = points1[self.hull0.vertices, :].T
        v2 = points2[self.hull1.vertices, :].T
        # value of objective function
        objective = torch.sum(torch.log(torch.matmul(self.n, v2) - self.d * torch.ones(1, v2.shape[1]))) + torch.sum(
            torch.log(self.d * torch.ones(1, v1.shape[1]) - torch.matmul(self.n, v1)))
        objective *= -1
        # The chain rule. Partial derivativ e of O with respect to n
        po_pn = torch.sum(v2 / (self.n.mm(v2) - self.d), dim=1) - torch.sum(v1 / (self.d - self.n.mm(v1)), dim=1)
        po_pn = po_pn.reshape(1, 3)
        # Partial derivative of O with respect to beta and phi
        po_pbeta = po_pn.mm(pn_pbeta)
        po_pphi = po_pn.mm(pn_pphi)
        # Partial derivative of O with respect to d
        po_pd = -torch.sum(1 / (self.n.mm(v2) - self.d)) + torch.sum(1 / (self.d - self.n.mm(v1)))
        # partial derivative of O with respect to v

        po_pv = self.n * torch.sum(1 / (self.n.mm(v2) - self.d)) - self.n * torch.sum(1 / (self.d - self.n.mm(v1)))
        po_ptheta = po_pv
        return {'objective': objective, 'po_pbeta': po_pbeta, 'po_pphi': po_pphi, 'po_pd': po_pd,
                "po_ptheta": po_ptheta}

    def obj(self):
        points0 = torch.Tensor(self.hull0.points)
        #     ch1 = ConvexHull(points0 + theta)
        #     points1 = torch.Tensor(ch1.points)
        points1 = points0 + self.theta
        points2 = torch.Tensor(self.hull1.points)
        # vertices of two convex hulls
        #     v1 = points1[ch1.vertices, :].T
        v1 = points1[self.hull0.vertices, :].T
        v2 = points2[self.hull1.vertices, :].T
        # value of objective function
        objective = torch.sum(torch.log(torch.matmul(self.n, v2) - self.d * torch.ones(1, v2.shape[1]))) + torch.sum(
            torch.log(self.d * torch.ones(1, v1.shape[1]) - torch.matmul(self.n, v1)))
        objective *= -1
        return objective