import torch
from convex_hulls import ConvexHulls


# define pi in torch
pi = torch.acos(torch.zeros(1)).item() * 2
data_type = torch.double


class ConvexHullSettings(object):
    def __init__(self, hull0: ConvexHulls, target: ConvexHulls):
        self.theta = torch.tensor([0, 0, 0], dtype=data_type).requires_grad_(True)
        self.hull0 = hull0
        self.hull2 = target
        self.distance, self.closest_pos0, self.closest_pos1 = hull0.distance_between_convex_hulls(target)
        self.centroid0, self.centroid2 = self._initialize_centroids()
        self.beta, self.phi = self._initialize_phi_beta()
        self.rotation_matrix = torch.eye(3, dtype=data_type)
        self.chart_reset()
        self.d = self._initialize_d().requires_grad_(True)

    def chart_reset(self, show_info=False):
        # reset beta and phi and n
        info = ''
        n = self.get_n()
        if show_info:
            info += (' n0=%s' % str(n.detach().numpy().T))
            info += (' beta=%s, phi=%s' % (str(self.beta.detach().numpy()), str(self.phi.detach().numpy())))
        rot_y_angle = -self.phi.detach()
        rot_z_angle = self.beta.detach()
        # print('rot_y_angle = %s, rot_z_angle = %s'%(rot_y_angle.detach().numpy(),rot_z_angle.detach().numpy()))
        rot_y = torch.tensor([[torch.cos(rot_y_angle), 0, torch.sin(rot_y_angle)],
                              [0, 1, 0],
                              [-torch.sin(rot_y_angle), 0, torch.cos(rot_y_angle)]
                              ], dtype=data_type)
        rot_z = torch.tensor([[torch.cos(rot_z_angle), -torch.sin(rot_z_angle), 0],
                              [torch.sin(rot_z_angle), torch.cos(rot_z_angle), 0],
                              [0, 0, 1]
                              ], dtype=data_type)
        self.beta.data.zero_()
        self.phi.data.zero_()
        self.rotation_matrix @= rot_z @ rot_y
        n = self.get_n()
        if show_info:
            info += (' beta=%s, phi=%s' % (str(self.beta.detach().numpy()), str(self.phi.detach().numpy())))
            info += (' n1=%s' % str(n.detach().numpy().T))
        return info

    def _initialize_phi_beta(self):
        closest_vec = torch.tensor(self.closest_pos1 - self.closest_pos0, dtype=data_type)
        closest_vec /= torch.norm(closest_vec)
        sin_phi = closest_vec[2]
        phi = torch.asin(sin_phi)
        beta = torch.atan2(closest_vec[1], closest_vec[0])
        return phi.detach().clone().requires_grad_(True), beta.detach().clone().requires_grad_(True)

    def _initialize_centroids(self):
        centroid0 = self._get_centroid(self.hull0)
        centroid2 = self._get_centroid(self.hull2)
        return centroid0, centroid2

    @staticmethod
    def _get_centroid(hull):
        return torch.tensor([torch.mean(torch.tensor(hull.points[hull.vertices, 0])),
                             torch.mean(torch.tensor(hull.points[hull.vertices, 1])),
                             torch.mean(torch.tensor(hull.points[hull.vertices, 2]))
                             ], dtype=data_type)

    def get_hulls(self):
        return self.hull0, self.hull2

    def get_vertices(self):
        points0 = torch.tensor(self.hull0.points, dtype=data_type)
        points2 = torch.tensor(self.hull2.points, dtype=data_type)
        v0 = points0[self.hull0.vertices, :]
        v2 = points2[self.hull2.vertices, :]
        return v0, v2

    def _initialize_d(self):
        v0, v2 = self.get_vertices()
        n = self.get_n()
        lower_bound = torch.min(v2.mm(n))
        upper_bound = torch.max(v0.mm(n))
        d = torch.mean(torch.tensor([lower_bound, upper_bound]))
        return d.detach().clone()

    def get_params(self):
        beta = self.beta.detach().clone().requires_grad_(True)
        phi = self.phi.detach().clone().requires_grad_(True)
        theta = self.theta.detach().clone().requires_grad_(True)
        d = self.d.detach().clone().requires_grad_(True)
        params = torch.cat([beta.reshape((1, 1)), phi.reshape(1, 1), theta.reshape(1, 3), d.reshape((1, 1))],
                           1).requires_grad_(True)
        v0, v2 = self.get_vertices()
        return [params, self.rotation_matrix, v0, v2, self.centroid0, self.centroid2]

    def get_n(self) -> torch.tensor:
        return self.rotation_matrix @ (
            torch.stack([torch.cos(self.beta) * torch.cos(self.phi),
                         torch.sin(self.beta) * torch.cos(self.phi),
                         torch.sin(self.phi)]).reshape(3, 1).requires_grad_(True))

    #
    def reset_parameters(self, params, chart_reset=True):
        params = params.reshape((-1))
        self.beta = params[0]
        self.phi = params[1]
        self.theta = params[2: 5]
        self.d = params[5]
        if chart_reset:
            self.chart_reset()
