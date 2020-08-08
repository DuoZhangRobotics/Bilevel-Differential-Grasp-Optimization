import matplotlib.pyplot as plt
from Optimizer import Optimizer
from copy import deepcopy


class Plot(object):

    def __init__(self, biOptimizer: Optimizer):
        self.biOptimizer = biOptimizer
        self.objectives = biOptimizer.objectives
        self.ch_settings = biOptimizer.ch_settings

    def plot_obj(self):
        fig = plt.figure()
        plt.plot(self.objectives, marker='d')
        plt.xlabel("iterations")
        plt.ylabel("value of objective function")
        fig.show()
        fig.savefig('Objective_function_values.png')

    def plot_convex_hulls(self):
        hull, hull2 = self.ch_settings.get_hulls()
        points1 = deepcopy(hull.points)
        points2 = deepcopy(hull2.points)

        theta = self.ch_settings.theta.detach().numpy()
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        ax.scatter3D(points1[:, 0], points1[:, 1], points1[:, 2], c='blue')
        for simplex1 in hull.simplices:
            ax.plot3D(points1[simplex1, 0], points1[simplex1, 1], points1[simplex1, 2], 'orange')

        points1 += theta
        ax.scatter3D(points1[:, 0], points1[:, 1], points1[:, 2], c='orchid')
        for simplex1 in hull.simplices:
            ax.plot3D(points1[simplex1, 0], points1[simplex1, 1], points1[simplex1, 2], 'dodgerblue')

        ax.scatter3D(points2[:, 0], points2[:, 1], points2[:, 2], c='tomato')
        for simplex2 in hull2.simplices:
            ax.plot3D(points2[simplex2, 0], points2[simplex2, 1], points2[simplex2, 2], 'lightblue')
        fig.show()
        fig.savefig('Convex_hulls.png')

    def plot_high_quality(self):
        import plotly.graph_objects as go
        hull, hull2 = self.ch_settings.get_hulls()
        theta = self.ch_settings.theta.detach().numpy()
        fig = go.Figure()
        points1 = deepcopy(hull.points)
        fig.add_trace(go.Mesh3d(x=points1[:, 0],
                                y=points1[:, 1],
                                z=points1[:, 2],
                                color='lightpink',
                                opacity=0.5,
                                legendgroup='Original Cube',
                                name='Original Cube',
                                alphahull=0,
                                showlegend=True))
        points2 = deepcopy(hull2.points)
        simplex2 = hull2.simplices
        for i in range(simplex2.shape[0]):
            if i == 0:
                fig.add_trace(go.Mesh3d(x=points2[simplex2[i], 0],
                                        y=points2[simplex2[i], 1],
                                        z=points2[simplex2[i], 2],
                                        color='tomato',
                                        opacity=0.5,
                                        legendgroup='Convex Hull2',
                                        name='Convex Hull2',
                                        showlegend=True))
            else:
                fig.add_trace(go.Mesh3d(x=points2[simplex2[i], 0],
                                        y=points2[simplex2[i], 1],
                                        z=points2[simplex2[i], 2],
                                        color='tomato',
                                        opacity=0.5,
                                        legendgroup='Convex Hull2',
                                        name='Convex Hull2',
                                        showlegend=False))

        points1 += theta
        fig.add_trace(go.Mesh3d(x=points1[:, 0],
                                y=points1[:, 1],
                                z=points1[:, 2],
                                color='dodgerblue',
                                opacity=0.5,
                                legendgroup='Cube',
                                name='Cube',
                                alphahull=0,
                                showlegend=True))

        fig.update_layout(scene_camera=dict(eye=dict(x=1.2, y=-1.6, z=1.0)),
                          margin=dict(t=0, r=10, l=10, b=10),
                          legend=dict(x=0.6, y=1),
                          # scene=dict(xaxis=dict(range=[0, 4]),
                          #            yaxis=dict(range=[0, 4]),
                          #            zaxis=dict(range=[0, 4])),
                          scene_aspectmode='cube'
                          )

        fig.show(renderer='browser')
