import plotly.graph_objects as go
import torch
from convex_hulls import ConvexHulls
import numpy as np

pi = torch.acos(torch.zeros(1)).item() * 2


def plotting(points1, points2):
    hull1 = ConvexHulls(points1)
    hull2 = ConvexHulls(points2)

    # dist, x1, x2 = hull1.distance_between_convex_hulls(hull2)
    # optimal_dots = np.array([x1, x2])
    # print(dist, x1, x2)

    simplex2 = hull2.simplices
    fig = go.Figure()

    fig.add_trace(go.Mesh3d(x=points1[:, 0],
                            y=points1[:, 1],
                            z=points1[:, 2],
                            color='lightpink',
                            opacity=0.5,
                            legendgroup='Convex Hull1',
                            name='Convex Hull1',
                            alphahull=0,
                            showlegend=True))

    for i in range(simplex2.shape[0]):
        if i == 0:
            fig.add_trace(go.Mesh3d(x=points2[simplex2[i], 0],
                                    y=points2[simplex2[i], 1],
                                    z=points2[simplex2[i], 2],
                                    color='lightblue',
                                    opacity=0.5,
                                    legendgroup='Convex Hull2',
                                    name='Convex Hull2',
                                    showlegend=True))
        else:
            fig.add_trace(go.Mesh3d(x=points2[simplex2[i], 0],
                                    y=points2[simplex2[i], 1],
                                    z=points2[simplex2[i], 2],
                                    color='lightblue',
                                    opacity=0.5,
                                    legendgroup='Convex Hull2',
                                    name='Convex Hull2',
                                    showlegend=False))

    # fig.add_trace(go.Scatter3d(x=optimal_dots[:, 0],
    #                            y=optimal_dots[:, 1],
    #                            z=optimal_dots[:, 2],
    #                            marker=dict(size=4,
    #                                        color='tomato',
    #                                        symbol='diamond'),
    #                            line=dict(color='tomato',
    #                                      width=4),
    #                            legendgroup="Smallest Distance",
    #                            name="Smallest Distance"))
    fig.update_layout(scene_camera=dict(eye=dict(x=1.2, y=-1.6, z=1.0)),
                      margin=dict(t=0, r=10, l=10, b=10),
                      legend=dict(x=0.6, y=1),
                      scene=dict(xaxis=dict(range=[0, 4]),
                                 yaxis=dict(range=[0, 4]),
                                 zaxis=dict(range=[0, 4])),
                      scene_aspectmode='cube'
                      )
    fig.show(renderer='browser')
