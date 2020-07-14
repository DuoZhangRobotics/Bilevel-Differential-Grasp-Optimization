import plotly.graph_objects as go
import torch
import torch.nn as nn
from pprint import pprint
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from convex_hulls import ConvexHulls
import numpy as np

pi = torch.acos(torch.zeros(1)).item() * 2

if __name__ == '__main__':
    points1 = np.random.rand(30, 3)
    points2 = np.random.rand(30, 3) + 3
    hull1 = ConvexHulls(points1)
    hull2 = ConvexHulls(points2)

    dist, x1, x2 = hull1.distance_between_convex_hulls(hull2)
    optimal_dots = np.array([x1, x2])
    print(dist, x1, x2)
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter3D(points1[:, 0], points1[:, 1], points1[:, 2], c='blue')
    for simplex1 in hull1.simplices:
        ax.plot3D(points1[simplex1, 0], points1[simplex1, 1], points1[simplex1, 2], 'orange')

    ax.scatter3D(points2[:, 0], points2[:, 1], points2[:, 2], c='red')
    for simplex2 in hull2.simplices:
        ax.plot3D(points2[simplex2, 0], points2[simplex2, 1], points2[simplex2, 2], 'green')

    fig.show()

    simplex1 = hull1.simplices
    simplex2 = hull2.simplices
    fig = go.Figure()
    for i in range(simplex1.shape[0]):
        if i == 0:
            fig.add_trace(go.Mesh3d(x=points1[simplex1[i], 0],
                                    y=points1[simplex1[i], 1],
                                    z=points1[simplex1[i], 2],
                                    color='lightpink',
                                    opacity=0.5,
                                    legendgroup='Convex Hull1',
                                    name='Convex Hull1',
                                    showlegend=True))
        else:
            fig.add_trace(go.Mesh3d(x=points1[simplex1[i], 0],
                                    y=points1[simplex1[i], 1],
                                    z=points1[simplex1[i], 2],
                                    color='lightpink',
                                    opacity=0.5,
                                    legendgroup='Convex Hull1',
                                    name='Convex Hull1',
                                    showlegend=False))

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

    fig.add_trace(go.Scatter3d(x=optimal_dots[:, 0],
                               y=optimal_dots[:, 1],
                               z=optimal_dots[:, 2],
                               marker=dict(size=4,
                                           color='tomato',
                                           symbol='diamond'),
                               line=dict(color='tomato',
                                         width=4),
                               legendgroup="Smallest Distance",
                               name="Smallest Distance"))
    fig.update_layout(scene_camera=dict(eye=dict(x=1.2, y=-1.6, z=1.0)),
                      margin=dict(t=0, r=10, l=10, b=10),
                      legend=dict(x=0.6, y=1))
    fig.show(renderer='browser')
