import plotly.graph_objects as go
import torch
import torch.nn as nn
from pprint import pprint
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from convex_hulls import Convex_Hull
import numpy as np
pi = torch.acos(torch.zeros(1)).item() *2


if __name__ == '__main__':
    points1 = np.random.rand(30, 3)
    points2 = np.random.rand(30, 3) + 3
    hull1 = Convex_Hull(points1)
    hull2 = Convex_Hull(points2)
    
    dist, x1, x2 = hull1.distance_between_convex_hulls(hull2)
    print(dist, x1, x2)
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter3D(points1[:,0], points1[:,1], points1[:,2], c='blue')
    for simplex1 in hull1.simplices:
        ax.plot3D(points1[simplex1, 0], points1[simplex1, 1], points1[simplex1, 2],  'orange')

    ax.scatter3D(points2[:,0], points2[:,1], points2[:,2], c='red')
    for simplex2 in hull2.simplices:
        ax.plot3D(points2[simplex2, 0], points2[simplex2, 1], points2[simplex2, 2] ,'green')
    
    fig.show()