import numpy as np,math
from scipy.spatial import KDTree
from mpl_toolkits.mplot3d import Axes3D
from scipy.spatial.distance import pdist, squareform

reflection = np.array([[0., -1.], [-1., 0.]])

def mesh_area(triangle_list):
    N = np.cross(triangle_list[:,1] - triangle_list[:,0], triangle_list[:,2] - triangle_list[:,0], axis = 1)
    N_norm = np.sqrt(np.sum(N ** 2, axis = 1))
    normals = N / np.expand_dims(N_norm, axis = 1)
    N_norm *= .5
    return N_norm, normals

def triangle_point_picking(triangle_list):
    # Compute uniform distribution over [0, 1]x[0, 1] lower triangle
    X = np.random.random((triangle_list.shape[0], 2))
    t = np.sum(X, axis = 1) > 1
    X[t] = np.dot(X[t], reflection) + 1.

    # Map the [0, 1]x[0, 1] lower triangle to the actual triangles
    ret = np.einsum('ijk,ij->ik', triangle_list[:,1:] - triangle_list[:,0,None], X) 
    ret += triangle_list[:,0]
    return ret

def uniform_sample_mesh(triangle_list, normal_list, triangle_area_list, sample_count):
    # Normalize the sum of area of each triangle to 1
    triangle_area = triangle_area_list / np.sum(triangle_area_list)
    '''
    For each sample
      * Pick a triangle with probability proportial to its surface area
      * pick a point on that triangle with uniform probability
    '''
    triangle_id_list = np.random.choice(triangle_list.shape[0], size = sample_count, p = triangle_area)
    return triangle_point_picking(triangle_list[triangle_id_list]), normal_list[triangle_id_list]

def blue_noise_sample_elimination(point_list, normal_list, mesh_surface_area, sample_count):
    # Parameters
    alpha = 8
    rmax = np.sqrt(mesh_surface_area / ((2 * sample_count) * np.sqrt(3.))) 

    # Compute a KD-tree of the input point list
    kdtree = KDTree(point_list)

    # Compute the weight for each sample
    D = np.minimum(squareform(pdist(point_list)), 2 * rmax)
    D = (1. - (D / (2 * rmax))) ** alpha

    W = np.zeros(point_list.shape[0])
    for i in range(point_list.shape[0]):
        W[i] = sum(D[i, j] for j in kdtree.query_ball_point(point_list[i], 2 * rmax) if i != j)

    # Pick the samples we need
    heap = sorted((w, i) for i, w in enumerate(W))

    id_set = set(range(point_list.shape[0]))
    while len(id_set) > sample_count:
        # Pick the sample with the highest weight
        w, i = heap.pop()
        id_set.remove(i)

        neighbor_set = set(kdtree.query_ball_point(point_list[i], 2 * rmax))
        neighbor_set.remove(i)
        heap = [(w - D[i, j], j) if j in neighbor_set else (w, j) for w, j in heap]                
        heap.sort()

    # Job done
    return point_list[sorted(id_set)], normal_list[sorted(id_set)]

def sample_convex_hulls(targets, counts, multiple=10):
    triangles = []
    for t in targets:
        m = t.mesh()
        triangles.append(np.array([[[m.vertices[fi][0],m.vertices[fi][1],m.vertices[fi][2]] for fi in f] for f in m.faces])) 
    
    tri_areas = []
    tri_normals =[]
    for triangle in triangles:
        area, normal = mesh_area(triangle)
        tri_areas.append(area)
        tri_normals.append(normal)
    
    tri_area_sums = [np.sum(tri_area) for tri_area in tri_areas]
    tri_area_sums_all = np.sum(tri_area_sums)
    
    point_normals = [uniform_sample_mesh(triangle, tri_normal, tri_area, int(math.ceil(multiple * counts * tri_area_sum / tri_area_sums_all))) 
                    for triangle, tri_normal, tri_area, tri_area_sum in zip(triangles, tri_normals, tri_areas, tri_area_sums)]
    for ti in range(len(targets)):
        pids = []
        points = point_normals[ti][0]
        for pid in range(points.shape[0]):
            contain = False
            for tj in range(len(targets)):
                if ti!=tj and targets[tj].contain(points[pid]):
                    contain=True
                    break
            if not contain:
                pids.append(pid)
        point_normals[ti] = [point_normals[ti][0][pids], point_normals[ti][1][pids]]

    points = np.vstack([pn[0] for pn in point_normals])
    normals = np.vstack([pn[1] for pn in point_normals])
    return blue_noise_sample_elimination(points, normals, tri_area_sums_all, counts)
