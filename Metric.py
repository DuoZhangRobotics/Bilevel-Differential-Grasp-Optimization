from ConvexHulls import ConvexHull
from Hand import Hand
import scipy, trimesh
import numpy as np


class Metric(object):
    def __init__(self, targets, count=100, friCoef=.7, res=3, M=None):
        self.targets = targets
        from PoissonDiskSampling import sample_convex_hulls
        self.points, self.normals = sample_convex_hulls(self.targets, count)
        self.mu = friCoef
        
        #sample directions
        from Directions import Directions
        self.dirs = Directions(res=res, dim=6).dirs
        self.dirs = np.array(self.dirs,dtype=np.float64)
        
        #compute gij
        f2w = np.zeros((self.points.shape[0],6,3), dtype=np.float64)
        f2w[:,0,0] = f2w[:,1,1] = f2w[:,2,2] = 1
        f2w[:,2,1] = self.points[:,0]
        f2w[:,1,2] =-self.points[:,0]
        f2w[:,0,2] = self.points[:,1]
        f2w[:,2,0] =-self.points[:,1]
        f2w[:,1,0] = self.points[:,2]
        f2w[:,0,1] =-self.points[:,2]
        if M is not None:
            f2w = np.einsum("ij,kjl->kil", M, f2w)
        df2w = np.einsum("ij,kjl->kil", self.dirs, f2w)
        w_perp = np.einsum("ijk,ik->ij", df2w, self.normals)
        w_para = np.linalg.norm(df2w - np.einsum("ij,il->ijl", w_perp, self.normals),axis=2)
        
        case1 = w_perp + w_para**2 / w_perp
        case2 = np.maximum(0 , w_perp + self.mu * w_para)
        choice = self.mu * w_perp > w_para
        self.gij = np.where(choice, case1, case2)
        
    def draw_samples(self, scale, length):
        # show hand
        import vtk
        from Hand import vtk_render, vtk_add_from_hand
        renderer = vtk.vtkRenderer()
        vtk_add_from_hand([t.mesh() for t in self.targets], renderer, 1.0, use_torch=True)
        vtk_add_metric_samples(renderer, self, scale, length)
        vtk_render(renderer, axes=True)

    def compute_exp_dist(self, link):
        if isinstance(link, Hand):
            return self.compute_exp_dist(link.palm)
        else:
            trans = link.joint_transform_torch.numpy()[0,:,:]
            pss = trans[:3,:3].T @ (self.points - trans[:3,3]).T
            ret = [[link.hull.distance_to(pss[:,i]) for i in range(pss.shape[1])]]
            if link.children:
                for c in link.children:
                    ret += self.compute_exp_dist(c)
            return ret

    def compute_metric(self, hand, alpha=.1):
        dists = self.compute_exp_dist(hand)
        dists_value = np.array([[d[0] for d in distsL] for distsL in dists])
        print(dists_value.shape)
        exp_dists_max = np.amax(np.exp(dists_value * -alpha), axis=0)
        Q_value = self.gij * np.expand_dims(exp_dists_max, axis=1)
        Q_value = np.mean(np.amax(Q_value, axis=1), axis=0)
        return Q_value

def vtk_add_metric_samples(renderer, metric, scale, length):
    import vtk
    for i in range(metric.points.shape[0]):
        #point
        sphere = vtk.vtkSphereSource()
        sphere.SetCenter(metric.points[i,0],
                         metric.points[i,1],
                         metric.points[i,2])
        sphere.SetRadius(scale)
        sphere.SetThetaResolution(24)
        sphere.SetPhiResolution(24)
        sphere_mapper = vtk.vtkPolyDataMapper()
        sphere_mapper.SetInputConnection(sphere.GetOutputPort())
        sphere_actor = vtk.vtkActor()
        sphere_actor.SetMapper(sphere_mapper)
        sphere_actor.GetProperty().SetColor(255. / 255, .0 / 255, 255.0 / 255)
        renderer.AddActor(sphere_actor)

        # normal
        if length <= 0.:
            continue
        normal = vtk.vtkArrowSource()
        normal.SetTipResolution(100)
        normal.SetShaftResolution(100)
        # Generate a random start and end point
        startPoint=[metric.points[i,0],
                    metric.points[i,1],
                    metric.points[i,2]]
        endPoint=[metric.points[i,0]+metric.normals[i,0]*length,
                  metric.points[i,1]+metric.normals[i,1]*length,
                  metric.points[i,2]+metric.normals[i,2]*length]
        rng = vtk.vtkMinimalStandardRandomSequence()
        rng.SetSeed(8775070)  # For testing.
        # Compute a basis
        normalizedX = [0 for i in range(3)]
        normalizedY = [0 for i in range(3)]
        normalizedZ = [0 for i in range(3)]
        # The X axis is a vector from start to end
        vtk.vtkMath.Subtract(endPoint, startPoint, normalizedX)
        length = vtk.vtkMath.Norm(normalizedX)
        vtk.vtkMath.Normalize(normalizedX)
        # The Z axis is an arbitrary vector cross X
        arbitrary = [0 for i in range(3)]
        for j in range(0, 3):
            rng.Next()
            arbitrary[j] = rng.GetRangeValue(-10, 10)
        vtk.vtkMath.Cross(normalizedX, arbitrary, normalizedZ)
        vtk.vtkMath.Normalize(normalizedZ)
        # The Y axis is Z cross X
        vtk.vtkMath.Cross(normalizedZ, normalizedX, normalizedY)
        matrix = vtk.vtkMatrix4x4()
        # Create the direction cosine matrix
        matrix.Identity()
        for j in range(0, 3):
            matrix.SetElement(j, 0, normalizedX[j])
            matrix.SetElement(j, 1, normalizedY[j])
            matrix.SetElement(j, 2, normalizedZ[j])
        # Apply the transforms
        transform = vtk.vtkTransform()
        transform.Translate(startPoint)
        transform.Concatenate(matrix)
        transform.Scale(length, length, length)
        # Transform the polydata
        transformPD = vtk.vtkTransformPolyDataFilter()
        transformPD.SetTransform(transform)
        transformPD.SetInputConnection(normal.GetOutputPort())
        # Create a mapper and actor for the arrow
        normalMapper = vtk.vtkPolyDataMapper()
        normalActor = vtk.vtkActor()
        USER_MATRIX = True
        if USER_MATRIX:
            normalMapper.SetInputConnection(normal.GetOutputPort())
            normalActor.SetUserMatrix(transform.GetMatrix())
        else:
            normalMapper.SetInputConnection(transformPD.GetOutputPort())
        normalActor.SetMapper(normalMapper)
        normalActor.GetProperty().SetColor(255.0 / 255, 0.0 / 255, 0.0 / 255)
        renderer.AddActor(normalActor)
