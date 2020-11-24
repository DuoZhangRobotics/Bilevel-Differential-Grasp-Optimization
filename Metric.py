from ConvexHulls import ConvexHull
from Hand import Hand
import scipy, trimesh
import numpy as np

class Metric(object):
    def __init__(self, targets, count=100, friCoef=.7, res=3, M=None):
        self.points = []
        self.normals = []
        self.targets = targets
        self.meshes = [t.mesh() for t in targets]
        samples = [trimesh.sample.sample_surface_even(t, count) for t in self.meshes]
        for tid, sample in enumerate(samples):
            for sid in range(sample[0].shape[0]):
                pt = sample[0][sid,:].tolist()
                fid = sample[1][sid]
                contain = False 
                for tid2, target in enumerate(targets):
                    if tid2!=tid and target.contain(pt):
                        contain=True
                        break
                if not contain:
                    self.points.append(pt)
                    self.normals.append(-self.meshes[tid].face_normals[fid])
        self.mu = friCoef
        self.points=np.array(self.points,dtype=np.float64).T
        self.normals=np.array(self.normals,dtype=np.float64).T
        
        #sample directions
        from Directions import Directions
        self.dirs = Directions(res=res, dim=6).dirs
        self.dirs=np.array(self.dirs,dtype=np.float64).T
        
        #compute gij
        self.gij = np.zeros((self.points.shape[1],self.dirs.shape[1]), dtype=np.float64)
        for i in range(self.gij.shape[0]):
            p=self.points[:,i]
            n=self.normals[:,i]
            for j in range(self.gij.shape[1]):
                d=self.dirs[:,j]
                f2w=np.concatenate((np.eye(3,3,dtype=np.float64),Metric.cross(p)))
                if M is not None:
                    f2w=np.matmul(M,f2w)
                w_perp=(np.mat(d)*f2w*np.mat(n).T)[0,0]
                w_para=np.linalg.norm(np.mat(d)*f2w*(np.identity(3)-np.mat(n).T*np.mat(n)))
                if self.mu*w_perp>w_para:
                    max_sw=(w_perp+w_para**2/w_perp)
                else: max_sw=max(0,w_perp+self.mu*w_para)
                self.gij[i][j]=max_sw
                
    @staticmethod
    def cross(p):
        ret=np.zeros((3,3),dtype=np.float64)
        ret[2,1]= p[0]
        ret[1,2]=-p[0]
        ret[0,2]= p[1]
        ret[2,0]=-p[1]
        ret[1,0]= p[2]
        ret[0,1]=-p[2]
        return ret

    def draw_samples(self, scale, length):
        # show hand
        import vtk
        from Hand import vtk_render,vtk_add_from_hand
        renderer = vtk.vtkRenderer()
        vtk_add_from_hand(self.meshes, renderer, 1.0, use_torch=True)
        vtk_add_metric_samples(renderer, self, scale, length)
        vtk_render(renderer, axes=True)

    def compute_exp_dist(self, link):
        if isinstance(link,Hand):
            return self.compute_exp_dist(link.palm)
        else:
            trans = link.joint_transform_torch.numpy()[0,:,:]
            pss = trans[:3,:3].T @ (self.points.T - trans[:3,3]).T
            ret = [[link.hull.distance_to(pss[:,i]) for i in range(pss.shape[1])]]
            if link.children:
                for c in link.children:
                    ret += self.compute_exp_dist(c)
            return ret

    def compute_metric(self, hand, alpha=.1):
        dists = self.compute_exp_dist(hand)
        dists_value = np.array([[d[0] for d in distsL] for distsL in dists])
        exp_dists_max = np.amax(np.exp(dists_value * -alpha),axis=0)
        Q_value = self.gij * np.expand_dims(exp_dists_max,axis=1)
        Q_value = np.mean(np.amax(Q_value,axis=1),axis=0)
        return Q_value

def vtk_add_metric_samples(renderer, metric, scale, length):
    import vtk
    for i in range(metric.points.shape[1]):
        #point
        sphere = vtk.vtkSphereSource()
        sphere.SetCenter(metric.points[0,i],
                         metric.points[1,i],
                         metric.points[2,i])
        sphere.SetRadius(scale)
        sphere.SetThetaResolution(24)
        sphere.SetPhiResolution(24)
        sphere_mapper = vtk.vtkPolyDataMapper()
        sphere_mapper.SetInputConnection(sphere.GetOutputPort())
        sphere_actor = vtk.vtkActor()
        sphere_actor.SetMapper(sphere_mapper)
        sphere_actor.GetProperty().SetColor(255./255,.0/255,255.0/255)
        renderer.AddActor(sphere_actor)

        #normal
        if length<=0.:
            continue
        normal=vtk.vtkArrowSource()
        normal.SetTipResolution(100)
        normal.SetShaftResolution(100)
        # Generate a random start and end point
        startPoint=[metric.points[0,i],
                    metric.points[1,i],
                    metric.points[2,i]]
        endPoint=[metric.points[0,i]+metric.normals[0,i]*length,
                  metric.points[1,i]+metric.normals[1,i]*length,
                  metric.points[2,i]+metric.normals[2,i]*length]
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
        arbitrary=[0 for i in range(3)]
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
        normalActor.GetProperty().SetColor(255.0/255, 0.0/255, 0.0/255)
        renderer.AddActor(normalActor)