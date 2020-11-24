from ConvexHulls import ConvexHull
import torch, scipy, trimesh
import numpy as np

data_type = torch.double

class Metric(object):
    def __init__(self, targets, count=100, friCoef=.7, res=3, M=None):
        obj = trimesh.boolean.union([t.mesh() for t in targets])
        self.samples, face_indices = trimesh.sample.sample_surface_even(obj, count)
        self.normals = -1 * obj.face_normals[face_indices]
        self.friCoef = friCoef
        
        #sample directions
        from Directions import Directions
        self.dirs = Directions(res=res, dim=6).dirs
        
        #compute gij
        self.gij = np.zeros((len(self.samples),len(self.dirs)),dtype=data_type)
        for i in range(self.gij.shape[0]):
            for j in range(self.gij.shape[1]):
                f2w=np.concatenate((np.eye(3,3,dtype=data_type),cross(self.samples[i])))
                if M is not None:
                    f2w=np.matmul(M,f2w)
                w_perp=(np.mat(self.dirs[j])*f2w*np.mat(self.normals[i]).T)[0,0]
                w_para=np.linalg.norm(np.mat(self.dirs[j])*f2w*(np.identity(3)-np.mat(self.normals[i]).T*np.mat(self.normals[i])))
                if mu*w_perp>w_para:
                    max_sw=(w_perp+w_para**2/w_perp)
                else: max_sw=max(0,w_perp+mu*w_para)
                self.gij[i][j]=max_sw
                