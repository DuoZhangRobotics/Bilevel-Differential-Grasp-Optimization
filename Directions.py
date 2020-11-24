import numpy as np

def cmpToKey(mycmp):
    """Convert a cmp= function into a key= function"""

    class K(object):
        def __init__(self, obj, *args):
            self.obj = obj

        def __lt__(self, other):
            return mycmp(self.obj, other.obj) < 0

        def __gt__(self, other):
            return mycmp(self.obj, other.obj) > 0

        def __eq__(self, other):
            return mycmp(self.obj, other.obj) == 0

        def __le__(self, other):
            return mycmp(self.obj, other.obj) <= 0

        def __ge__(self, other):
            return mycmp(self.obj, other.obj) >= 0

        def __ne__(self, other):
            return mycmp(self.obj, other.obj) != 0

    return K

class Directions:
    def __init__(self, res=4, dim=6):
        self.dirs = []
        self.res = res
        self.dim = dim
        n = np.array([0.0] * dim, dtype=np.float64)
        for d in range(dim):
            n[d] = 1
            self.addDirs(0, d, n)
            n[d] = -1
            self.addDirs(0, d, n)

        # sort dir
        def cmp(A, B):
            for d in range(dim):
                if A[d] < B[d]:
                    return -1
                elif A[d] > B[d]:
                    return 1
            return 0

        self.dirs = sorted(self.dirs, key=cmpToKey(cmp))
        # make compact
        j = 0
        for i in range(len(self.dirs)):
            if i > 0 and (self.dirs[i] == self.dirs[i - 1]).all():
                continue
            else:
                self.dirs[j] = self.dirs[i]
                j += 1
        self.dirs = self.dirs[0:j]
        # normalize
        for i in range(len(self.dirs)):
            self.dirs[i] /= np.linalg.norm(self.dirs[i])

    def addDirs(self, d0, d, n):
        if d0 == self.dim:
            self.dirs.append(n.astype(np.float64))
        elif d0 == d:
            self.addDirs(d0 + 1, d, n)
        else:
            for i in range(self.res):
                n[d0] = -1 + 2 * i / float(self.res - 1)
                self.addDirs(d0 + 1, d, n)

    def printDirs(self):
        print('res=%d dim=%d #dir=%d' % (self.res, self.dim, len(self.dirs)))
        for d in self.dirs:
            print(d)

if __name__ == '__main__':
    Directions(res=3, dim=3).printDirs()
