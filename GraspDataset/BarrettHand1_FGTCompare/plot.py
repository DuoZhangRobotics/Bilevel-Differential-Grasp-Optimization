import matplotlib.pyplot as plt
import pickle,glob,scipy
from PIL import Image

plt.rc('lines', linewidth=2)
plt.rc('pdf', fonttype=42)
plt.rc('ps', fonttype=42)
plt.rc('font', size=15)

def readTime(file,iters=100):
    file=open(file,'r')
    prefix="[INFO] OptimizeSQP %d iterations, average time="%iters
    for l in file.readlines():
        if l.startswith(prefix):
            return float(l[len(prefix):])
    assert False

def readData(a=1,b=8):
    x=[i for i in range(a,b+1)]
    FGT=[readTime("FGT%d.dat"%(i*100)) for i in range(a,b+1)]
    noFGT=[readTime("noFGT%d.dat"%(i*100)) for i in range(a,b+1)]
    return x,FGT,noFGT

if __name__=="__main__":
    x,FGT,noFGT=readData()
    fig=plt.figure()
    
    ax=fig.add_subplot(1,1,1)
    ax.plot(x,FGT,'-o',label='with FGT')
    ax.plot(x,noFGT,'-o',label='without FGT')
    ax.set_xlabel('Rel. Point Cloud Density')
    ax.set_ylabel('Avg. SQP Iteration Cost (s)')
    
    plt.legend(loc='best')
    plt.title('Cost - Point Cloud Density')
    plt.savefig('plot.pdf',bbox_inches='tight',pad_inches=0)
    plt.show()