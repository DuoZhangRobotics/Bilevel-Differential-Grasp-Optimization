from sys import *
import re

def read_data(file):
    f=open(file,'r')
    inOurs=False
    inGraspIt=False
    oursAll=[]
    graspItAll=[]
    for line in f.readlines():
        if line.startswith('Ours'):
            inOurs=True
            ours=[]
        elif line.startswith('$Q_1$-\cite{GraspIt}'):
            inGraspIt=True
            graspIt=[]
        else:
            m=re.search(r'[+\-]?(?:0|[1-9]\d*)(?:\.\d*)?(?:[eE][+\-]?\d+)',line)
            if m is not None:
                if inOurs:
                    ours.append(float(m.group(0)))
                if inGraspIt:
                    graspIt.append(float(m.group(0)))
            else:
                if inOurs:
                    inOurs=False
                    oursAll.append(ours)
                if inGraspIt:
                    inGraspIt=False
                    graspItAll.append(graspIt)
    return oursAll,graspItAll
                
if __name__=='__main__':
    oursAll,graspItAll=read_data('resultsTable.tex')
    for ours,graspIt in zip(oursAll,graspItAll):
        improve=[a/b for a,b in zip(ours,graspIt)]
        print("%.2f/%.2f/%.2f"%(min(improve),sum(improve)/len(improve),max(improve)))