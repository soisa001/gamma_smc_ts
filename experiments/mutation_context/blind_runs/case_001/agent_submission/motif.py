import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS']:os.environ[k]='4'
import numpy as np,json
from model import load,features,counts
from explore import ROOT
ctx,y,c=load();f,t,rev=features(ctx,y)
eq=((ctx//1024)%4)==((ctx//4)%4)
base=f['cpg'];f['equal2']=base*2+eq
print('equal2 counts',counts(f['equal2'],t,6).tolist())
for j in range(6):
    cc=counts(f['equal2'],t,6)[j];print('group',j,cc.sum(),'probabilities',cc/cc.sum())
for d in [1,2,3]:
    equal=((ctx//(4**(3+d)))%4)==((ctx//(4**(3-d)))%4)
    print('equality',d,'counts',counts(base*2+equal,t,6).tolist())
