# -*- coding: utf-8 -*-
"""
Created on Sun Apr 17 04:11:07 2022

@author: wwww
"""
from parallel.sparse import *
from jsonrpc.json_io import *
from utils import *
from LIPA2.qp_solver import *
import numpy as np
from scipy.sparse.linalg import splu
norm=np.linalg.norm

import numba


@numba.jit(nopython=True, parallel=True)
def csrMult_parallel(x,Adata,Aindices,Aindptr,Ashape,xout=None): 

    numRowsA = Ashape[0]    
    
    Ax = np.zeros(numRowsA) if xout is None else xout

    for i in numba.prange(numRowsA):
        Ax_i = 0.0        
        for dataIdx in range(Aindptr[i],Aindptr[i+1]):

            j = Aindices[dataIdx]
            Ax_i += Adata[dataIdx]*x[j]

        Ax[i] = Ax_i            

    return Ax


z=7+17j


d=decode(r'O:\__ss\matrix\FEM128k.json',1)

Az=d.K+z*d.G+z**2*d.M

'''
ac=decode('V:\\work\\MID3D\\AC.json',1)

z=7+17j
#z=z/1000
Az=ac.A+z*ac.C
#Az=Az.tocsc()

#Az=decode(r'O:\__ss\matrix\Az-wave1D-1M.json',1).tocsc()
'''


Az=Az.tocsr()
N=Az.shape[0]

x0=np.random.rand(N);
rAz=Az.real

rp=10
tic()
for k in range(rp):
    y=rAz.dot(x0)
ts=toc('s')



xout=0*np.random.rand(N).real;

yc=  csrMult_parallel(x0,rAz.data,rAz.indices,rAz.indptr,rAz.shape,xout)
yc=  csrMult_parallel(x0,rAz.data,rAz.indices,rAz.indptr,rAz.shape)


tic()
for k in range(rp):
    yc=  csrMult_parallel(x0,rAz.data,rAz.indices,rAz.indptr,rAz.shape,xout)
tp=toc('sc')


print('err_diff=',norm(y-yc)/norm(y))
print('perf=',ts/tp)



'''
print('err_diff=',norm(y-yc)/norm(y))

print('err_c=',norm(x0-Az.dot(yc))/norm(x0))

print('err=',norm(x0-Az.dot(y))/norm(x0))



print(st)
'''