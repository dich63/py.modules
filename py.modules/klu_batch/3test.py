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

import pypardiso

norm=np.linalg.norm


z=7+17j

fn=r'O:\__ss\matrix\sfFEM128k.json'
fn=r'O:\__ss\matrix\sfFEM24k.json'
fn=r'O:\__ss\matrix\sfFEM1024k.json'
#fn=r'O:\__ss\matrix\sfFEM128k.json'
fn=r'O:\__ss\matrix\sfFEM24k.json'
d=decode(fn,1)

Az=d.K+z*d.G+z**2*d.M
K=d.K;
qq=K.tocsr()
K,G,M=[m.tocsr() for m in [d.K,d.G,d.M]]

QQ=K+G;


'''
ac=decode('V:\\work\\MID3D\\AC.json',1)

z=7+17j
#z=z/1000
Az=ac.A+z*ac.C
Az=Az.tocsc()
K,G,M=Az,Az,Az
#Az=decode(r'O:\__ss\matrix\Az-wave1D-1M.json',1).tocsc()
'''


Az=Az.tocsc()
N=Az.shape[0]


def sp_LU_f(A):    
    #lu=sp.linalg.splu(A
    options=dict(Equil=False,PivotGrowth=True
                                     ,PrintStat=False
                                     ,DiagPivotThresh=0.1
                                     ,ColPerm='MMD AT PLUS A'#'MMD_ATA'#
                                     ,ConditionNumber=False,IterRefine='NOREFINE')
    #    options={}
    #options=dict(Equil=False, IterRefine='SINGLE')
    lu=sp.linalg.splu(A,permc_spec='MMD AT PLUS A',options=options)
    
    return lambda x: lu.solve(x) ;

#step=sp_LU_f(Az)
x0=np.random.randn(N)+1j*np.random.randn(N)

tic()
si=sparse_solver_invoker(Az)
st=si.analyze()
toc('analy:')
#tic()
st=si.factorize()
toc('lu ctx')

tic()
for k in range(1):
    si.y=x0; 
    si.solve()
    yc=si.y

toc('ctx tic')



tic()
step=sp_LU_f(Az)
toc('lu ')

tic()
for k in range(1):
    y=step(x0)    
toc('tic')
Az=Az.tocsc()
Az=Az.tocsr()
tic()
ycc = pypardiso.spsolve(Az,x0)
toc('pardiso:')

print('err_diff=',norm(y-ycc)/norm(y))
print('err_diff=',norm(y-yc)/norm(y))

print('err_c=',norm(x0-Az.dot(yc))/norm(x0))

print('err=',norm(x0-Az.dot(y))/norm(x0))
print('err_p=',norm(x0-Az.dot(ycc))/norm(x0))


print(st)

nd=3;
tic()
DC=[ m.tocsr() for m in  [K, G,M]]
toc('tocsr:')

Czxx_out=np.zeros((nd,N),dtype=np.complex128);
xx0=np.zeros((nd,N),dtype=np.complex128);
xxz_out=np.zeros((nd,N),dtype=np.complex128);
buf=np.zeros(N,dtype=np.complex128);
#Az=d.K+z*d.G+z**2*d.M
#DC[-1]=None
tic()
ts.AzC0(xx0,DC,Czxx_out);
toc('jet0:');ts.AzC1(xx0,DC,z,xxz_out,Czxx_out);toc('jet:')

