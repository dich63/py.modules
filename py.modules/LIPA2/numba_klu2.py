# -*- coding: utf-8 -*-
"""
Created on Sun Apr 17 04:11:07 2022

@author: wwww
"""
import lipa.pade_exp_poles_res
from parallel.sparse import *
from jsonrpc.json_io import *
from utils import *
from LIPA2.qp_solver import *
import numpy as np
from scipy.sparse.linalg import splu
norm=np.linalg.norm
#from klu.klu_lib import *
from klu.iklu import *
import klu.klu_lib as kl
import numba  
import numba as nb
from numba.experimental import jitclass

spec = [
    ('number', numba.int64),
]

@jitclass(spec)
class some_class:
    def __init__(self, something):
        self.number = something
    def get_num(self):
        return self.number

#@numba.vectorize(["float32(float32)", "float64(float64)"])
@numba.vectorize(["float64(float64,float64)"], target='parallel')
def fsin(r,g):
    return r/(r+g)
    #return np.sin(x);


#@numba.jit(numba.int64(numba.int64[:]),nopython=True,
@numba.jit("int32[:](int32[:],int32,int32)",nopython=True,
           parallel=True,
           nogil=True,
           fastmath=True,
           cache=True
           )
def ppsleep(x,t=100,nn=8):
    for i in numba.prange(nn):
            x[i]=klu_sleep(t)
    return x

@numba.jit("int32[:](int32[:],int32,int32)",nopython=True,
           parallel=False,
           nogil=True,
           fastmath=True
           )
def spsleep(x,t=100,nn=8):
    for i in numba.prange(nn):
            x[i]=klu_sleep(t)
    return x


#@numba.jit(forceobj=True,nogil=True, parallel=True,cache=True)
@numba.jit(nogil=True, parallel=True,cache=True)
def opsleep(x,t=100,nn=8):
    for i in numba.prange(nn):
            x[i]=klu_sleep(t)
    return x


@numba.jit(nogil=True, parallel=True)
def klu_numerics(ispms,ihsymbolic,iphfactor):
    #st=_klu_context_numeric(spm.iptr,0,hsymbolic,None,byref(hfactor))
    nn=ispms.size
    states=np.empty(nn,dtype=np.int32)
    for i in numba.prange(nn):
        kl._klu_context_numeric(ispms,0)
            
    return states


z=7+17j


d=decode(r'O:\__ss\matrix\FEM1024k.json',1)

pr=lipa.pade_exp_poles_res.poles_res(2,2)
ml=10;
#zz=np.array([ml*p[0] for p in pr])
zz=np.array([p[0] for p in pr])
Azz=[ d.K+z*d.G+z**2*d.M for z in zz]
Azz=[ m.tocsc() for m in Azz]
smzz=[ csc_2_sp_t(m) for m in Azz]
      

raise SystemExit();     
'''
ac=decode('V:\\work\\MID3D\\AC.json',1)

z=7+17j
#z=z/1000
Az=ac.A+z*ac.C
#Az=Az.tocsc()

#Az=decode(r'O:\__ss\matrix\Az-wave1D-1M.json',1).tocsc()
'''



Az[0,100]=100;
Az[0,1]=1;
Az[1,0]=-1;
Az=Az.tocsr()



N=Az.shape[0]

x0=np.random.rand(N)+1j*np.random.rand(N);
rAz=Az#.real
y=x0
rep=10
tic()
for k in range(rep):
    y=rAz.dot(x0)
ts=toc('s')



xout=0*x0;
y2=csrMult_parallel_c(x0,rAz,xout);
yc=  csrMult_parallel(x0,rAz.data,rAz.indices,rAz.indptr,rAz.shape,xout)
#yc=  csrMult_parallel(x0,rAz.data,rAz.indices,rAz.indptr,rAz.shape)
#


tic()
for k in range(rep):
    yc=  csrMult_parallel(x0,rAz.data,rAz.indices,rAz.indptr,rAz.shape,xout)
tp=toc('p')
xout=0*x0;

tic()
for k in range(rep):
    ycc=  csrMult_parallel_c(x0,rAz,xout)
tp2=toc('pp')


print('err_diff=',norm(y-yc)/norm(y))
print('err_diff2=',norm(y-ycc)/norm(y))
print('perf=',ts/tp)
print('perf2=',ts/tp2)



'''
print('err_diff=',norm(y-yc)/norm(y))

print('err_c=',norm(x0-Az.dot(yc))/norm(x0))

print('err=',norm(x0-Az.dot(y))/norm(x0))



print(st)
'''