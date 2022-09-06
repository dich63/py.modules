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

@numba.jit(nopython=True)
def complex_mul(x,y):
    xr,xi=x.real,x.imag
    yr,yi=y.real,y.imag
    ri=xr*yr-xi*yi
    rr=xr*yi+xi*yr
    return  rr,ri


@numba.jit(nopython=True,
           parallel=True,
           nogil=True,
           fastmath=True,
           cache=True)
def csrMult_parallel(x,Adata,Aindices,Aindptr,Ashape,xout=None): 

    numRowsA = Ashape[0]    
    
    Ax = np.zeros(numRowsA,dtype=np.complex128) if xout is None else xout

    for i in numba.prange(numRowsA):
    
        Ax_i = 0.0
        for dataIdx in range(Aindptr[i],Aindptr[i+1]):

            j = Aindices[dataIdx]
            Ax_i += Adata[dataIdx]*x[j]
            
        Ax[i] = Ax_i            

    return Ax


@numba.jit(nopython=True, parallel=True)
def csrMult_parallel2(x,Adata,Aindices,Aindptr,Ashape,xout=None): 

    numRowsA = Ashape[0]    
    
    
    Ax = np.zeros(2*numRowsA,dtype=np.double) if xout is None else xout

    for i in numba.prange(numRowsA):
    #for i in range(numRowsA):   
        rAx_i = 0.0
        iAx_i = 0.0
        for dataIdx in range(Aindptr[i],Aindptr[i+1]):

            j = Aindices[dataIdx]
            xr=x[j<<1]
            xi=x[(j<<1)+1]
            
            rA=Adata[dataIdx<<1]
            iA=Adata[(dataIdx<<1)+1]
            
            rAx_i += rA*xr-iA*xi
            iAx_i += rA*xi+iA*xr
            
            #Ax_i += Adata[dataIdx]*x[j]
            '''
            rr,ri=complex_mul(Adata[dataIdx],x[j])
            Ax_i.real += rr
            Ax_i.imag += ri
            '''
            

        Ax[i<<1] = rAx_i            
        Ax[(i<<1)+1] = iAx_i            

    return Ax


#def csrMult_parallel_c(x,Adata,Aindices,Aindptr,Ashape,xout=None): 
@numba.jit(forceobj=True)
def csrMult_parallel_c(x,A,xout=None): 
    
    N2=A.shape[0]
    N2=N2<<1
    ND2=A.data.shape[0]*2
    data2=np.frombuffer(A.data.data,count=ND2,dtype=np.double);
    if  xout is None:
        xout=np.zeros(A.shape[0],dtype=np.complex128)    
        
    xout2=np.frombuffer(xout.data,count=N2,dtype=np.double);
    x2=np.frombuffer(x.data,count=N2,dtype=np.double);
    
    csrMult_parallel2(x2,data2,A.indices,A.indptr,A.shape,xout2);
    
    return xout;
    


z=7+17j


d=decode(r'O:\__ss\matrix\FEM1024k.json',1)

Az=d.K+z*d.G+z**2*d.M

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