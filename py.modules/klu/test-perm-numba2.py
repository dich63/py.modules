# -*- coding: utf-8 -*-
"""
Created on Mon May 23 21:35:37 2022

@author: wwww
"""

import numba  
import numba as nb

import numpy as np
import scipy
import scipy.sparse as sp
from utils import *
from jsonrpc.json_io import *

norm=np.linalg.norm


@numba.jit("int32(int32,complex128[:],complex128[:],int32[:])"
           ,nopython=True,
           parallel=True,
           nogil=True,
           fastmath=True           
           )
def permute(n,s,d,ip):
    for i in numba.prange(n):
        d[i]=s[ip[i]];
    return n;


@numba.jit("int32(int32,int32,complex128[:,:],complex128[:,:],int32[:])"
           ,nopython=True,
           parallel=True,
           nogil=True,
           fastmath=True           
           )
def permute2(m,n,s,d,ip):
    for i in numba.prange(n):
        for k in range(m):
            d[k,i]=s[k,ip[i]];
    return n;

@numba.jit("int32(int32,int32,complex128[:,:],complex128[:,:],int32[:,:])"
           ,nopython=True,
           parallel=True,
           nogil=True,
           fastmath=True           
           )
def permute3(m,n,s,d,ip):
    for i in numba.prange(n):
        for k in range(m):
            d[k,i]=s[k,ip[k,i]];
    return n;

@numba.jit("int32(int32,int32,complex128[:,:],complex128[:,:],int32[:,:],float64[:,:])"
           ,nopython=True,
           parallel=True,
           nogil=True,
           fastmath=True           
           )
def permute3_rescale(m,n,s,d,ip,rs):
    for i in numba.prange(n):
        for k in range(m):
            j=ip[k,i]
            d[k,i]=s[k,j]/rs[k,j];
    return n;



#@numba.jit(nopython=True,nogil=True,cache=True)
@numba.jit("int32[:](int32[:])",nopython=True,nogil=True,cache=True)
def iperm(i):
    n=len(i);
    ie=np.arange(n,dtype=np.int32);
    ie[i]=np.arange(n);
    return ie;


@numba.jit("int32(int32,int32[:],int32[:],int32[:],int32[:],int32[:],int32[:],int32[:],int32[:])",          
         nopython=True,nogil=True)
def _csr_permute(n,Ap,Ai,Ax,Cp,Ci,Cx,iP,Q):
    
    nz=0;
    
    for k in range(n):
        Cp[k]=nz;
        j=Q[k];
        for t in range(Ap[j],Ap[j+1]):
            Cx[nz]=Ax[t]
            Ci[nz]=iP[Ai[t]]
            nz+=1
        Cp[n]=nz;
    
    return nz;
    
def csr_permute(A,q,p):
    
    #A=sp.csc_matrix(A)
    n=A.shape[0];
    Ap,Ai,Ax=A.indptr,A.indices,A.data    
    [Cp,Ci,Cx]=[np.empty_like(s) for s in (Ap,Ai,Ax)]
    p=np.array(p,dtype=np.int32)
    q=np.array(q,dtype=np.int32)
    ip=iperm(p)      
    _csr_permute(n,Ap,Ai,Ax,Cp,Ci,Cx,ip,q);
    
    '''
    iq=iperm(q)    
    _csc_permute(np.int32(n),Ap,Ai,Ax,Cp,Ci,Cx,p,iq);    
    '''    
    return A.__class__((Cx,Ci,Cp),shape=A.shape);
    #return sp.csc_matrix((Cx,Ci,Cp),shape=A.shape);
    
    
    
    


#@numba.jit("int32[:](int32[:])",nopython=True,nogil=True,cache=True)
def iperms(i):
    n=len(i);
    ie=np.arange(n);
    ie[i]=np.arange(n);
    return ie;

@numba.jit("int32[:](int32[:])",nopython=True,nogil=True,cache=True)
def jlen(i):    
    return i

m=8
n=1024*128
#n=1024*1024*2

fn=r'O:\__ss\matrix\sfFEM1024k.json'
#fn=r'O:\__ss\matrix\sfFEM128k.json'
#fn=r'O:\__ss\matrix\sfFEM24k.json'
d=decode(fn,1)
C=d.K
C=C.tocsr();
C=C+sp.tril(C)
data=np.arange(1,C.nnz+1)
A=sp.csr_matrix((data,C.indices,C.indptr));



n=A.shape[0];i=np.random.permutation(np.arange(n));k=np.random.permutation(np.arange(n));

tic();c=csr_permute(A,i,k);toc('sp')

tic();t=A[:,k]; cc=t[i,:];toc(':')

ee=c-cc
print('nnz:',ee.nnz)


#raise SystemExit()



i=np.random.permutation(np.arange(n));#ie=np.arange(n);ie[i]=np.arange(n)

iperm(np.arange(3))

'''
s=np.arange(n,dtype=complex)
d=np.arange(n,dtype=complex)
dp=np.arange(n,dtype=complex)
print('start')
tic();d[:]=s[i];toc('perm2')

tic();permute(n,s,dp,i);toc('ppperm2')

'''

s=np.random.randn(m,n)*(1+1j)
d=np.random.randn(m,n)*2j
dp=np.random.randn(m,n)*(1+2j)



print('start')
tic();d[:,:]=s[:,i];t1=toc('perm2')

tic();permute2(m,n,s,dp,i);t2=toc('ppperm2')

print('perf=',t1/t2)



'''
n=5;i=np.random.permutation(np.arange(n));ie=np.arange(n);ie[i]=np.arange(n)
xx=np.eye(n);xp=xx[:,i]
sx=sp.csc_matrix(xx)

px=sx[i,:]
print(px.todense())
'''


