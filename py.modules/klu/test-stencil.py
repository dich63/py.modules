# -*- coding: utf-8 -*-
"""
Created on Mon Jun 13 20:39:50 2022

@author: DICH
"""

import numba
import numba as nb

@numba.stencil
def _smooth(x):
    return (x[-1, -1] + x[-1, 0] + x[-1, 1] +
            x[ 0, -1] + x[ 0, 0] + x[ 0, 1] +
            x[ 1, -1] + x[ 1, 0] + x[ 1, 1]) / 2.0

import numpy as np
x = np.ones((100, 100),dtype=np.int32)
y = np.ones((100, 100),dtype=np.float64)

#%timeit _smooth(x)

@numba.jit("(int32[:,:],f8[:,:])",nopython=True,nogil=True,cache=True)
def test(x,y):
    y[:,:]= _smooth(x)
    return 11;

def ttt(x,y):
    y[:,:]= _smooth(x)
    return 11;


def _xcopy(N,x,y):
    for i in numba.prange(N):
        y[i]=x[i]
    
xcopy=numba.jit("(i4,f8[:],f8[:])",nopython=True,nogil=True,cache=True)(_xcopy)


def _xcopys(N,x,y):
    for i in range(N):
        y[i]=x[i]
    
xcopys=numba.jit("(i4,f8[:],f8[:])",nopython=True,nogil=True,cache=True)(_xcopys)


jit_proto={
    'signature_or_function':'(i4,f8[:],f8)',
    'nopython':True,
    'nogil':True,
    'cache':True
    }

@numba.jit(**jit_proto)
def op_set(i,x,s):
    x[i]=s;
    
#@numba.jit('(i4,f8[:],f8)',nopython=True,nogil=True,cache=True)
@numba.jit(**jit_proto)
def op_plus(i,x,s):
    x[i]+=s;

#@numba.jit((numba.typeof(op_set),nb.i4,nb.i4,nb.f8[:,:],nb.f8[:]),parallel=True,nopython=True,nogil=True,cache=True)

#

prange=range

def _xsum_op(op,N,M,x,y):
    
    for i in prange(N):
        s=0.0        
        for m in range(M):
            s+=x[m,i]
            op(i,y,s);

#xsum_op=numba.jit(parallel=True,nopython=True,nogil=True,cache=True)(_xsum_op)
#xsum_op=numba.jit((nb.typeof(op_set),nb.i4,nb.i4,nb.f8[:,:],nb.f8[:]),parallel=True,nopython=True,nogil=True,cache=True)(_xsum_op)
#xsum_op=numba.jit("(numba.typeof(x),i4,i4,f8[:,:],f8[:])",parallel=True,nopython=True,nogil=True,cache=True)(_xsum_op)

prange=numba.prange

xsum_op=numba.jit([(nb.typeof(op_set),nb.i4,nb.i4,nb.f8[:,:],nb.f8[:]),
                   (nb.typeof(op_plus),nb.i4,nb.i4,nb.f8[:,:],nb.f8[:])
                   ],parallel=True,nopython=True,nogil=True,cache=True)(_xsum_op)

def _xsum(N,M,x,y):
    
    for i in numba.prange(N):
        s=0.0        
        for m in range(M):
            s+=x[m,i]
            y[i]=s

xsum=numba.jit("(i4,i4,f8[:,:],f8[:])",parallel=True,nopython=True,nogil=True,cache=True)(_xsum)
xsumv=numba.jit(parallel=True,nopython=True,nogil=True,cache=True)(_xsum)


def _xsums(N,M,x,y):
    
    for i in range(N):
        s=0.0        
        for m in range(M):
            s+=x[m,i]
            y[i]=s

xsums=numba.jit("(i4,i4,f8[:,:],f8[:])",parallel=False,nopython=True,nogil=True,cache=True)(_xsums)



#@numba.njit
#@numba.jit("int32[:](int32[:])",nopython=True,nogil=True,cache=True)
#@numba.jit(nopython=True,nogil=True,cache=True)
#@numba.jit("(int32[:])",nopython=True,nogil=True,cache=True)
@numba.jit(nopython=True,nogil=True,cache=True)
def smooth(x):
    y= _smooth(x)
    return y

%timeit smooth(x)

s="(int32[:,:],f8[:,:])"
s="(i4[:,:],f8[:,:])"
qq=numba.jit(s,nopython=True,nogil=True,cache=True)(ttt)

N=4000000;M=30;x=np.random.rand(M,N);y=-1*x;

%timeit xsum(N,M,x,y[0])