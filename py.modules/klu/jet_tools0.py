# -*- coding: utf-8 -*-
"""
Created on Thu Jun 16 20:50:06 2022

@author: DICH
"""

import numpy as np


import numba  
import numba as nb

'''    
except ImportError or ModuleNotFoundError:
    class nb(object):
        pass
    nb=nb();
    nb.prange=range
'''

def jetx0z_np(N,D,M,zp,x00,x0d,x0zd_out):
    
    for m in range(M):
        r=x00;        
        zm=zp[m]
        for d in range(D):
            x0zd_out[m,d]=r=zm*r+x0d[d]
        



# -Az^-1
def _jetx0z(N,D,M,zp,x00,x0d,x0zd_out):
    
    for n in nb.prange(N):
        
        for m in nb.prange(M):
            r=x00[n]
            zm=zp[m]
            for d in range(D):
                x0zd_out[m,d,n]=r=zm*r+x0d[d,n]
    pass

#xsums=numba.jit("(i4,i4,f8[:,:],f8[:])",parallel=False,nopython=True,nogil=True,cache=True)(_xsums)

def _jetx0zs(N,D,M,zp,x00,x0d,x0zd_out):
    
    for n in range(N):
        
        for m in range(M):
            r=x00[n]
            zm=zp[m]
            for d in range(D):
                x0zd_out[m,d,n]=r=zm*r+x0d[d,n]
    pass


@numba.jit('(i4,i4,c16[:,:],c16)',nopython=True,nogil=True,cache=True)
def op_set(d,n,x,s):
    x[d,n]=-s;
    

@numba.jit('(i4,i4,c16[:,:],c16)',nopython=True,nogil=True,cache=True)
def op_sub(d,n,x,s):
    x[d,n]-=s;

@numba.jit('(i4,i4,c16[:,:],c16)',nopython=True,nogil=True,cache=True)
def op_negsub(d,n,x,s):
    x[d,n]=-x[d,n]-s;


def _assembly_jetz_op(op,N,D,M,res,ZD,xz,xxzD,xtD):
    #xt=xtD[0];
    for n in nb.prange(N):
        for d in range(D):
            s=0j;
            zd=ZD[d];
            for m in range(M):
                s+=res[m]*(zd[m]*xz[m][n]-xxzD[m][d])
            op(d+1,n,xtD,s)
            
        s=0j        
        for m in range(M):
                s+=res[m]*xz[m][n]
        op(0,n,xtD,s)
        
    pass


assembly_jetz_op=numba.jit(nopython=True,nogil=True,cache=True)(_assembly_jetz_op)



jetx0z32ss=numba.jit(
    '(i4,i4,i4,c16[:],c16[:],c16[:,:],c16[:,:,:])',
    parallel=False,nopython=True,
    nogil=True,cache=True
    )(_jetx0z)



jetx0z32=numba.jit(
    '(i8,i8,i8,c16[:],c16[:],c16[:,:],c16[:,:,:])',
    parallel=True,nopython=True,
    nogil=True,cache=True
    )(_jetx0z)

jetx0z=numba.jit( 
    [
     '(i4,i4,i4,c16[:],c16[:],c16[:,:],c16[:,:,:])',
     '(i8,i8,i8,c16[:],c16[:],c16[:,:],c16[:,:,:])'
     ],
    parallel=True,nopython=True,
    nogil=True,cache=True
    )(_jetx0z)

jetx0zv=numba.jit(     
    parallel=True,nopython=True,
    nogil=True,cache=True
    )(_jetx0z)


'''        
jetx0z32s=numba.jit(
    '(i4,i4,i4,c16[:],c16[:],c16[:,:],c16[:,:,:])',
    parallel=False,nopython=True,
    nogil=True,cache=True
    )(_jetx0zs)
'''
jetx0z32s=numba.jit(    
    parallel=False,nopython=True,
    nogil=True,cache=True
    )(_jetx0zs)

if __name__=='__main__': 
    
    
    
    #from lipa.parallel import LIPA_solver_st,LIPA_solver_pp,lu_solver_factory,solver_test_pp,Tic
    from utils import *
    from utils.c_structs import *
    #from LIPA2.qp_solver import *
    #from klu_lipa import *
    from jsonrpc.json_io import *
    import os
    import scipy
    import scipy.sparse as sp
    norm=np.linalg.norm
    normm= lambda x:  norm(x.reshape(-1),ord=np.inf)
    crand=lambda *ls : np.random.randn(*ls)+1j*np.random.randn(*ls)   
    
    
    
    N=4*1024*1024;
    N=128*1024;
    M=16;
    D=4;
    
    #def _jetx0z(N,D,M,zp,x00,x0d,x0zd_out):
    print('alloc')    
    zp=crand(M)    
    #zD=[zp**(k+1) for k in range(D) ]
    resp=crand(M)    
    x00=crand(N);
    x0d=crand(D,N);
    x0zd_out=0*crand(M,D,N);
    x0zd_out2=0*crand(M,D,N);
    print('end alloc')
    
    
    tic()
    jetx0z32(N,D,M,zp,x00,x0d,x0zd_out)
    toc('jetx0z32')
    
        