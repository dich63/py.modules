# -*- coding: utf-8 -*-
"""
Created on Thu Jun 16 20:50:06 2022

@author: DICH
"""

import numpy as np


import numba  
import numba as nb
import jet.tools  as ts 
import jet.csr_tools  as cts 


cts.prange=ts.prange=nb.prange

op_set=numba.jit(nopython=True,nogil=True,fastmath=False)(ts.op_set)
op_sub=numba.jit(nopython=True,nogil=True,fastmath=False)(ts.op_sub)
op_negsub=numba.jit(nopython=True,nogil=True,fastmath=False)(ts.op_negsub)


minus_invAz=numba.jit(
    parallel=True,
    nopython=True,
    nogil=True,
    fastmath=False,
    cache=True)(ts.minus_invAz)

minus_invAz_s=numba.jit(
    parallel=False,
    nopython=True,
    nogil=True,
    cache=True)(ts.minus_invAz_s)

#assembly_jetz_ops=numba.jit(nopython=True,parallel=False,nogil=True,cache=True)(_assembly_jetz_op)



assembly_jetz_op2_D=numba.jit(nopython=True,
                           parallel=True,
                           fastmath=False,
                           nogil=True)(ts.assembly_jetz_op2_D)

assembly_jetz_op2_0=numba.jit(nopython=True,
                           parallel=True,
                           fastmath=False,
                           nogil=True)(ts.assembly_jetz_op2_0)

assembly_jetz_op2_D_real=numba.jit(nopython=True,
                           parallel=True,
                           fastmath=False,
                           nogil=True)(ts.assembly_jetz_op2_D_real)

assembly_jetz_op2_0_real=numba.jit(nopython=True,
                           parallel=True,
                           fastmath=False,
                           nogil=True)(ts.assembly_jetz_op2_0_real)


assembly_jetz_op=numba.jit(nopython=True,
                           parallel=True,
                           nogil=True,
                           cache=True)(ts.assembly_jetz_op)

assembly_jetz_op2_D_scale=numba.jit(nopython=True,
                           parallel=True,
                           fastmath=False,
                           nogil=True)(ts.assembly_jetz_op2_D_scale)


'''
assembly_jetz_op2=ts.assembly_jetz_op2
'''
#assembly_jetz_set=numba.jit(nopython=True,parallel=True,nogil=True,cache=True)(_assembly_jetz_set)
#assembly_jetz_sets=numba.jit(nopython=True,parallel=False,nogil=True,cache=True)(_assembly_jetz_sets)


makeHz=numba.jit(nopython=True,
                           parallel=True,
                           nogil=True,
                           fastmath=False,
                           cache=True)(cts.makeHz)

csr_gaxpy_jetz=numba.jit(nopython=True,
                           parallel=True,
                           nogil=True,
                           fastmath=False,
                           cache=True)(cts.csr_gaxpy_jetz)

csr_gaxpy_jetz2=numba.jit(nopython=True,
                           parallel=True,
                           nogil=True,
                           fastmath=False,
                           cache=True)(cts.csr_gaxpy_jetz2)


csr_gaxpy_NM=numba.jit(nopython=True,
                           parallel=True,
                           nogil=True,
                           cache=True)(cts.csr_gaxpy_NM)

csr_gaxsy_NM=numba.jit(nopython=True,
                           parallel=True,
                           nogil=True,
                           cache=True)(cts.csr_gaxsy_NM)

set_zeros=numba.jit(nopython=True,
                           parallel=True,
                           nogil=True,
                           cache=True)(ts.set_zeros)

set_zeros2=numba.jit(nopython=True,
                           parallel=True,
                           nogil=True,
                           cache=True)(ts.set_zeros2)
set_zeros2t=numba.jit(nopython=True,
                           parallel=True,
                           nogil=True,
                           cache=True)(ts.set_zeros2t)


if __name__=='__main__': 
    
    def repf(fu,n=1000):
        for k in range(n):
            fu();
            
    @numba.jit(nopython=True,nogil=True,cache=True)    
    def repfn(fu,ls,n=1000):
        for k in range(n):
            fu(*ls);
    #rr=numba.njit()(repfn)
    #from lipa.parallel import LIPA_solver_st,LIPA_solver_pp,lu_solver_factory,solver_test_pp,Tic
    from utils import *
    from utils.c_structs import *
    from LIPA2.qp import *
    #from klu_lipa import *
    from jsonrpc.json_io import *
    import os
    import scipy
    import scipy.sparse as sp
    norm=np.linalg.norm
    normm= lambda x:  norm(x.reshape(-1),ord=np.inf)
    crand=lambda *ls : np.random.randn(*ls)+1j*np.random.randn(*ls)   
    
    
    
    N=16*1024*1024;
    N=128*1024;
    nnz=5*N
    N=8*128*1024;
    #N=8*1024;
    N=4*1024*1024;
    #N=4*128*1024;
    M=8;
    D=2;
    
    #M=4;
    #D=2;
    
    #def _jetx0z(N,D,M,zp,x00,x0d,x0zd_out):
    zp=crand(M)    
    zD=np.array([zp**(k+1) for k in range(D) ])
    xz=crand(M,N)    
    dataD=crand(D,nnz)    
    dataz=crand(M,nnz)    
    xtD=crand(D+1,N)
    resp=crand(M)    
    x00=crand(N);
    x0d=crand(D,N);
    x0zd_out=0*crand(M,D,N);
    x0zd_out2=0*crand(M,D,N);
    xxzD=x0zd_out
    
    A=np.random.rand(8,10,10);
    
    inv=nb.jit()(np.linalg.inv)
    
    @numba.jit#( parallel=True)
    def invB(A,M,N):
        B=np.zeros((M,N,N),dtype=np.float64);
        for m in range(M):
            b=np.linalg.inv(A[m]);
            B[m,:,:]=b;
        return B;
       
    def invBs(A,M,N):
        B=np.zeros((M,N,N),dtype=np.float64);
        for m in range(M):
            b=np.linalg.inv(A[m]);
            B[m,:,:]=b;
        return B;
        
    
    
    
    
    tic()
    #minus_invAz(N,D,M,zp,x00,x0d,x0zd_out)
    minus_invAz(2,2,2,zp,x00,x0d,x0zd_out)
    toc('minus_invAz 0')
    tic()
    minus_invAz(N,D,M,zp,x00,x0d,x0zd_out)
    toc('minus_invAz 1')
    
    tic()
    assembly_jetz_op(op_set,N,D,M,resp,zD,xz,xxzD,xtD)
    toc('assembly_jetz_op')
    
    tic()
    #def makeHz(nnz,D,M,zD,datas,outHz)
    makeHz(nnz,D,M,zD.T,dataD,dataz)     

    toc('makeHz')    
    Ns=7
    F=sp.random(N,Ns,density=1e-4)+1j*sp.random(N,Ns,density=1e-4);
    F=F.tocsr();
    FT=F.T
    cfz=crand(M,Ns)
    
    fz=crand(M,N);
    fzs=crand(M,N);
    
    def smul(F,cfz,fz):
        M,N=fz.shape
        for m in range(M):
            fz[m]=F@cfz[m]
        return fz
    
    def pmul(F,cfz,fz):
        M,N=fz.shape
        csr_gaxsy_NM(N,M,F.indptr,F.indices,F.data,cfz,fz)
        return fz;
    p=pmul(F,cfz,fz)            
    #%timeit p=pmul(F,cfz,fz)
    #%timeit q=smul(F,cfz,fz)
    
    
    import jet.tools  as ts
    '''
    %timeit makeHz(nnz,D,M,zD.T,dataD,dataz)
    %timeit assembly_jetz_op(op_set,N,D,M,resp,zD,xz,xxzD,xtD)
    '''
    
        