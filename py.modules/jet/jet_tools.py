# -*- coding: utf-8 -*-
"""
Created on Thu Jun 16 20:50:06 2022

@author: DICH
"""

import numpy as np


import numba  
import numba as nb


def minus_invAz_np(N,D,M,zp,x00,x0d,x0zd_out):
    
    for m in range(M):
        r=x00;        
        zm=zp[m]
        for d in range(D):
            x0zd_out[m,d]=r=zm*r+x0d[d]
        


# -Az^-1
def _minus_invAz(N,D,M,zp,x00,x0d,x0zd_out):
    
    for n in nb.prange(N):
        
        for m in nb.prange(M):
            r=x00[n]
            zm=zp[m]
            for d in range(D):
                x0zd_out[m,d,n]=r=zm*r+x0d[d,n]
    pass



def _minus_invAz_s(N,D,M,zp,x00,x0d,x0zd_out):
    
    for n in range(N):
        
        for m in range(M):
            r=x00[n]
            zm=zp[m]
            for d in range(D):
                x0zd_out[m,d,n]=r=zm*r+x0d[d,n]
    pass


@numba.jit(nopython=True,nogil=True,cache=True)
def op_set(d,n,x,s):
    x[d,n]=-s;
    

@numba.jit(nopython=True,nogil=True,cache=True)
def op_sub(d,n,x,s):
    x[d,n]-=s;

@numba.jit(nopython=True,nogil=True,cache=True)
def op_negsub(d,n,x,s):
    x[d,n]=-x[d,n]-s;


prange=range
#@numba.jit(nopython=True,parallel=True,nogil=True,cache=True)
def _assembly_jetz_op(op,N,D,M,res,ZD,xz,xxzD,xtD):
    
    for n in prange(N):
        for d in range(D):
            s=0j;            
            for m in range(M):
                s+=res[m]*(ZD[d,m]*xz[m][n]-xxzD[m][d][n])
            op(d+1,n,xtD,s)
            #xtD[d-1,n]=-s
            
        s=0j        
        for m in range(M):
                s+=res[m]*xz[m][n]
        op(0,n,xtD,s)
        #xtD[0,n]=-s
        
    pass



def assembly_jetz_np(N,D,M,res,ZD,xz,xxzD,xtD):
    #xt=xtD[0];
    
    for d in range(D):
        s=0j;            
        for m in range(M):
            s+=res[m]*(ZD[d,m]*xz[m]-xxzD[m][d])
        #op(d+1,n,xtD,s)
        xtD[d+1]=-s
        
    s=0j        
    for m in range(M):
            s+=res[m]*xz[m]
    #op(0,n,xtD,s)
    xtD[0]=-s
        
    pass



def _assembly_jetz_set(N,D,M,res,ZD,xz,xxzD,xtD):
    #xt=xtD[0];
    for n in nb.prange(N):
        for d in range(0,D):
            s=0j;            
            for m in range(M):
                s+=res[m]*(ZD[d,m]*xz[m][n]-xxzD[m][d][n])
            #op(d+1,n,xtD,s)
            xtD[d+1,n]=-s
            
        s=0j        
        for m in range(M):
                s+=res[m]*xz[m][n]
        #op(0,n,xtD,s)
        xtD[0,n]=-s
        
    pass

def _assembly_jetz_sets(N,D,M,res,ZD,xz,xxzD,xtD):
    #xt=xtD[0];
    for n in range(N):
        for d in range(D):
            s=0j;            
            for m in range(M):
                s+=res[m]*(ZD[d,m]*xz[m][n]-xxzD[m][d][n])
            #op(d+1,n,xtD,s)
            xtD[d+1,n]=-s
            
        s=0j        
        for m in range(M):
                s+=res[m]*xz[m][n]
        #op(0,n,xtD,s)
        xtD[0,n]=-s
        
    pass


prange=nb.prange

minus_invAz=numba.jit(
    parallel=True,
    nopython=True,
    nogil=True,
    cache=True)(_minus_invAz)

minus_invAz_s=numba.jit(
    parallel=True,
    nopython=True,
    nogil=True,
    cache=True)(_minus_invAz_s)

#assembly_jetz_ops=numba.jit(nopython=True,parallel=False,nogil=True,cache=True)(_assembly_jetz_op)

assembly_jetz_op=numba.jit(nopython=True,parallel=True,nogil=True,cache=True)(_assembly_jetz_op)
assembly_jetz_set=numba.jit(nopython=True,parallel=True,nogil=True,cache=True)(_assembly_jetz_set)
assembly_jetz_sets=numba.jit(nopython=True,parallel=False,nogil=True,cache=True)(_assembly_jetz_sets)







if __name__=='__main__': 
    
    def repf(fu,n=1000):
        for k in range(n):
            fu();
            
    @numba.jit(nopython=True,nogil=True,cache=True)    
    def repfn(fu,ls,n=1000):
        for k in range(n):
            fu(*ls);
    
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
    
    
    
    N=16*1024*1024;
    N=128*1024;
    N=2*128*1024;
    #N=8*1024;
    #N=1*1024*1024;
    M=8;
    D=4;
    
    #def _jetx0z(N,D,M,zp,x00,x0d,x0zd_out):
    zp=crand(M)    
    zD=np.array([zp**(k+1) for k in range(D) ])
    xz=crand(M,N)    
    xtD=crand(D+1,N)
    resp=crand(M)    
    x00=crand(N);
    x0d=crand(D,N);
    x0zd_out=0*crand(M,D,N);
    x0zd_out2=0*crand(M,D,N);
    xxzD=x0zd_out
    
    
    tic()
    #minus_invAz(N,D,M,zp,x00,x0d,x0zd_out)
    minus_invAz(2,2,2,zp,x00,x0d,x0zd_out)
    toc('minus_invAz 0')
    tic()
    minus_invAz(N,D,M,zp,x00,x0d,x0zd_out)
    toc('minus_invAz 1')
    
    tic()
    #_assembly_jetz_op(op_set,N,D,M,resp,zD,xz,xxzD,xtD)
    toc()    
        