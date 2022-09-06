# -*- coding: utf-8 -*-
"""
Created on Sun Jul  3 18:37:01 2022

@author: DICH
"""
import numba 
#numba.config.THREADING_LAYER='threadsafe'

import numpy as np
from lipa.pade_exp_poles_res import get_poles_res_array
from jet.symb_pade_exp import exp_poles_res,exp_poles_res_dt_ZD,pade_exp_poles_res
#from klu.iklu import *
#import klu.iklu  as iklu
from klu.klu_numba import *
from jet.nb_tools import * 
from FEM.FEM2sparse import *


def set_zeros_MN(y):
    M,N=y.shape
    set_zeros2(M,N,y);

class jet_csr_klu_t(object):
    
     def __init__(this,spmb,maskDiff=None,Dmax=None,freal=False,dtype=complex,**opts):
         
         def _find1i(im):
             try:
                 return tuple(im).index(1);
             except:
                 return -1;
                 
         this.freal=freal;
         this.spmb=spmb;
         imask=np.arange(len(spmb.data))  if maskDiff is None else np.array(maskDiff)
         
         this.ind1=_find1i(imask);
         
         this.Dm=Dm=max(imask)-1
         D =  Dm if not Dmax else max(Dm,Dmax);
         
             
         this.D,this.imask=D,imask;
         
         #this.common=common=common_update(**opts);
         this.common=common=common_update(ordering=0,scale=0,btf=0,xtol=0.1)
         
         tclass=np.int64(0x0109)
         N=spmb.shape[0]
         
         this.hsymbolic=hsymbolic=iklu_symbolic_batch(N,spmb.indptr,spmb.indices,tclass,common)
         this.state=hsymbolic.state
         
         
         this._xx=np.zeros((D+1,N),dtype=dtype);
         
         if freal:
             this.assembly_jetz_op= assembly_jetz_op2_D_real if Dmax else assembly_jetz_op2_0_real
         else:
             this.assembly_jetz_op= assembly_jetz_op2_D if Dmax else assembly_jetz_op2_0
             
         this.dt=None;
         
         this.singular_source,this.source=None,None
         
         
         
     def preset(this,dt=1,LM=(4,4)):    
         L,M=this.LM=LM;
         D,imask=this.D,this.imask;
         
         
         this.dt=dt;                
         
         
         
         
         poles,res=pade_exp_poles_res(LM,t=dt,fhalf=this.freal)
         
         this.Mp=Mp=len(poles)
         #poles,res=get_poles_res_array(LM); poles,res=poles/dt,res/dt
         
         zD=np.array([poles**deg for deg in imask]).T
         
         this.ZD=ZD=np.array([poles**d for d in range(1,D+1)]) 
         
         #this.ZD1=np.array([poles1**d for d in range(1,D+1)])
         this.scaleD=np.array([dt**d for d in range(1,D+1)])         
         
         
         
         this.poles_res=poles,res
         nnz=this.spmb.indptr[-1];
         
         Hz=this.Hz=np.empty((Mp,nnz),dtype=np.complex128);         
         
         datas=this.spmb.data;  
         this.Dm=Dm=len(datas)
         
         makeHz(nnz,Dm,Mp,zD,datas,Hz)
         
         this.N=N=this.spmb.shape[0]
         
         this.states=-np.ones((Mp,),dtype=np.int32);
         this.xx0zd=np.empty((Mp,D,N),dtype=np.complex128)
         this.yz=np.zeros((Mp,N),dtype=np.complex128)
         
         
         if L<M:
             this.asm_op=op_set
         elif M&1==0:
             this.asm_op=op_sub
         else:
             this.asm_op=op_negsub
                 
         
         this.reset_singular_source(this.singular_source)
         this.reset_source(this.source);          
         
         
         
         return this;
         
         
         
     def reset(this,dt=1,LM=(4,4)):
         return this.preset(dt=dt,LM=LM).factorize();
         
         
         
     def factorize(this): 
         
         N=this.N
         Hz=this.Hz
         
         hsymbolic=this.hsymbolic
         (d,i,p)=this.spmb.triplet
         
         tclass=np.int64(0x0109)
         
         this.hfactors=iklu_numeric_batch(N,p,i,Hz,tclass,hsymbolic,parallel=True)                       
         
         return this;
     
     def _step(this,parallel=True):
         
         N,D,Mp=this.N,this.D,this.Mp
         
         asm_op=this.asm_op
         
         zp,Bz=this.poles_res;
         
         
         
         factors=this.hfactors;
         
         '''
         if Fz is None:
             Fz=np.zeros((N,M),dtype=complex)
         '''
         
         xx,xx0zd,yz=this._xx,this.xx0zd,this.yz
         x00=xx[0]
         x0D=xx[1:]
         
         imask=this.imask
         Dm,ind1=this.Dm,this.ind1;         
         
         (d,i,p)=this.spmb.triplet         
         
         minus_invAz(N,D,Mp,zp,x00,x0D,xx0zd)
         
         csr_gaxpy_jetz2(N,Dm,Mp,p,i,d,ind1,imask,x00,xx0zd,yz)                  
         
         iklu_solve_ts_batch(factors,yz,parallel=parallel,states=this.states);
         
         this.assembly_jetz_op(asm_op,N,D,Mp,Bz,this.ZD,yz,x00,xx0zd,xx);
         
         
     def step(this,parallel=True):
         #this.yz[:]=0;         
         
         # set breakpoint         
         #import pdb; pdb.set_trace()    
         
         if this.source:
             this.source.fill_fz(this.yz).tic();
         else:
             set_zeros_MN(this.yz)
             
         
         
         if this.singular_source:
             # set breakpoint
             #import pdb; pdb.set_trace()
             this.singular_source.add_fz(this.yz).tic();
               
         
         this._step(parallel=parallel)
         
     
         
     def _reset_source(this,source=None,F=None):
         if source:
             if this.dt:
                 source.set_zt(this.poles_res[0],this.dt);                 
             if F is not None :
                 source.set_F(F);
         return source;
         
     def reset_source(this,source=None,F=None):
         this.source=this._reset_source(source,F)
         return this
         
             
     def reset_singular_source(this,source=None,F=None):
         this.singular_source=this._reset_source(source,F)
         return this
         
     reset_J=reset_source    
     reset_singular_J=reset_singular_source
     
     def __call__(this,rep=1,parallel=True):
         
         for k in range(rep):
             this.step(parallel)
             
         return this._xx;
     
     def dump(this,N,rep=1,parallel=True):
         
         xx=this._xx
         xxn=np.empty((N,)+xx.shape,dtype=xx.dtype);
         
         for n in range(N):
             xxn[n]=this(rep);
         return xxn;
         
        
             
     
     @property    
     def xx(this):
         return this._xx;
    
     @xx.setter
     def xx(this,v):
         this._xx[:]=v;    
         
     @property    
     def x(this):
         return this._xx[0];
    
     @x.setter
     def x(this,v):
         this._xx[1:]=0
         this._xx[0]=v;    
         

import scipy
import scipy.sparse as sp
import LIPA2.qp as qps

class source_csr_t(object):
    
    def __init__(this,source_matrix_t=qps.qp_matrix_t,dtype=np.complex128):                
        this.dtype=dtype        
        this.source_matrix_t=source_matrix_t
        this.set_qp();
        
    def __bool__(this):
        return not not this.qpm
    
    def set_qp(this,qp=None,g=None,t0=None):        
        this.qpm=this.source_matrix_t(qp=qp,g=g,t0=t0,dtype=this.dtype)
        return this
    
    def set_F(this,F,N=None):        
        this.F=sp.csr_matrix(F,dtype=this.dtype);
        s=this.F.shape
        if (N is not None) and (s[1]!=N):
            raise Exception("source.F.shape[1]=%d mismatch with N=%d"%(s[1],N))                           
            
        return this
        
    def set_zt(this,z,dt):
        this.z,this.t=z,dt;               
        this.qpl=this.qpm.reset(dt).laplace(z,dtype=this.dtype)
        return this
        
    def tic(this):    
        this.qpm()
        return this;
    
    def fill_fz(this,fz):
        M,N=fz.shape
        F=this.F;
        cfz=this.qpl.imz;
        #csr_gaxsy_NM(N,M,F.indptr,F.indices,F.data,cfz,fz)
        fz[:]=cfz@F
        return this;
    
    def add_fz(this,fz):
        M,N=fz.shape
        F=this.F;
        cfz=this.qpl.imz;
        #csr_gaxpy_NM(N,M,F.indptr,F.indices,F.data,cfz,fz)
        fz[:]=fz+cfz@F
        return this;



def jet_csr_klu_number(DC,D=None,FF=None,qp=None,g=None,t0=None):
    
    import FEM.FEM2sparse as fem
    import scipy
    import scipy.sparse as sp
    
    def toa(x):
        return sp.coo_matrix(x,dtype=complex);
    
    DC=fem.sparse2coo_matrix_batch(*[toa(m) for m in DC]).tocsr()
    source=source_csr_t().set_qp(qp,g,t0).set_F(FF);    
    jet=jet_csr_klu_t(DC,Dmax=D);
    return jet.reset_J(source)
    

def source_csr(qp=None,g=None,t0=None,dtype=np.complex128):
    return source_csr_t(source_matrix_t=qps.qp_matrix_t,dtype=dtype).set_qp(qp,g,t0)
    
def source_singular_csr(qp=None,g=None,t0=None,dtype=np.complex128):
    return source_csr_t(source_matrix_t=qps.schz_matrix_t,dtype=dtype).set_qp(qp,g,t0)
        
         
if __name__=='__main__':
    
    from utils import *
    
    from jet.tools import *
    import LIPA2.tools as l2ts
    from FEM.FEM2sparse import *
    
    from jsonrpc.json_io import *
    import os
    import scipy
    import scipy.sparse as sp
    
    
    import jet.csr_tools
    
    #jet.csr_tools.csr_matrix =sp.csr_matrix
    
    from jet.csr_tools import *
    from jet.nb_tools import *
    #import jet.nb_tools as nbt
    
    
    
    nan=np.nan
    cnan=nan+1j*nan
    norm=np.linalg.norm
    normm= lambda x:  norm(x.reshape(-1),ord=np.inf)
    crand=lambda *ls : np.random.randn(*ls)+1j*np.random.randn(*ls)
    
    fn=r'O:\__ss\matrix\sfFEM1024k.json'
    fn=r'O:\__ss\matrix\sfFEM128k.json'
    #fn=r'O:\__ss\matrix\KGM.json'
    fn=r'O:\__ss\matrix\sfFEM24k.json'
    #fn=r'O:\__ss\matrix\sfFEM64k.json'
    d=decode(fn,1)
    mK,mG,mM=d.K,d.G,d.M
    mG=0*mG;
    '''
    mK,mG,mM=[m.tocoo() for m in [mK,mG,mM] ]
    datas=[mK.data,mG.data,mM.data]
    
    #coo_matrix_batch_t(FEM_datas,row,col,shape=shape); 
    
    CD=coo_matrix_batch_t(datas,mK.row,mK.col)
    CDcsr=CD.tocsr(dtype=complex)
    N=CD.shape[0]
    '''
    #q=mK.tocsr().copy();q[111,1122]=7;mK=q.tocoo()
    
    tic();CD=sparse2coo_matrix_batch(mK,mG,mM);mtoc('sparse2coo_matrix_batch:')
    tic();CDcsr=CD.tocsr(dtype=complex);mtoc('CD.tocsr(dtype=complex):')
    
    
    
    j=jet_csr_klu_t(CDcsr)
    #%timeit
    j.reset(LM=(3,4))
    j.step(1)

