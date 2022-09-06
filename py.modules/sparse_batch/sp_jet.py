# -*- coding: utf-8 -*-
"""
Created on Mon Apr 25 03:53:05 2022

@author: wwww
"""
from sparse_batch.sp_types import *
from utils import *



    
_spjet_path=r"O:\__ss2022\sp_jet\x64\Release\sp_jet.dll"    

#_spjet_path=r"O:\__ss2022\sp_jet\x64\debug\sp_jet.dll"    
#_klu_path=r"O:\__ss2022\build\x64\debug\klu_batch.dll" 
#_spjet_path=__file__+r"/../sp_jet.dll"    
_spjet=cdll.LoadLibrary(_spjet_path)




_csr_mulsum=_spjet.csr_mulsum;
_csr_mulsum.rstype=c_int32
_csr_mulsum.argtypes = (p_spm_t,c_void_p,c_void_p,c_long,c_long)


_csr_sort_indices=_spjet.csr_sort_indices_batch;
_csr_sort_indices.rstype=c_int32
_csr_sort_indices.argtypes = (p_spm_t,c_long,c_long)




_csr_mul_zz=_spjet.csr_mul_zz;
_csr_mul_zz.rstype=c_int32
_csr_mul_zz.argtypes = (p_spm_t,c_void_p,p_spm_t,c_long,c_long)



_jetz=_spjet.jetz;
_jetz.rstype=c_int32
_jetz.argtypes = (c_long,c_void_p,c_long,c_long,c_void_p,c_void_p,c_long,c_long,c_void_p)

_csr_mulsum_add_jetz=_spjet.csr_mulsum_add_jetz;
_csr_mulsum_add_jetz.rstype=c_int32
_csr_mulsum_add_jetz.argtypes = (c_long,p_spm_t,c_long,c_void_p,c_void_p,c_long,c_long)

#assembly_jetz(int nz, int D, size_t N, void** ZD, void** xz, void*** xxzD, void** xtD, long nthread, size_t chunk) {

_assembly_jetz=_spjet.assembly_jetz;
_assembly_jetz.rstype=c_int32
_assembly_jetz.argtypes = (c_long,c_long,c_long,c_void_p,c_void_p,c_void_p,c_void_p,c_void_p,c_long,c_long)



_assembly_jetz_op=_spjet.assembly_jetz_op;
_assembly_jetz_op.rstype=c_int32
_assembly_jetz_op.argtypes = (c_long,c_long,c_long,c_long,c_void_p,c_void_p,c_void_p,c_void_p,c_void_p,c_long,c_long)

import numpy as np
import scipy
import scipy.sparse as sp


def make_crs_smzz(spm,z,chunk=11111,nthread=0,smzz_out=None,dtype=np.complex128):
    #def new_spm_like(spm,nbatch=-1,datas=None,dtype=np.complex128)
    #tt=tic('smzz',1)
    sm=spm.ptr.contents;
    nZ=len(z);
    nD=sm.nbatch;
    zz=[ np.array(z,dtype=dtype)**k for k in range(nD) ]
    
    ppzz=ptr_ptr_array(zz);
    
    smzz= new_spm_like(spm,nZ,dtype=dtype) if smzz_out is None else smzz_out
    #toc('smzz',t=tt)
    _csr_mul_zz(spm.ptr,ppzz.ptr,smzz.ptr,nthread,chunk)
    #toc('_csr_mul_zz',t=tt)
    
    
    return smzz;
        
    
         
class xxz_buffer_t(object):
    
     def reset(self,dt=1):
         
         res,zz,zzD,polusres=self.res,self.zz,self.zzD,self.polusres 
         zt,rt=polusres.T
         zz[:]=zt/dt;
         res[:]=rt/dt;
         zzD[:]=np.array([zz**k for k in range(0,self.D)])
         
         return self;
     

         
     #def __init__(self,zz,res,xx,dt=1.0,dtype=np.complex128):
     def __init__(self,polusres,xx,dt=1.0,dtype=np.complex128):
         
         self.polusres=polusres=np.array(polusres,dtype=dtype);
         
         self.D=D=xx.shape[0];   
         self.mp=mp=len(polusres)
         
         self.zz=zz=np.empty(mp,dtype=dtype);
         self.res=res=np.empty(mp,dtype=dtype);
         self.zzD=zzD=np.empty([D,mp],dtype=dtype);
         
               
         self.reset(dt);
         
         
         x0=xx[0];
         self.N=N=len(x0);
         
         
         
         
         self.xx=xx;
         self._ppxx=ppxx=np.empty((mp,D),dtype=object);
         
         for m in range(mp):
             ppxx[m][0]=x0;
             for d in range(1,D):
                 ppxx[m][d]=np.empty(N,dtype=dtype);
                 ppxx[m][d][:]=np.nan
                 
         
         self._lpres=lpres=link_ptr_t(res.ctypes.data,res)
         self._lpz=lpz=link_ptr_t(zz.ctypes.data,zz)
         self._lpxx=lpxx=ptr_ptr_array(xx)
         self._lpxz=lpxz=ptr_ptr_ptr_array(ppxx)
         
         self._lpzzD=lpzzD=ptr_ptr_array(zzD)
         
         self.Bzptr=lpres.ptr;         
         self.ZDptr=lpzzD.ptr;
         self.zptr=lpz.ptr;
         self.xxptr=lpxx.ptr;
         self.xzptr=lpxz.ptr;
         self._nt=c_long(-1)
     
     @property    
     def nt(self):
         return self._nt.value
     
     def __call__(self,chunk=290,nthread=0):
         
         return _jetz(self.mp,self.zptr,self.D,self.N,
                      self.xxptr,self.xzptr,nthread,chunk,byref(self._nt));
     def gaxpyjet(self,sm,ppy,chunk=1,nthread=0,offset=1):
         #HRESULT csr_mulsum_add_jetz(spm_t* spm, int mz, void*** ppxxD, void** ppy, long nthread, size_t chunk)
         mz=len(self.zz)
         return _csr_mulsum_add_jetz(offset,sm.ptr,
                                    mz,self.xzptr,ppy.ptr,
                                    nthread,chunk);
     @property                 
     def xxz(self):
         return np.array([ list(p)  for p in  self._ppxx]);
     
     def assembly(self,pxz,chunk=1,nthread=0):
         D,nz=self.zzD.shape;
         D=D-1;
         N=self.N;         
         st=_assembly_jetz(nz,D,N,self.Bzptr,self.ZDptr,pxz,self.xzptr,self.xxptr, nthread,  chunk);
         return st; 
     
     def assembly_op(self,op,pxz,chunk=1,nthread=0):
         D,nz=self.zzD.shape;
         D=D-1;
         N=self.N;         
         st=_assembly_jetz_op(op,nz,D,N,self.Bzptr,self.ZDptr,pxz,self.xzptr,self.xxptr, nthread,  chunk);
         return st; 
         #_assembly_jetz
         

class sp_view_t(object):
    def __init__(self,spm,*lp):
        self.ptr=spm.ptr;
        self.sm=spm.ptr.contents
        self._lp=lp
        
    @property
    def nbatch(self):
        return self.sm.nbatch;
    def __getitem__(self,n):   
        fmt=self.sm.fmt
        p,i,d=spm_triplet(self.sm,n);
        
        if fmt==b'coo':
            s=sp.coo_matrix((d,(p,i)));
        elif fmt==b'csr':
            s=sp.csr_matrix((d,i,p));
        else:
            s=sp.csc_matrix((d,i,p));    
            
        return s;
         
    
class sp_qp_t(object):
    pass
    

if __name__=='__main__':
    
    
    def sp_mulsum(CC,xx,yout):
        yout[:]=0;
        for k in range(len(CC)):
            yout+=CC[k].dot(xx[k])
        return yout;
    
    
    from lipa.parallel import LIPA_solver_st,LIPA_solver_pp,lu_solver_factory,solver_test_pp,Tic
    import lipa.pade_exp_poles_res
    from utils import *
    from utils.c_structs import *
    #from LIPA2.qp_solver import *
    from LIPA2.tools import *
    #from klu_lipa import *
    from jsonrpc.json_io import *
    import os
    import scipy
    import scipy.sparse as sp
    
    LIPA_solver=LIPA_solver_st;
    
    norm=np.linalg.norm
    normm= lambda x:  norm(x,ord=np.inf)
    
    randnc= lambda *lp:  np.random.randn(*lp)+1j*np.random.randn(*lp)
    
    pr=lipa.pade_exp_poles_res.poles_res(8,8)
    pr=lipa.pade_exp_poles_res.poles_res(4,4)
    zz=np.array([p[0] for p in pr])
    
    A=[[1,2],[3,4]]
    A=[[1,1],[1,1]]
    K,G,M=[(0+1)*sp.csr_matrix(A,dtype=complex) for k in range(3)]
    
    z=np.array([1,1,1,1,1,1,1,1],dtype=complex);
    z=np.array(range(1,7+1),dtype=complex);
    nD=5
    zz0=np.array([z**k  for k in range(nD) ])
    Jet=[(0+1)*sp.csr_matrix(A,dtype=complex) for k in range(nD)]
    
    #zz=list(zz);
    '''
    zz=[1,1]
    zz[0]=np.array([z[0]**k  for k in range(3) ])
    zz[1]=np.array([z[1]**k  for k in range(3) ])
    '''
    zz=[1,1]
    zz[0]=np.array(zz0[0]);
    zz[1]=np.array(zz0[1]);
    
    zz=np.array(zz)
    
    zz=np.ascontiguousarray(zz0.T)
    zz=np.ascontiguousarray(zz0)
    #zz=np.transpose(zz0)
    
    sm=make_spm(*Jet);   
    smzz=new_spm_like(sm,nbatch=len(z))
    ppzz=ptr_ptr_array(zz);
    #_csr_mul_zz.argtypes = (p_spm_t,c_void_p,p_spm_t,c_long,c_long)
    print('pid=',os.getpid())
    _csr_mul_zz(sm.ptr,ppzz.ptr,smzz.ptr,0,11)
    sv=sp_view_t(smzz)
    print(sv[0].todense())
    print(sv[-1].todense())
    #raise SystemExit(0);
    fi=__file__
    fn=r'O:\__ss\matrix\sfFEM1024k.json'
    #fn=r'O:\__ss\matrix\sfFEM128k.json'
    #fn=r'O:\__ss\matrix\sfFEM24k.json'
    d=decode(fn,1)
    
    
    print(fn)
    k=d.K
    k=k.tocsr();
    
    
    
    K,G,M=[sp.csr_matrix(m,dtype=complex) for m in [d.K,d.G,d.M]]
    
    sm=make_spm(K,G,M);   
    
    
    
    
    
    
    '''
    for m in [K,G,M]:
        m.sum_duplicates()
    '''
    N=K.shape[0];
    xx=randnc(3,N);
    y=randnc(N);
    y2=randnc(N);
    
    
    
    print('pid=',os.getpid())
    
    tic(); sp_mulsum((K,G,M),xx,y);t1=toc('sp_mulsum:')
    
    
    
    
    '''
    K.data[:]=1;
    G.data[:]=2;
    M.data[:]=3;
    '''
    sm=make_spm(K,G,M);   
    print('pid=',os.getpid())
    
    tic();_csr_sort_indices(sm.ptr,1,111);toc('_csr_sort_indice:')
    
    raise SystemExit(0);
    
    ppx=ptr_ptr_array(xx);
    py=y2.ctypes.data
    
    nt,chunk=0,290
    tic(); st=_csr_mulsum(sm.ptr,ppx.ptr,py,nt,chunk);t2=toc('csr_mulsum:')
    tic(); st=_csr_mulsum(sm.ptr,ppx.ptr,py,nt,chunk);t2=toc('csr_mulsum:')
    print("nt,chunk=",nt,chunk)
    err=2*norm(y-y2)/(norm(y)+norm(y2))
    print('err',err)
    print('perf',t1/t2)
    
    
    
    
    
    
    '''
    xz=xxz_buffer_t(zz,xx)
    #tic();xz(nthread=8,chunk=5);toc('xz')
    tic();xz(nthread=0,chunk=2000);toc('xz')
    tic();xz(nthread=0,chunk=2000);toc('xz')
    
    zz=0*zz
    
    
    sh=(len(zz),)+xx.shape
    xxz_out=np.empty(sh,dtype=complex);
    yyz_out=np.zeros((len(zz),xx.shape[1]),dtype=complex);
    yyz_outp=np.zeros((len(zz),xx.shape[1]),dtype=complex);
    ppyz=ptr_ptr_array(yyz_out);
    
    tic();st=xz.gaxpyjet(sm,ppyz);toc('gaxpyjet')
    DC=[K,G,M]
    tic()
    AzC0(xx,DC,yyz_outp);
    for k in range(len(zz)):
        xxz_out[k][0]=xx[0]
        AzC1(xx,DC,zz[k],xxz_out[k],yyz_outp);
        
    toc('jet:')
    #errz=2*norm(xz.xxz(0)-xxz_out)/(norm(xxz_out)+norm(xz.xxz(0)))
    #print('errz',errz)


    
    tic(); st=_csr_mulsum(sm.ptr,ppx.ptr,py,nt,chunk);t2=toc('csr_mulsum:')
    tic(); st=_csr_mulsum(sm.ptr,ppx.ptr,py,nt,chunk);t2=toc('csr_mulsum:')
    tic(); st=_csr_mulsum(sm.ptr,ppx.ptr,py,nt,chunk);t2=toc('csr_mulsum:')
    tic(); st=_csr_mulsum(sm.ptr,ppx.ptr,py,nt,chunk);t2=toc('csr_mulsum:')
    tic(); st=_csr_mulsum(sm.ptr,ppx.ptr,py,nt,chunk);t2=toc('csr_mulsum:')
    '''
    
#    raise SystemExit(0);
    