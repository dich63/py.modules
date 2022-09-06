# -*- coding: utf-8 -*-
"""
Created on Mon Apr 25 03:53:05 2022

@author: wwww
"""

from sparse_batch.sp_types import *

class klu_common_t(Structure):
    _pack_=8
    _fields_ = [
    ("tol",c_double),
    ("memgrow",c_double),
    ("initmem_amd",c_double),
    ("initmem",c_double),
    ("maxwork",c_double),

    ("btf",c_long),
    ("ordering",c_long),  
    ("scale",c_long),
    #("nz",c_longlong),

    ("user_order",c_void_p),  
    ("user_data",c_void_p),  

    ("halt_if_singular",c_long),

    ("status",c_long),  
    ("nrealloc",c_long),

    ("structural_rank",c_long),  

    ("numerical_rank",c_long),

    ("singular_col",c_long),

    ("noffdiag",c_long),

    ("flops",c_double),
    ("rcond",c_double),
    ("condest",c_double),
    ("rgrowth",c_double),
    ("work",c_double),

    ("memusage",c_longlong),
    ("mempeak",c_longlong)

    ];
    
class klu_handle_t(Structure):
      _fields_ = []
  

class klu_factor_t(Structure):
    _pack_=8
    _fields_ = [
    ("Symbolic",POINTER(klu_handle_t)),
    ("Numeric",POINTER(klu_handle_t)),    
    ("Common",klu_common_t)
    ];
    
class klu_batch_t(Structure):
    _pack_=8
    _fields_ = [
    ("Symbolic",POINTER(klu_handle_t)),
    ("pfactors",POINTER(klu_factor_t)),  
    ("nbatch",c_longlong),
    ("Common",klu_common_t)
    ];
    
p_klu_batch_t  = POINTER(klu_batch_t)
pp_klu_batch_t  = POINTER(p_klu_batch_t)


    
_klu_path=r"O:\__ss2022\build\x64\Release\klu_batch.dll"    
_klu_path=r"O:\__ss2022\build\x64\debug\klu_batch.dll"    
_klu=cdll.LoadLibrary(_klu_path)

_klub_create=_klu.klub_create;
_klub_create.rstype=c_int32
_klub_create.argtypes = (pp_klu_batch_t,)


_klub_release=_klu.klub_release;
_klub_release.rstype=c_int32
_klub_release.argtypes = (p_klu_batch_t,)

_klub_addref=_klu.klub_addref;
_klub_addref.rstype=c_int32
_klub_addref.argtypes = (p_klu_batch_t,)

_klub_symbolic=_klu.klub_symbolic;
_klub_symbolic.rstype=c_int32
_klub_symbolic.argtypes = (p_klu_batch_t,p_spm_t)

_klub_numeric=_klu.klub_numeric;
_klub_numeric.rstype=c_int32
_klub_numeric.argtypes = (p_klu_batch_t,p_spm_t,c_int32)


_klub_solve=_klu.klub_solve;
_klub_solve.rstype=c_int32
_klub_solve.argtypes = (p_klu_batch_t,p_rhs_t,c_int32)

_csr2csc_batch=_klu.csr2csc_batch;
_csr2csc_batch.rstype=c_int32
_csr2csc_batch.argtypes = (p_spm_t,p_spm_t,c_int32)



class klu_batch(object):
    
    def __init__(self,**opts):
        self._st=-1;
        self._phase=-1;
        self._hkb=hkb=p_klu_batch_t();
        _klub_create(byref(hkb))       
        
        o=dict(ordering=0,scale=-1,btf=0)
        o.update(opts)
        c_update(hkb.contents.Common,**o)       
        
        
    def __del__(self):
        _klub_release(self._hkb);
        
    @property
    def phase(self):
        return self._phase;
    @property
    def opts(self):
        return c_2_dict(self.hkb.contents.Common);
    
    def symbolic(self,smp):
        self._smp=smp;
        self._st=st=_klub_symbolic(self._hkb,smp.ptr)
        if st==0:
            self._phase=1;
        return self;
    
    def factorize(self,smp,nthreads=0):
        if self._phase<=1:
            self.symbolic(smp);
        if self._st==0:
            self._st=_klub_numeric(self._hkb,smp.ptr,nthreads);
            if self._st==0:
                self._phase==2               
            
        return self;   
    
    def __call__(self,bb,nthreads=0):
        st=self._st
        if st==0:
            st=_klub_solve(self._hkb,bb.ptr,nthreads); 
        return st;
    
            
        
        
        
            
        
    



    

if __name__=='__main__':
    
    
    
    
    
    
    from lipa.parallel import LIPA_solver_st,LIPA_solver_pp,lu_solver_factory,solver_test_pp,Tic
    from utils import *
    from utils.c_structs import *
    #from LIPA2.qp_solver import *
    #from klu_lipa import *
    from jsonrpc.json_io import *
    import os
    import scipy
    import scipy.sparse as sp
    
    LIPA_solver=LIPA_solver_st;
    
    norm=np.linalg.norm
    normm= lambda x:  norm(x,ord=np.inf)
    
    
    
    
    
    fn=r'O:\__ss\matrix\sfFEM1024k.json'
    #
    fn=r'O:\__ss\matrix\sfFEM128k.json'
    #fn=r'O:\__ss\matrix\sfFEM24k.json'
    d=decode(fn,1)
    k=d.K
    k=k.tocsr();
    k=k.tocsr();
    q=k.has_sorted_indices
    
    #K,G,M=[m.tocsr() for m in [d.K,d.G,d.M]]
    
    zz=[7+17j,2+1j,2+4j,2-4j,1+1j,2+1j,2+4j,2-4j]
    zz=[7+17j,2+1j,2+4j,2-4j]
    zz=[7+17j for k in range(8)]
    nt=1
    #zz=[1+1j,2+1j]
    #    zz=[7+17j,]
    
    N=d.M.shape[0]
    di=np.ones(N-1)
    
    #G=sp.tril(G,format='csc')
    #ab=[K+z*G+z**2*M for z in zz ];
    K,G,M=[m for m in [d.K,d.G,d.M]]
    daz=np.empty_like(K.data,dtype=np.complex128)
    ab=[sp.coo_matrix((copy.copy(daz),(K.row,K.col)),shape=K.shape)  for z in zz ]
    
    ds=[K.data+z*G.data+z**2*M.data for z in zz]
    
    for k in range(len(ds)):
        ab[k].data[:]=ds[k];
    
    #ab=[K+z*G+z**2*M for z in zz ];
    
    ab=[m.tocsr() for m in ab ];
    ab=[m.tocsc() for m in ab ];
    
    tic()
    for m in ab:
        m.sum_duplicates()
    toc('sum')
    
    #ab=[m.T+m for m in ab ];
    
    
    #ab=[m.tocsr() for m in ab ];
    Az=ab[0].tocsr()
    
    #    ab=[m.tocoo() for m in ab ];
    #ab=[m.tocsc() for m in ab ];
    
    tic();sm=make_spm(*ab);toc('make_spm')
    '''
    abo=[m.tocoo() for m in ab ];
    
    ab[0].data[:]=np.nan;
    
    smo=make_spm(*abo)
    st=-111
    st=_csr2csc_batch(smo.ptr,sm.ptr)
    
    print('_csr2csc_batch st=',st)
    '''
    
    #ab=[m.tocsc() for m in ab ];
    #$ab=[m.tocoo() for m in ab ];
    #ab[0].data[:]=np.nan;
    
    
    
    
    
    
    
    
    
    print('pid=',os.getpid())
    
    h=p_klu_batch_t()
    
    _klub_create(byref(h))
    
    #c_update(h.contents.Common,ordering=0,scale=-1,btf=0)
    c_update(h.contents.Common,ordering=0,scale=1,btf=1)
    #print(c_2_dict(h.contents.Common))
    
    print(c_2_dict(sm.ptr.contents))
    tic();sts=_klub_symbolic(h,sm.ptr);toc('_klub_symbolic:')
    
    psm=byref(sm.ptr)
    
    #tic();st=_klub_numeric(h,byref(sm.ptr),1);toc('_klub_numeric:')
    #tic()
    
    
    st=_klub_numeric(h,sm.ptr,nt);toc('_klub_numeric:')
    
    print('sts=',sts,'st=',st)
    
    bb=xx_buffer_t(Az.shape[0],len(ab));
    z0=88+77j
    N=Az.shape[0]
    z0=np.random.randn(N)+1j*np.random.randn(N)
    
    bb[0][0,:]=z0;
    
        
    '''
    print('bb._rhs.nrhs=',bb._rhs.nrhs)
    st=_klub_solve(h,byref(bb._rhs),1)
    print('bb._rhs.nrhs=',bb._rhs.nrhs)
    '''
    
    print('pid=',os.getpid())
    tic();st=_klub_solve(h,bb.ptr,nt);toc('_klub_solve:')
    print('st=',st)
    tic()
    err=norm(Az*bb[0][0,:]-z0)/norm(Az*bb[0][0,:])
    toc('err')
    print('fmt=',sm.ptr.contents.fmt)
    print('err=',err)
    
    bb2=xx_buffer_t(Az.shape[0],len(ab),nrhs=2);
    bb2[0][0,:]=z0;
    bb2[0][1,:]=z0;
    
    tic();st=_klub_solve(h,bb2.ptr,nt);toc('_klub_solve2:')
    print('st=',st)
    tic()
    
    tic();
    for k in range(10):
        st=_klub_solve(h,bb.ptr,nt);
    toc('_klub_solve10:')
    #tic();st=_klub_numeric(h,byref(sm.ptr),4);toc('_klub_numeric:')
    #del K,G,M,sm,d     
    
    '''
    tic();st=_klub_numeric(h,byref(sm.ptr),8);toc('_klub_numericpp:')
    print('st=',st)
    tic();st=_klub_numeric(h,byref(sm.ptr),8);toc('_klub_numericpp:')
    print('st=',st)
    tic();st=_klub_numeric(h,byref(sm.ptr),8);toc('_klub_numericpp:')
    print('st=',st)
    '''
    tic();q=Az*bb[0][0,:];toc('Az*bb[0][0,:]:')
    # 
    _klub_release(h)

