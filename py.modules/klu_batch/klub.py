# -*- coding: utf-8 -*-
"""
Created on Mon Apr 25 03:53:05 2022

@author: wwww
"""
from ctypes import *
from utils.c_structs import *

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
    
    
    
class spm_t(Structure):
    _pack_=8
    _fields_ = [
    ("indptr",POINTER(c_int32)),
    ("indices",POINTER(c_int32)),
    ("pdata",c_void_p),  
    ("n",c_long),
    ("nbatch",c_long),
    ("nnz",c_long),
    ("fmt",c_char*4)
    ];


    
class rhs_t(Structure):
    _pack_=8
    _fields_ = [    
    ("xx",c_void_p),  
    ("nrhs",c_long),        
    ]


p_klu_batch_t  = POINTER(klu_batch_t)
pp_klu_batch_t  = POINTER(p_klu_batch_t)

p_spm_t  = POINTER(spm_t)
p_rhs_t  = POINTER(rhs_t)

null_spm=p_spm_t()
    
_klu_path=r"O:\__ss2022\build\x64\Release\klu_batch.dll"    
#_klu_path=r"O:\__ss2022\build\x64\debug\klu_batch.dll"    
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





import numpy as np
import scipy as sp


def _sp_info(s):
    
    fmt=s.format;
    if fmt!='coo':
        return (fmt,s.nnz,s.indptr,s.indices)
    else:
        return (fmt,s.nnz,s.col,s.row)
        
    
def is_sparse_pattern_equ(*lsp):
    
    if len(lsp):    
        
        fmt0,nnz0,p0,i0=_sp_info(lsp[0])        
        
        for s in lsp[1:]:
            
            fmt,nnz,p,i=_sp_info(s)
            
            if not (  (fmt0==fmt) and (nnz0==nnz) 
                    and  np.array_equal(i0, i) 
                    and np.array_equal(p0,p)
                    ):
                return False;
        return True;
    else:
        return False;
            


class link_ptr_t(object):
    def __init__(self,p,*link):
        self._ptr,self._link=p,link
    @property
    def ptr(self):
        return self._ptr;

def ptr_ptr_array(pp,offset=0):    
    L=[p.ctypes.data for p in pp]
    for k in range(offset):
        L.insert(0,0);
    ip=np.array(L,dtype=np.uint64);
        
    return link_ptr_t(ip.ctypes.data,ip,pp);
    
    
def make_spm(*lsp):
    if is_sparse_pattern_equ(*lsp):
        nbatch=len(lsp);        
        b=lsp[0]
        n=b.shape[0]
        nnz=b.nnz
        
        
        fmt,nnz,p,i = _sp_info(b)
        
        indptr =  p.ctypes.data_as(POINTER(c_int32))
        indices = i.ctypes.data_as(POINTER(c_int32))        
        
        fmt=bytes(fmt,'utf8');
        
        pdata=ptr_ptr_array([ m.data for m in lsp ]); 
        
        '''
        fmt=bytes(b.format,'utf8');
        
        
        pdata=ptr_ptr_array([ m.data for m in lsp ]); 
        
        
        if fmt!=b'coo':
            indptr=b.indptr.ctypes.data_as(POINTER(c_int32))
            indices=b.indices.ctypes.data_as(POINTER(c_int32))        
        else:
            indptr=b.col.ctypes.data_as(POINTER(c_int32))
            indices=b.row.ctypes.data_as(POINTER(c_int32))        
        
        '''
        spm=spm_t(n=n,nbatch=nbatch,
                  pdata=pdata.ptr,
                  indptr=indptr,
                  indices=indices,fmt=fmt,nnz=nnz)
        
        return link_ptr_t(pointer(spm),spm,pdata,indptr,indices,lsp);
        
    else:
        raise Exception('not is_sparse_pattern_equ ')
        
        

    
class xx_buffer_t(object):   
    
    def __init__(self,n,nbatch,nrhs=1,dtype=np.complex128):
        
        #xxl=[np.ascontiguousarray(np.zeros((nrsh,n),dtype=dtype)) for k in range(nbatch)];
        
        xx=np.empty(nbatch,dtype=np.object);     
        for k in range(nbatch):
            xx[k]=np.ascontiguousarray(np.zeros((nrhs,n),dtype=dtype))       
        self._xx=xx;        
        self._ppxx=ppxx=ptr_ptr_array(xx);        
        self._rhs=rhs_t(xx=ppxx.ptr,nrhs=nrhs)
        
    @property
    def ptr(self):
        return pointer(self._rhs);
        
    def __getitem__(self,i):        
        return self._xx[i];
    
    def __setitem__(self,i,value):
        self._xx[i][:]=value;
        
    def __iter__(self):    
        return iter(self._xx)
    
    def __len__(self):
        return len(self._xx)
    
         
        
    
    

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
    d=decode(fn,1)
    K,G,M=[m.tocsr() for m in [d.K,d.G,d.M]]
    
    zz=[7+17j,2+1j,2+4j,2-4j,1+1j,2+1j,2+4j,2-4j]
    zz=[7+17j,2+1j,2+4j,2-4j]
    #zz=[1+1j,2+1j]
    #    zz=[7+17j,]
    
    N=M.shape[0]
    di=np.ones(N-1)
    
    G=sp.tril(G,format='csc')
    ab=[K+z*G+z**2*M for z in zz ];
    
    
    ab=[m.tocsr() for m in ab ];
    
    
    
    #ab=[m.T+m for m in ab ];
    
    
    ab=[m.tocsr() for m in ab ];
    Az=ab[0]
    
    ab=[m.tocoo() for m in ab ];
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
    nt=16
    #c_update(h.contents.Common,ordering=0,scale=-1,btf=0)
    c_update(h.contents.Common,ordering=0,scale=-1,btf=0)
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
    
    
    tic();st=_klub_solve(h,bb.ptr,nt);toc('_klub_solve:')
    print('st=',st)
    
    err=norm(Az*bb[0][0,:]-z0)/norm(Az*bb[0][0,:])
    print('fmt=',sm.ptr.contents.fmt)
    print('err=',err)
    
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
    
    #  
    _klub_release(h)

