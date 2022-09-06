#

from ctypes import *
from utils.p23 import *
from utils.c_pointers import *


class sp_t(Structure):
    _pack_=8
    _fields_ = [
    ("indptr",c_void_p),
    ("indices",c_void_p),
    ("ppdata",c_void_p),  
    ("tclass",c_uint64),             
    ("count",c_int64),    
    ("n",c_int64),
    ("m",c_int64),
    ("nnz",c_int64),    
    ("fmt",c_char*8)       
    ];

p_sp_t  = POINTER(sp_t)
null_sp=p_sp_t()
 
class rh_t(Structure):    
    _pack_=8
    _fields_ = [    
    ("xx",c_void_p),  
    ("n",c_int64),        
    ("nrhs",c_int64)        
    ]
    
    

p_rh_t  = POINTER(rh_t)
null_rh=p_rh_t()

class context_t(Structure):
    _fields_ = []    

h_context_t  = POINTER(context_t)
ph_context_t  = POINTER(h_context_t)


class klu_common_t(Structure):
    _pack_=8
    _fields_ = [
    ("tol",c_double),
    ("memgrow",c_double),
    ("initmem_amd",c_double),
    ("initmem",c_double),
    ("maxwork",c_double),

    ("btf",c_int64),
    ("ordering",c_int64),  
    ("scale",c_int64),
    #("nz",c_longlong),

    ("user_order",c_void_p),  
    ("user_data",c_void_p),  

    ("halt_if_singular",c_int64),

    ("status",c_int64),  
    ("nrealloc",c_int64),

    ("structural_rank",c_int64),  

    ("numerical_rank",c_int64),

    ("singular_col",c_int64),

    ("noffdiag",c_int64),

    ("flops",c_double),
    ("rcond",c_double),
    ("condest",c_double),
    ("rgrowth",c_double),
    ("work",c_double),

    ("memusage",c_longlong),
    ("mempeak",c_longlong)

    ];


p_klu_common_t  = POINTER(klu_common_t)
null_common=p_klu_common_t()


_klu_path=r"O:\__ss2022\build\x64\Release\klu_lib.dll"    
#_klu_path=r"O:\__ss2022\build\x64\debug\klu_lib.dll"    
_klu=cdll.LoadLibrary(_klu_path)


_klu_context_release=_klu.klu_context_release;
_klu_context_release.rstype=c_long
_klu_context_release.argtypes = (h_context_t,)

_context_release=_klu.klu_context_release;
_context_release.rstype=c_long
_context_release.argtypes = (c_void_p,)


_context_refcount=_klu.klu_context_refcount;
_context_refcount.rstype=c_long
_context_refcount.argtypes = (c_void_p,)


_klu_context_symbolic=_klu.klu_context_symbolic;
_klu_context_symbolic.rstype=c_long
_klu_context_symbolic.argtypes = (p_sp_t,p_klu_common_t,ph_context_t)

_klu_context_numeric=_klu.klu_context_numeric;
_klu_context_numeric.rstype=c_long
#_klu_context_numeric.argtypes = (p_sp_t,c_int64,h_context_t,p_klu_common_t,ph_context_t)
_klu_context_numeric.argtypes = (c_void_p,c_int64,c_void_p,c_void_p,c_void_p)



_klu_context_solve=_klu.klu_context_solve;
_klu_context_solve.rstype=c_long
_klu_context_solve.argtypes = (h_context_t,c_void_p,c_int64,p_klu_common_t)

_klu_common_defaults=_klu.klu_common_defaults;
_klu_common_defaults.rstype=c_long
_klu_common_defaults.argtypes = (p_klu_common_t,)


_klu_get_common=_klu.klu_get_common;
_klu_get_common.rstype=c_long
_klu_get_common.argtypes = (h_context_t,p_klu_common_t)

'''
klu_sleep=_klu.klu_sleep;
klu_sleep.rstype=c_long
klu_sleep.argtypes = (c_void_p,)
'''


import numpy as np
import scipy
import scipy.sparse as sp

def csc_2_sp_t(A,dtype=np.complex128):   
    
    A=sp.csc_matrix(A,dtype=dtype);
    
    indptr =  A.indptr.ctypes.data_as(c_void_p)
    indices = A.indices.ctypes.data_as(c_void_p)    
      
    n,m=A.shape;
    
    ppdata=ptr_ptr_array([ A.data]); 
    
    spm=sp_t(n=n,m=m,count=1,tclass=0x0109,
              ppdata=ppdata.ptr,
              indptr=indptr,
              indices=indices,fmt=to_bytes(A.format),nnz=-1)
    
    return link_ptr_t(pointer(spm),spm,ppdata,indptr,indices,A);





def _common_update(**opts):
    common=klu_common_t();
    _klu_common_defaults(byref(common));
    c_update(common,ordering=0,scale=-1,btf=0)
    return c_update(common,**opts);
    
class klu_t(object):
    
    def __init__(self,A,**opts): 
        
        self.spm=spm=csc_2_sp_t(A);
        self.hsymbolic=hsymbolic=h_context_t()
        self.hfactor=h_context_t();
        
        self._phase=-1;            
        self.common=common=_common_update(**opts);
        
        
        status=_klu_context_symbolic(spm.ptr,byref(common),byref(hsymbolic))
        if status==0:
            self._phase=1;
        self._status=status;
        
    def __del__(self):
        _klu_context_release(self.hsymbolic)
        _klu_context_release(self.hfactor)
        '''
        try:
            print('del _context')             
            _context_release(self.hsymbolic)
            _context_release(self.hfactor)
        except:
            print('+++_context_release+++') 
            pass
        '''
        
    @property
    def status(self):
        return self._status;    
    @property
    def phase(self):
        return self._phase;
    @property
    def opts(self):
        return c_2_dict(self.common);   
    
    @property
    def info(self,fnum=True):
        common=klu_common_t();
        h=self.hfactor if fnum else self.hsymbolic;
        _klu_get_common(h,byref(common));
        return c_2_dict(common);   

        
        
    def factorize(self):   
        if self._status==0:
            
            if  self._phase==2:
                return 0;
            
            self._status=_klu_context_numeric(
                self.spm.ptr,0,
                self.hsymbolic,null_common,
                byref(self.hfactor)
                )
            if self._status==0:
                self._phase=2              
            
        return self._status;
    
    def __call__(self,B,nrhs=1):
        
        if self._status==0:   
            if self.phase<2:
                if self.factorize():
                    return self._status;
            bp=ptr_array(B);
            st=_klu_context_solve(self.hfactor,bp.ptr,nrhs,null_common);
            
        return self._status
       
    
        



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
    
    norm=np.linalg.norm
    normm= lambda x:  norm(x,ord=np.inf)
    
    A=np.array([[1,1j],[2,-2]])
    x=np.array([1,11j],dtype=complex)
    
    b=1*x;
    
    
    
    spm=csc_2_sp_t(A);
    
    hsymbolic=h_context_t()
    hfactor=h_context_t()
    print('hex(addressof(spm.ptr))=',hex(addressof(spm.ptr)))
    print('spm.ptr',spm.ptr)
    print('hex(addressof(spm.ptr.contents))',hex(addressof(spm.ptr.contents)))
    print('hex(addressof(hsymbolic))',hex(addressof(hsymbolic)))
    print('pid=',os.getpid())
    st=_klu_context_symbolic(spm.ptr,null_common,byref(hsymbolic))
    print('klu_context_symbolic=',st);
    #st=_klu_context_numeric(spm.ptr,0,hsymbolic,null_common,byref(hfactor))
    st=_klu_context_numeric(spm.iptr,0,hsymbolic,None,byref(hfactor))
    print('klu_context_numeric=',st);
    
    
    
    
    #raise SystemExit()
    common=klu_common_t() 
    bp=ptr_array(b);
    st=_klu_context_solve(hfactor,bp.ptr,1,byref(common))
    print('klu_context_solve=',st);
    err=norm(A@b-x)
    print('err=',err)
    '''
    
    klu=klu_t(A)
    print('klu_context_symbolic=',klu.status);
    klu.factorize();    print('klu.factorize=',klu.status);
    print('pid=',os.getpid())
    klu(b)
    print('klu.solve=',klu.status);
    err=norm(A@b-x)
    print('err=',err)
    print('common=',klu.info)
    
    #del klu
    '''
