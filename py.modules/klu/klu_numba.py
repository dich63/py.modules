#

import numpy as np
import numba  
import numba as nb
from klu.iklu import *
from klu.numba_permutations import *
import functools
#st=_iklu_numeric(n,indptr.ptr,indices.ptr,pdata.ptr,tclass,hsymbolic.iptr,None,hfactor.ibyref)

#@numba.jit("int32[:](int64,int64,int64,int64,int64[:],int64,int64,int64)"
@numba.jit("int32[:](int64,int64,int64,int64,int64[:],int64,int64,int64,int32[:])"
           ,nopython=True,
           parallel=True,
           nogil=True
           )
def _iklu_numeric_batch(nbatch,n,indptr,indices,pdata,tclass,hsymbolic,ppfactors,states):    
    
    
    for i in numba.prange(nbatch):        
        states[i]=iklu_numeric(n,indptr,indices,pdata[i],tclass,hsymbolic,0,ppfactors+8*i)
        
    return states 

'''
@numba.jit("int32[:](int64,int64,int64,int64,int64[:],int64,int64,int64,int32[:])"
           ,nopython=True,
           parallel=False,
           nogil=True,
           fastmath=True           
           )
'''
def _iklu_numeric_batch_s(nbatch,n,indptr,indices,pdata,tclass,hsymbolic,ppfactors,states):        
    
    for i in range(nbatch):        
        states[i]=iklu_numeric(n,indptr,indices,int(pdata[i]),tclass,hsymbolic,0,ppfactors+8*i)
        
    return states 

    

@numba.jit("int32[:](int64,int64[:],int64[:],int32[:],int64)"
           ,nopython=True,
           parallel=True,
           nogil=True
           )
def _iklu_solve_ts_batch(nbatch,factors,xx,states,pcommon):     
    
    for i in numba.prange(nbatch):        
        states[i]=iklu_solve_ts(factors[i],xx[i],1,pcommon)
        
    return states 

def _iklu_solve_ts_batch_s(nbatch,factors,xx,states,pcommon):     
    
    for i in range(nbatch):        
        states[i]=iklu_solve(int(factors[i]),int(xx[i]),1,int(pcommon))
        
    return states 


@numba.jit("int32[:](int64,int64[:],int64[:],int32[:],int64)"
           ,nopython=True,
           parallel=True,
           nogil=True
           )
def _iklu_tsolve_ts_batch(nbatch,factors,xx,states,pcommon):     
    
    for i in numba.prange(nbatch):        
        states[i]=iklu_tsolve_ts(factors[i],xx[i],1,pcommon)
        
    return states 

def _iklu_tsolve_ts_batch_s(nbatch,factors,xx,states,pcommon):     
    
    for i in range(nbatch):        
        states[i]=iklu_tsolve(int(factors[i]),int(xx[i]),1,int(pcommon))
        
    return states 



def to_handles(xx=None,nbatch=0,release=None,addref=None):
    if type(xx)== handle_array_t:
        return xx;
    else:
        return handle_array_t(xx,nbatch,release,addref);


    
    
def iklu_numeric_batch(n,indptr,indices,datas,tclass,hsymbolic,parallel=True):
    
    datas=to_handles(datas)
    nbatch=len(datas);
    pdata=datas.addresses 
    
    hfactors=to_handles(nbatch=nbatch,release=iklu_release,addref=iklu_addref);
    hfactors.states=states=np.full(nbatch,-1,dtype=np.int32);
    ppfactors=hfactors.ptr
    p=indptr.ctypes.data
    i=indices.ctypes.data
    
    if parallel:
        _iklu_numeric_batch(nbatch,n,p,i,pdata,tclass,hsymbolic.iptr,ppfactors,states)
    else:
        _iklu_numeric_batch_s(nbatch,n,p,i,pdata,tclass,hsymbolic.iptr,ppfactors,states)
        
    return hfactors

def iklu_solve_ts_batch(factors,xx,states=None,common=None,parallel=True):   
    
    nbatch=len(factors);
    pcommon=0 if common is None else addressof(common)     
    if states is None:
        states=np.full(nbatch,-1,dtype=np.int32);        
    
    ifactors=factors.addresses    
    ixx=to_handles(xx).addresses;
    if parallel:
        _iklu_solve_ts_batch(nbatch,ifactors,ixx,states,pcommon);
    else:
        _iklu_solve_ts_batch_s(nbatch,ifactors,ixx,states,pcommon);
    
    
    return states;


def iklu_tsolve_ts_batch(factors,xx=None,states=None,common=None,parallel=True,ixx=None):   
    
    nbatch=len(factors);
    pcommon=0 if common is None else addressof(common)     
    if states is None:
        states=np.full(nbatch,-1,dtype=np.int32);        
    
    ifactors=factors.addresses    
    if ixx is None:
        ixx=to_handles(xx).addresses;
    if parallel:
        _iklu_tsolve_ts_batch(nbatch,ifactors,ixx,states,pcommon);
    else:
        _iklu_tsolve_ts_batch_s(nbatch,ifactors,ixx,states,pcommon);
    
    
    return states;




#st=iklu_symbolic(n,indptr.ptr,indices.ptr,tclass,addressof(common),hsymbolic.ibyref)

def iklu_symbolic_batch(n,indptr,indices,tclass,common=None):
    
    p=indptr.ctypes.data
    i=indices.ctypes.data       
    
    hsymbolic=smart_ptr_t(h_context_t())
    pcommon=0 if common is None else addressof(common) 
    hsymbolic.state=iklu_symbolic(n,p,i,tclass,pcommon,hsymbolic.ibyref)
    return hsymbolic

def iklu_numeric_permuts(factors):
    M=len(factors);
    PP,RR=[],[]
    fNone=False;
    for m in range(M):
        err,p,r=get_number_PRS(factors[m]);
        PP+=[p]
        RR+=[r]
        
    RR=[r for r in filter(lambda x: not x is None,RR)];
    lR=len(RR)
    if lR==M:
        return [np.array(PP),np.array(RR)]
    elif lR==0:
        return [np.array(PP),None]
    else:
        raise Exception('factors contain not equal permuts')

        
def iklu_tsolve_ts_batch_pp(factors,xx,xxbuf,QPnRs,states=None,xxout=None):   
    
    nbatch=len(factors);
    Q,Pn,Rs=QPnRs
    n=len(Q)
     
    if states is None:
        states=np.full(nbatch,-1,dtype=np.int32);        
    
    common=klu_common_t(user_data=0xbabaeb);    
    pcommon=addressof(common)

    ifactors=factors.addresses    
    
    ixx=to_handles(xxbuf).addresses;
    permute2(nbatch,n,xx,xxbuf,Q)
    
    _iklu_tsolve_ts_batch(nbatch,ifactors,ixx,states,pcommon);    
    
    if not (xxout is None):
        xx=xxout
    
    if Rs is None:
        permute3(nbatch,n,xxbuf,xx,Pn)
    else:
        permute3_rescale(nbatch,n,xxbuf,xx,Pn,Rs)
    
    
    return states;
    
    


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
    
    
    common=common_update();
    
    A=np.array([[1,1j],[2,1j]])
    B=np.array([[1,1j],[2,-2j]])
    x=np.array([1,11j],dtype=complex)
    
    b=1*x;
    
    
    sA=sp.csc_matrix(A,dtype=complex)
    sB=sp.csc_matrix(B,dtype=complex)
    
    indptr=ptr_array(sA.indptr)
    indices=ptr_array(sA.indices)
    pdata=ptr_array(sA.data)
    
    bp=ptr_array(b);
    #spm=csc_2_sp_t(A);
    
    
    
    hsymbolic=smart_ptr_t(h_context_t())
    hfactor=smart_ptr_t(h_context_t())
    
    #hsymbolic=link_ptr_t(h_context_t())
    
    n=A.shape[0];
    tclass=np.int64(0x0109)
    
    
    #st=iklu_symbolic(n,indptr.ptr,indices.ptr,tclass,addressof(common),hsymbolic.ibyref)
    #print('iklu_symbolic=',st);
    hsymbolic=iklu_symbolic_batch(n,sA.indptr,sA.indices,tclass,common)
    print('iklu_symbolic_batch=',hsymbolic.state);
    
    #datas=ptr_ptr_array([sA.data,sA.data])
    datas=[m.data.ctypes.data for m in [sA,sB,sA,sB]]
    
    datas=[m.data for m in [sA,sB,sA,sB]]
    bb=np.array([b,b,b,b])
    hbb=to_handles(bb)
    
    dataso=np.empty(4,dtype=np.object);
    for m in range(4):
        dataso[m]=datas[m]
    
    print('pid=',os.getpid())
    factors=iklu_numeric_batch(n,sA.indptr,sA.indices,datas,tclass,hsymbolic,parallel=0)   
    
    
    print('iklu_numeric=',factors.states);
    
    st=iklu_solve_ts_batch(factors,bb,parallel=0)   
    print('iklu_solve=',st)
    
    err=norm(A@bb[0]-x);    print('err=',err)
    raise SystemExit()
    
    
    st=iklu_solve_ts(factors[2],bp.ptr,1,addressof(common))
    print('iklu_solve=',st)
    print('common=',str(c_2_dict(common)));
    
    
    err=norm(A@b-x)
    print('err=',err)
    
    raise SystemExit()
    st=iklu_numeric(n,indptr.ptr,indices.ptr,pdata.ptr,tclass,hsymbolic.iptr,0,hfactor.ibyref)
    print('iklu_numeric=',st);
    st=iklu_solve_ts(hfactor.iptr,bp.ptr,1,addressof(common))
    print('iklu_solve=',st)
    print('common=',str(c_2_dict(common)));
    
    
    err=norm(A@b-x)
    print('err=',err)

        
    