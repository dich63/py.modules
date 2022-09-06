#

from ctypes import *
from utils.p23 import *
from utils.c_pointers import *

    
    
    
class spm_t(Structure):
    _pack_=8
    _fields_ = [
    ("indptr",POINTER(c_int32)),
    ("indices",POINTER(c_int32)),
    ("pdata",c_void_p),  
    ("tdata",c_long),  
    ("n",c_long),
    ("m",c_long),
    ("nbatch",c_long),
    ("nnz",c_long),    
    ("fmt",c_char*4)
    ];


    




p_spm_t  = POINTER(spm_t)
null_spm=p_spm_t()


def _sp_info(s):
    
    fmt=s.format;
    if fmt!='coo':
        return (fmt,s.nnz,s.indptr,s.indices)
    else:
        return (fmt,s.nnz,s.col,s.row)
        
    
def is_sparse_pattern_equ(*lsp):
    
    lsp=list(filter(lambda s : not s is None,lsp));
    
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
            

def make_spm_csr_empty(n,nnz,nbatch,m=-1,dtype=np.complex128):
    
    if m<0:
        m=n;
        
    cc=cache_t()
    
    p = cc(np.empty(n+1,dtype=np.int32));
    i = cc(np.empty(nnz,dtype=np.int32));
    indptr =  cc(p.ctypes.data_as(POINTER(c_int32)))
    indices = cc(i.ctypes.data_as(POINTER(c_int32))) 
    
    datas=np.empty(nbatch,dtype=object); 
    for k in range(nbatch):
        datas[k]=np.empty((nnz,),dtype=dtype);
    
    pdata=cc(ptr_ptr_array(datas))
    
    
    
    spm=spm_t(n=n,m=m,nbatch=nbatch,
                  pdata=pdata.ptr,
                  indptr=indptr,
                  indices=indices,fmt=b'csr',nnz=nnz)
    
    return link_ptr_t(pointer(spm),spm,cc);


def new_spm_like(spm,nbatch=-1,datas=None,dtype=np.complex128):    
    
    
    sm=spm.ptr.contents
    if datas is None:
        datas=np.zeros((nbatch,sm.nnz),dtype=dtype);
        #datas=np.nan*np.ones((nbatch,sm.nnz),dtype=dtype);
    else:
        nbatch,n=datas.shape;
        
    ppdata=ptr_ptr_array(datas);
    
    smn=c_copy(sm,nbatch=nbatch,pdata=ppdata.ptr);
    
    return link_ptr_t(pointer(smn),smn,ppdata)
     
    
    
def make_spm(*lsp,dtype=np.complex128,fnocheck=False):
    if fnocheck or is_sparse_pattern_equ(*lsp):
        nbatch=len(lsp);        
        b=lsp[0]
        if b.dtype!=dtype:
            raise Exception('Matrix type not is '+str(dtype))
            
        n,m=b.shape
        
        nnz=b.nnz
        
        
        fmt,nnz,p,i = _sp_info(b)
        
        indptr =  p.ctypes.data_as(POINTER(c_int32))
        indices = i.ctypes.data_as(POINTER(c_int32))        
        
        fmt=bytes(fmt,'utf8');
        
        pdata=ptr_ptr_array([ m.data for m in lsp ]); 
        
        spm=spm_t(n=n,m=m,nbatch=nbatch,
                  pdata=pdata.ptr,
                  indptr=indptr,
                  indices=indices,fmt=fmt,nnz=nnz)
        
        return link_ptr_t(pointer(spm),spm,pdata,indptr,indices,lsp);
        
    else:
        raise Exception('not is_sparse_pattern_equ ')
        

def spm_triplet(sm,nb=0,dtype=complex):
    
    n=sm.n;
    nbatch=sm.nbatch;
    
    if nb<0:
        nb=nbatch+nb;
    
    if sm.fmt==b'coo':
        nnz=sm.nnz
        buf=(c_long*nnz).from_address(addressof(sm.indptr.contents))
        indptr=np.frombuffer(buf,dtype=np.int32);
    else:
        buf=(c_long*(n+1)).from_address(addressof(sm.indptr.contents))
        indptr=np.frombuffer(buf,dtype=np.int32);
        nnz=indptr[n];
    
    
    buf=(c_long*nnz).from_address(addressof(sm.indices.contents))
    
    indices=np.frombuffer(buf,dtype=np.int32);
    
    szof=ptr_sizeof()
    buf=(c_byte*(szof*nbatch)).from_address(sm.pdata)
    
    ipdata=np.frombuffer(buf,dtype=iptr_type());        
    
    szof=np.empty(1,dtype=dtype).nbytes
    buf=(c_byte*(szof*nnz)).from_address(int(ipdata[nb]))
    data= np.frombuffer(buf,dtype=dtype);         
    
    return indptr,indices,data
        
        

    
        
        
    
