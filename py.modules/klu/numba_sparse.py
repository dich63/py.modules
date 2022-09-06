#
import numba  
import numba as nb
import numpy as np


        
    




#@numba.jit("int32(int32,int32,int32,int32[:],int32[:],int32[:],int32[:],int32[:],int32[:])",
@numba.jit(           
           nopython=True,
           parallel=False,
           nogil=True,
           fastmath=True,
           cache=True           
           )
def _icoo_tocsr(n_row,
               n_col,
               nnz,
               Ai,
               Aj,
               Ax,
               Bp,
               Bj,
               Bx):
    #Bp must be zeros
    for n in range(nnz):
        Bp[Ai[n]]+=1
        
    cumsum=0;
    
    for i in range(n_row):
        tmp=Bp[i]
        Bp[i]=cumsum;
        cumsum+=tmp;
    
    Bp[n_row]=nnz;
    
    for n in range(nnz):
        row=Ai[n];
        dest=Bp[row]
        Bj[dest]=Aj[n]
        Bx[dest]=Ax[n]
        Bp[row]+=1;
        
    last=0;
    for i in range(n_row):
        tmp=Bp[i]
        Bp[i]=last;
        last=tmp;
        
    return nnz

    
    
def icoo_tocsr(n_row,n_col,cols,rows):    
    
    Ai,Aj=[np.array(m,dtype=np.int32,copy=False) for m in (cols,rows)]
    nnz=len(Ai);
    Ax=np.arange(nnz,dtype=np.int32);
    Bx=np.empty(nnz,dtype=np.int32);
    Bj=np.empty(nnz,dtype=np.int32);
    Bp=np.zeros(n_row+1,dtype=np.int32);
    _icoo_tocsr(n_row,n_col,nnz,
               Ai,Aj,Ax,
               Bp,Bj,Bx);
    
    return (Bx,Bj,Bp)

def icoo_tocsr_np(n_row,n_col,cols,rows):
    from scipy.sparse.sparsetools import coo_tocsr
    Ai,Aj=[np.array(m,dtype=np.int32,copy=False) for m in (cols,rows)]
    nnz=len(Ai);
    Ax=np.arange(nnz,dtype=np.int32);
    Bx=np.empty(nnz,dtype=np.int32);
    Bj=np.empty(nnz,dtype=np.int32);
    Bp=np.zeros(n_row+1,dtype=np.int32);
    
    coo_tocsr(n_row,n_col,nnz,
               Ai,Aj,Ax,
               Bp,Bj,Bx);
    
    return (Bx,Bj,Bp)
    



    


#@numba.jit("int32(int32,int32[:],int32[:],int32[:])",    
'''
@numba.jit(                      
           nopython=True,
           parallel=False,
           nogil=True,
           fastmath=True,
           cache=True           
           )

'''
@numba.njit
def _csr_sum_duplicates_mask(n_row,Ap,Aj,Bx):
    
    nnz=0;    
    row_end=0;
    for i in range(n_row):
        jj=row_end;
        row_end=Ap[i+1]
        while jj<row_end:
            j=Aj[jj]
            Bx[jj]=nnz;
            jj+=1;
            while jj<row_end and Aj[jj]==j:
                Bx[jj]=nnz;                
                jj+=1;
            Aj[nnz]=j
            nnz+=1
            
        Ap[i+1]=nnz
        
    return nnz;

def csr_sum_duplicates_mask(n,indptr,indices):
    '''
    @numba.jit(                      
           nopython=True,
           parallel=False,
           nogil=True,
           fastmath=True,
           cache=True           
           )
    def _csr_sum_duplicates_mask(n_row,Ap,Aj,Bx):
        
        nnz=0;    
        row_end=0;
        for i in range(n_row):
            jj=row_end;
            row_end=Ap[i+1]
            while jj<row_end:
                j=Aj[jj]
                Bx[jj]=nnz;
                jj+=1;
                while jj<row_end and Aj[jj]==j:
                    Bx[jj]=nnz;                
                    jj+=1;
                Aj[nnz]=j
                nnz+=1
                
            Ap[i+1]=nnz
            
        return nnz;
    '''
    #_csr_sum_duplicates_mask=numba.njit()(__csr_sum_duplicates_mask)
    Ap=np.array(indptr,dtype=np.int32,copy=True)
    Aj=np.array(indices,dtype=np.int32,copy=True)
    nnz=Ap[-1];    
    Bx=-np.ones(nnz,dtype=np.int32)    
    
    nnz=_csr_sum_duplicates_mask(n,Ap,Aj,Bx)
    Aj.resize(nnz)
    return (Bx,Aj,Ap);
        

@numba.jit("int32(int32,complex128[:],complex128[:],int32[:])"
           ,nopython=True,
           parallel=False,
           nogil=True,
           fastmath=True,
           cache=True                                 
           )
def _permute_sum(n,s,d,p):
    for i in range(n):
        d[p[i]]+=s[i];
    return n;

prange=numba.prange

#@numba.jit("int32(int32,int32,complex128[:,:],complex128[:,:],int32[:])",
@numba.jit(
           nopython=True,
           parallel=True,
           nogil=True,
           fastmath=True           
           )
def _permute_sum2(m,n,s,d,ip):
    
    for i in prange(m):
        for k in range(n):
            d[i,ip[k]]+=s[i,k];
    return n;



'''
@numba.jit("int32(int32,int32,complex128[:,:],complex128[:,:],int32[:])"
           ,nopython=True,
           parallel=True,
           nogil=True,
           fastmath=True           
           )
def _permute_sum2(M,N,s,d,p):    
    
    for m in numba.prange(M):        
        for n in range(N):
            d[m,p[n]]=s[m,n];
            
    return n;
'''


def permute_sum2(s,d,p):
    m,n=s.shape;
    return _permute_sum2(m,n,s,d,p);

#@numba.jit
def permute_sum2_np(s,d,ip):  
    m,n=s.shape;
    for k in range(n):
        d[:,ip[k]]+=s[:,k];
    
    return n;

def permute_sum2_np2(s,d,ip):  
    m,n=s.shape;
    dt=d.T
    st=s.T
    for k in range(n):
        dt[ip[k],:]+=st[k,:];
    d[:,:]=dt.T
    return n;


def permute_sum(s,d,p):
    nz=len(s);
    return _permute_sum(nz,s,d,p)


#@numba.jit("int32[:](int32[:])",nopython=True,nogil=True,cache=True)
@numba.njit
def iperm(i):
    n=len(i);
    ie=np.arange(n,dtype=i.dtype);
    ie[i]=np.arange(n,dtype=i.dtype);    
    return ie;


def iperms(i):
    n=len(i);
    ie=np.arange(n,dtype=np.int32);
    ie[i]=np.arange(n);
    return ie;
    
'''    
    


@numba.jit("int32(int32,complex128[:],complex128[:],int32[:])"
           ,nopython=True,
           parallel=True,
           nogil=True,
           fastmath=True           
           )
def permute(n,s,d,ip):
    for i in numba.prange(n):
        d[i]=s[ip[i]];
    return n;


@numba.jit("int32(int32,int32,complex128[:,:],complex128[:,:],int32[:])"
           ,nopython=True,
           parallel=True,
           nogil=True,
           fastmath=True           
           )
def permute2(m,n,s,d,ip):
    for i in numba.prange(n):
        for k in range(m):
            d[k,i]=s[k,ip[i]];
    return n;

@numba.jit("int32(int32,int32,complex128[:,:],complex128[:,:],int32[:,:])"
           ,nopython=True,
           parallel=True,
           nogil=True,
           fastmath=True           
           )
def permute3(m,n,s,d,ip):
    for i in numba.prange(n):
        for k in range(m):
            d[k,i]=s[k,ip[k,i]];
    return n;

@numba.jit("int32(int32,int32,complex128[:,:],complex128[:,:],int32[:,:],float64[:,:])"
           ,nopython=True,
           parallel=True,
           nogil=True,
           fastmath=True           
           )
def permute3_rescale(m,n,s,d,ip,rs):
    for i in numba.prange(n):
        for k in range(m):
            j=ip[k,i]
            d[k,i]=s[k,j]/rs[k,j];
    return n;
'''


def is_sparse_pattern_equ(*lsp): 
    
    def _sp_info(s): 
        fmt=s.format;
        if fmt!='coo':
            return (fmt,s.nnz,s.indptr,s.indices)
        else:
            return (fmt,s.nnz,s.col,s.row)
    

    
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


import scipy
import scipy.sparse as sp

def unicoo(AA):
    
    
    
    nnzs=[m.nnz for m in AA]
    nnz=np.sum(nnzs);  
    
    
    zdatas=[ np.zeros_like(m.data) for m in AA ]
    
    rows=np.concatenate([ np.array(m.row) for m in AA ])    
    cols=np.concatenate([ np.array(m.col) for m in AA ])
    rAA=[]
    for m in range(len(AA)):
        d=zdatas.copy();
        d[m][:]=AA[m].data;
        d=np.concatenate(d)
        rAA+=[sp.coo_matrix((d,(rows,cols)),shape=AA[m].shape)]
    
    return rAA;
    

class sp_batch_maitix(object):
    def __init__(self,AA):
        if AA is not None:
            self.feq=feq=is_sparse_pattern_equ(*AA);
            if feq:
                pass
                
            
            #indices,indptr=
            #AA=[sp.coo_matrix(m,copy=False) for m in AA];
            
            



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
    
    fn=r'O:\__ss\matrix\sfFEM1024k.json'
    #
    fn=r'O:\__ss\matrix\sfFEM128k.json'
    #fn=r'O:\__ss\matrix\sfFEM24k.json'
    #fn=r'O:\__ss\matrix\KGM.json'
    #fn=r'O:\__ss\matrix\sfFEM0k.json'
    d=decode(fn,1)
    
    K,G,M=[m for m in [d.K,d.G,d.M]]
    K,G,M=[sp.csr_matrix(m,dtype=complex) for m in [d.K,d.G,d.M]]
    
    K=K+sp.tril(K);
    
    K,G,M=[m.tocoo() for m in [K,G,M]]
    
    
    
    K0,G0,M0=K,G,M
    tic(); feq= is_sparse_pattern_equ(K,G,M);mtoc('is_sparse_pattern_equ=%d:'%feq)
    if 0 and feq:
        uK,uG,uM=K,G,M
    else:
        tic(); uK,uG,uM=unicoo([K,G,M]);mtoc('unicoo')
    
    K,G,M=uK,uG,uM
    
    kkk=K.tocsr();
    
    tic();t=icoo_tocsr(K.shape[0],K.shape[0],K.row,K.col);mtoc('icc')
    n=K.shape[0]
    idata=np.arange(K.nnz,dtype=np.int32)
    tic();
    iK=sp.coo_matrix((idata,(K.row,K.col)),shape=K.shape)
    mtoc('iK')
    tic();
    q=sp.csr_matrix(t,shape=K.shape);
    mtoc('crs')    
    q.sort_indices();
    mtoc('crs+sort')
    kK=iK.tocsr()
    tic()
    qq=q.copy()
    qq.sum_duplicates()
    mtoc('sum_np')    
    
    tic()
    (Bx,Aj,Ap)=csr_sum_duplicates_mask(n,q.indptr,q.indices)
    mtoc('sum')    
    nnz=Ap[-1]
    nnzs=K.nnz
    
    cdata=K.data+0j
    p=q.data
    ip=iperm(p)
    cdp=cdata[p]
    cD=sp.csr_matrix(K,dtype=complex)
    d=np.zeros(nnz,dtype=complex);permute_sum(cdp,d,Bx);err=norm(cD.data-d);err
    print('err=',err)
    Bxp=Bx[ip]
    d=np.zeros(nnz,dtype=complex);permute_sum(cdata,d,Bxp)
    
    d=np.zeros(nnz,dtype=complex);
    tic();permute_sum(cdata,d,Bxp);t1=mtoc('permute_sum s')
    err=norm(cD.data-d);err
    print('err=',err)
    
    mp=16
    
    cds=np.empty((mp,cdata.size),dtype=complex)
    cds[:]=cdata;
    ds=np.zeros((mp,nnz),dtype=complex);permute_sum2(cds,ds,Bxp)
    ds=np.zeros((mp,nnz),dtype=complex)
    tic();permute_sum2(cds,ds,Bxp);t2=mtoc('permute_sum2 pp')
    
    print("perf=",mp*t1/t2)

    

    

    


