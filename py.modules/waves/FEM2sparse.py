# -*- coding: utf-8 -*-
"""
Created on Sun Jun 26 00:41:31 2022

@author: DICH
"""

import numba  
import numba as nb

import numpy as np

from scipy import sparse as sp
from scipy.sparse import coo_matrix,csc_matrix,csr_matrix
from scipy.sparse.sparsetools import coo_tocsr


@numba.njit
def iperm(i):
    n=len(i);
    ie=np.arange(n,dtype=i.dtype);
    ie[i]=np.arange(n,dtype=i.dtype);    
    return ie;


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
    

    Ap=np.array(indptr,dtype=np.int32,copy=True)
    Aj=np.array(indices,dtype=np.int32,copy=True)
    nnz=Ap[-1];    
    Bx=-np.ones(nnz,dtype=np.int32)    
    
    nnz=_csr_sum_duplicates_mask(n,Ap,Aj,Bx)
    Aj.resize(nnz)
    return (Bx,Aj,Ap);


def icoo_tocsr(n_row,n_col,cols,rows,Ax=None,dtype=np.int32):
    
    Ai,Aj=[np.array(m,dtype=dtype,copy=False) for m in (cols,rows)]
    nnz=len(Ai);
    if Ax is None:
        Ax=np.arange(nnz,dtype=dtype);
    
    Bx=np.empty(nnz,dtype=dtype);
    Bj=np.empty(nnz,dtype=dtype);
    Bp=np.empty(n_row+1,dtype=dtype);
    
    coo_tocsr(n_row,n_col,nnz,
               Ai,Aj,Ax,
               Bp,Bj,Bx);
    
    return (Bx,Bj,Bp)



def iFEM2sparse_coo(trs,lbound=0,dtype=np.int32):
    
    trs=np.array(trs,dtype=dtype,copy=False);
    
    N,nd=trs.shape
    
    lbound=dtype(lbound)    
    if lbound :
        trs=trs-lbound;
        
    Nf=N*nd
    ndndN=nd*nd*N;
    e=np.ones(nd,dtype=dtype);
    row=np.kron(e,trs).reshape((ndndN,));
    col=np.kron(trs,e).reshape((ndndN,));
    idata=np.arange(ndndN,dtype=dtype)
    
    return coo_matrix( (idata,(row,vol)),shape=(Nf,Nf),copy=False);
    
def icoo2pp_csr(A,dtype=np.int32):
    
    (idata,indices,indptr)=icoo_tocsr(A.shape[0],A.shape[1],A.row,A.col,dtype=dtype);
    q=sp.csr_matrix((idata,indices,indptr),shape=A.shape)
    q.sort_indices();     
    p=q.data
    ip=iperm(p)     
    (mask,indices,indptr)=csr_sum_duplicates_mask(q.shape[0],q.indptr,q.indices)
    pp=mask[ip]
    #ip=iperm(p)
    return (pp,indices,indptr)
     
def iFEM2pp_csr(trs,lbound=0,dtype=np.int32):
    icoo=iFEM2sparse_coo(trs,lbound=lbound,dtype=dtype)
    (pp,indices,indptr)=icoo2pp_csr(icoo,dtype=dtype)
    return (pp,indices,indptr);
    
    