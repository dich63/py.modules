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


@numba.jit(
           nopython=True,
           parallel=True,
           nogil=True
           )

def _permute_sum2(m,n,s,d,ip):
    
    for i in nb.prange(m):
        for k in range(n):               
            d[i,ip[k]]+=s[i,k];
    return n;


@numba.njit
def iperm(i):
    n=len(i);
    ie=np.arange(n,dtype=i.dtype);
    ie[i]=np.arange(n,dtype=i.dtype);    
    return ie;


@numba.njit(cache=True)
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


def permute_sum2(s,d,p):
    m,n=len(s),len(p);
    return _permute_sum2(m,n,s,d,p);


            
def csr_sum_duplicates_mask(n,indptr,indices):
    

    Ap=np.array(indptr,dtype=np.int32,copy=True)
    Aj=np.array(indices,dtype=np.int32,copy=True)
    nnz=Ap[-1];    
    Bx=-np.ones(nnz,dtype=np.int32)    
    
    nnz=_csr_sum_duplicates_mask(n,Ap,Aj,Bx)
    Aj.resize(nnz,refcheck=False)
    #Aj=np.resize(Aj,nnz)
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


class csr_matrix_batch_t(object):
    def __init__(self,datas,indices,indptr,shape=None):        
        self._datas,self.indices,self.indptr=datas,indices,indptr        
        self.mm=[csr_matrix((d,indices,indptr),copy=False,shape=shape) for d in datas ]
        
    def __getitem__(self,i):
        return self.mm[i]
    def __iter__(self):
        return iter(self.mm)
    
    def __len__(self):
        return len(self.mm)
    
    def __iadd__(self,o):
        datas= o.datas if type(o) == type(self) else o
        self._datas+=datas;
    
    @property
    def data(self):
        return self._datas;
    
    @property
    def shape(self):
        return self.mm[0].shape
    @property
    def nnz(self):
        return self.indptr[-1]
    @property
    def dtype(self):
        return self.mm[0].dtype
    @property
    def triplet(self):
        return self._datas,self.indices,self.indptr    

class coo_matrix_batch_t(object):
    
    def __init__(self,ds,row,col,shape=None):
        self.row,self.col=row,col
        '''
        self.datas,self.row,self.col=datas,row,col
        self.mm=[coo_matrix((d,(row,col)),shape=shape,copy=False) for d in datas ]
        self._csr_scheme=None;
        '''
        self._csr_scheme=None;
        self._init_data(ds,shape)
        
    @property
    def data(self):
        return self._datas;
    
    @data.setter
    def data(self,ds):
        self._init_data(ds,self._shape);
        
    def _init_data(self,ds,shape):
        row,col,nnz=self.row,self.col,self.nnz
        if (ds is not None)  and len(ds):
            ds=np.array(ds,copy=False).reshape(len(ds),-1)[:,:nnz];
        self._datas=ds;
            
        self.mm=mm=[coo_matrix((d,(row,col)),shape=shape,copy=False) for d in ds ]        
        self._shape = mm[0].shape if len(mm) else shape;
        
    
     
    def __getitem__(self,i):
        return self.mm[i]
    def __iter__(self):
        return iter(self.mm)
    
    def set_csr_scheme(self):
        csr_scheme=self._csr_scheme
        
        if csr_scheme is None:            
            tint=self.row.dtype
            idata=np.arange(self.nnz,dtype=tint);  
            row,col=self.row,self.col
            icoo=coo_matrix((idata,(row,col)),copy=False,shape=self.shape);
            #iFEM2pp_csr(trs=None,icoo=None,lbound=0,dtype=np.int32)
            self._csr_scheme=csr_scheme=iFEM2pp_csr(icoo=icoo,dtype=tint);
        return csr_scheme;
        
    
    def tocsr(self,datas=None,dtype=None):               
            
        if dtype is None:
            dtype=self.dtype;
            
        scheme=self.set_csr_scheme()
        if datas is None:
            datas=self._datas
            
        if dtype is None:
            dtype=datas.dtype;
            
    
        return FEM2data_csr(scheme,datas,dtype=dtype);
            
            
        
    @property
    def shape(self):
        return self._shape
    @property
    def nnz(self):
        return len(self.row)
    @property
    def dtype(self):
        if len(self.mm):
            return self.mm[0].dtype

        
        


def rcsz_iFEM2sparse_coo(trs,lbound=0,N=None,dtype=np.int32,fcoo=True):
    
    
    trs=np.array(trs,dtype=dtype,copy=False);
    
    if N is None:
        N,nd=trs.shape
    else:
        nd=trs.size//N;
        trs=trs.reshape((N,nd));
        
    #trs=trs.reshape(-1)
    lbound=int(lbound)    
    if lbound :
        trs=trs-lbound;
    
        
    Nf=N*nd
    ndndN=nd*nd*N;    
    e=np.ones(nd,dtype=dtype);
    row=np.kron(e,trs);
    col=np.kron(trs,e);
    
    if fcoo:
        [row,col]=[ v.reshape(-1) for v in [row,col]]
    
    #shape= None if Nd is None else  (Nd,Nd)
    
    return row,col,ndndN   
    




def iFEM2sparse_coo(trs,Nd=None,shape=None,lbound=0,dtype=np.int32):
    
    row,col,nnz = rcsz_iFEM2sparse_coo(trs,lbound,dtype=trs.dtype);    
    idata=np.arange(nnz,dtype=dtype);  
    return coo_matrix((idata,(row,col)),copy=False,shape=shape);
    
    


    
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
     
def iFEM2pp_csr(trs=None,icoo=None,lbound=0,dtype=np.int32):
    
    if icoo is None:
        icoo=iFEM2sparse_coo(trs,lbound=lbound,dtype=dtype)
        
    (pp,indices,indptr)=icoo2pp_csr(icoo,dtype=dtype)
    return (pp,indices,indptr);



    
def FEM2data_csr(csr_scheme,FEM_datas,dtype=None,shape=None):
    
    pp,indices,indptr=csr_scheme;    
    nnz=indptr[-1];
    
    FEM_datas=np.array(FEM_datas,copy=False).reshape(len(FEM_datas),-1);        
    K,N=FEM_datas.shape        
    
    if dtype is None:
        dtype=FEM_datas[0].dtype;
        
    datas=np.zeros((K,nnz),dtype=dtype);    
    
    
    permute_sum2(FEM_datas,datas,pp);
    
    return csr_matrix_batch_t(datas,indices,indptr,shape=shape) 
  

def FEM2data_coo(trs,FEM_datas=[],shape=None,lbound=0,dtype=None,copy=False):
      
    row,col,nnz = rcsz_iFEM2sparse_coo(trs,lbound,dtype=trs.dtype);
    
    if (FEM_datas is not None) and len(FEM_datas):
        if dtype is None:
            dtype=FEM_datas[0].dtype;
        FEM_datas=np.array(FEM_datas,copy=copy,dtype=dtype).reshape(len(FEM_datas),-1); 
        #FEM_datas=FEM_datas[:,:nnz]
    else: 
        FEM_datas=[];
    
    
    
    return coo_matrix_batch_t(FEM_datas,row,col,shape=shape); 


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

def coo_disjunct_sum(AA):    
    
    
    def _zero_datas():
        return [ np.zeros_like(m.data) for m in AA ]
    
    nnzs=[m.nnz for m in AA]
    nnz=np.sum(nnzs);  
    
    
    #zdatas=[ np.zeros_like(m.data) for m in AA ]
    
    rows=np.concatenate([ np.array(m.row) for m in AA ])    
    cols=np.concatenate([ np.array(m.col) for m in AA ])
    rAA=[]
    for m in range(len(AA)):
        #d=zdatas.copy();
        d=_zero_datas();
        d[m][:]=AA[m].data;
        d=np.concatenate(d)
        rAA+=[sp.coo_matrix((d,(rows,cols)),shape=AA[m].shape)]

    return rAA;

def sparse2coo_matrix_batch(*AA):
    
    AA=[coo_matrix(m) for m in AA]
    '''    
    AA=coo_disjunct_sum(AA)
    '''
    if not is_sparse_pattern_equ(*AA):
        AA=coo_disjunct_sum(AA)
    
    datas=[m.data for m in AA]        
    m=AA[0] 
    return coo_matrix_batch_t(datas,m.row,m.col,m.shape)
    
    
    
    
    

'''     
def dddo(datas):
    d=[d for d in datas  ]
    return np.array(d,dtype=object);
'''
        
if __name__=='__main__': 
    
    from utils import *
    from FEM.trs import *
    C=1
    D=2
    N=2**23
    N=2**16
    N=8
    
    trs=make_trs_1D(7,fcycle=1)
    datas=make_trs2datas_FEM(trs)
    coo=FEM2data_coo(trs,datas)
    
    
    r,c,nn=rcsz_iFEM2sparse_coo(trs)
    
    #N=4
    crand=lambda *ls: np.random.randn(*ls)+1j*np.random.randn(*ls)
    #datas=crand(D,N);
    datas=np.ones((D,N,2*C,2*C));
    
    
    tic()
    trs=make_trs_1D(N,C=C,fcycle=1);mtoc('make_trs_1D:')    
    tic()
    coo=FEM2data_coo(trs)
    ccsr=coo.tocsr(datas);
    
    coo=FEM2data_coo(trs,datas)
    mtoc('coo=FEM2data_coo(trs,datas):')
    nnz=coo[0].nnz
    tic()
    ccsr=coo.tocsr();
    mtoc('ccsr=coo.tocsr()')    
    coo=FEM2data_coo(trs,datas)
    tic()
    ccsr=coo.tocsr();
    mtoc('ccsr=coo.tocsr()')    
    tic()
    ccsr=coo.tocsr();
    mtoc('ccsr=coo.tocsr()[2]')    
    #print(coo[0].todense().real)
    
    tic()
    icoo=iFEM2sparse_coo(trs)   
    mtoc('iFEM2sparse_coo:')
    sh=iFEM2pp_csr(icoo=icoo);toc('iFEM2pp_csr:')
    
    
    
    
    #d=[d for d in datas  ]
    #datas=np.array(d,dtype=object);
    #sh=iFEM2pp_csr(trs);toc('iFEM2pp_csr:')
    
    tic();
    
    cc=FEM2data_csr(sh,datas);mtoc('FEM2data_csr:')
    
    print('icoo.nnz/N=',icoo.nnz/N)
    print('cc[0].nnz/N=',cc[0].nnz/N)
    tic();
    row,col,nnz = rcsz_iFEM2sparse_coo(trs,0,dtype=trs.dtype);
    mtoc('row,col,nnz = rcsz_iFEM2sparse_coo(trs,lbound,trs.dtype):')