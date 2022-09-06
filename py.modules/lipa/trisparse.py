# -*- coding: utf-8 -*-
"""
Created on Tue Jan 12 07:41:53 2016

@author: dich6_000
"""
import numpy as np
#import ctypes
import copy
from scipy import sparse as sp
from scipy.sparse import coo_matrix,csc_matrix,csr_matrix


class coo_scheme:
    def __init__(self,shape=None,row=np.array([], dtype=np.int32),col=np.array([], dtype=np.int32),data=None):
        self.format='coo'
        self.shape=shape
        self.row=row
        self.col=col
        self.data=data;
    @property
    def nnz(self):
        return np.size(self.row);
    @property
    def dtype(self):
        d=self.data;
        return None if d is None else d.dtype;
    def tocoo(self):
        return coo_matrix((self.data,(self.row,self.col)),shape=self.shape)
    def tocsc(self):
        return self.tocoo().tocsc();
    def tocsr(self):
        return self.tocoo().tocsr();
    def todense(self):
        return self.tocoo().todense();



def trimesh2sparse(trs,nvxs,data=None,nd=3,lbound=0):    
    
    
        
    trs=np.ascontiguousarray(trs,dtype=np.int32);
    
    if nd is None:
        nd=trs.shape[1];

    if not data is None:
        data=np.ascontiguousarray(data);
        data=data.reshape([data.size]);

    if lbound :        
        trs=trs-lbound;
    Nf=int(np.prod(trs.shape))
    N=int(Nf/nd);
    trs=trs.reshape(N,nd);
    ndndN=nd*nd*N;
    e=np.ones(nd,dtype=np.int32);
    row=np.kron(e,trs).reshape((ndndN,));
    col=np.kron(trs,e).reshape((ndndN,));
    return coo_scheme((nvxs,nvxs),row,col,data)

