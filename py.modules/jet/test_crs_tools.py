# -*- coding: utf-8 -*-
"""
Created on Sat Jul  2 15:12:58 2022

@author: DICH
"""
import numpy as np



from jet.tools import *
import LIPA2.tools as l2ts
from FEM.FEM2sparse import *

from jsonrpc.json_io import *
import os
import scipy
import scipy.sparse as sp


import jet.csr_tools

#jet.csr_tools.csr_matrix =sp.csr_matrix

from jet.csr_tools import *
from jet.nb_tools import *
#import jet.nb_tools as nbt



nan=np.nan
cnan=nan+1j*nan
norm=np.linalg.norm
normm= lambda x:  norm(x.reshape(-1),ord=np.inf)
crand=lambda *ls : np.random.randn(*ls)+1j*np.random.randn(*ls)   

M=8;
D=4;

fn=r'O:\__ss\matrix\sfFEM0k.json'
fn=r'O:\__ss\matrix\sfFEM128k.json'
#fn=r'O:\__ss\matrix\KGM.json'
d=decode(fn,1)
mK,mG,mM=d.K,d.G,d.M
mK,mG,mM=[m.tocoo() for m in [mK,mG,mM] ]
datas=[mK.data,mG.data,mM.data]
#coo_matrix_batch_t(FEM_datas,row,col,shape=shape); 
CD=coo_matrix_batch_t(datas,mK.row,mK.col)
CDcsr=CD.tocsr(dtype=complex)
N=CDcsr.shape[0]
xx0=crand(1+D,N);
xx0[:]=1
xx0[2]=1
x0=xx0[0];
x0D=xx0[1:];

z=10j
zp=np.zeros(M,dtype=complex)
zp[:]=z;

xxz_out1=cnan*np.ones([M,D+1,N])
#xxz_out1[0]=x0
xxz_out2=cnan*np.ones([M,D,N])
Czxx_out1=cnan*np.ones([M,N])

Czxx_out2=1j*np.zeros([M,N])

#xxz_out,Czxx_out=
for m in range(M):
    l2ts.AzC(xx0,CDcsr,zp[m],xxz_out1[m],Czxx_out1[m])

minus_invAz(N,D,M,zp,x0,x0D,xxz_out2)
datas,indices,indptr=CDcsr.triplet
datas=datas[1:]
Dm=len(datas)
idatas=np.arange(Dm)
xxzD=xxz_out2;
yz=Czxx_out2

#csr_gaxpy_jetz_np(N,Dm,M,indptr,indices,datas,idatas,x0,xxzD,yz)    
csr_gaxpy_jetz(N,Dm,M,indptr,indices,datas,idatas,x0,xxzD,yz)    
#csr_gaxpy(N,indptr,indices,datas,idatas,x0,yz[0])

#yz0=CDcsr[1].dot(x0)

print('err=',normm(Czxx_out1-yz))

#csr_gaxpy_jetz(

