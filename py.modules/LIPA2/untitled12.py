# -*- coding: utf-8 -*-
"""
Created on Wed Apr 27 14:04:41 2022

@author: wwww
"""


import numpy as np
import scipy.sparse as sp
import copy,os
import klu_batch.klub as klub

norm=np.linalg.norm
normm= lambda x:  norm(x,ord=np.inf)

N=5;A = sp.rand(N, N, density=0.3, format='csr')

print(A.todense())
cA=(1+0j)*A

cA.data[:]=range(1,cA.nnz+1)
cAc=copy.deepcopy(cA)
cAc.data[:]=np.nan

print(cA.todense().real)
print(cAc.todense().real)

Acoo=cA.tocoo()

sA=klub.make_spm(Acoo)
sAc=klub.make_spm(cAc)
print('pid=',os.getpid())

st=klub._csr2csc_batch(sA.ptr,sAc.ptr,0)

print(cA.todense().real)
print(cAc.todense().real.T)

print('err=',norm(cA.todense()-cAc.todense().T))

print(st)

raise SystemExit(0)
i,p=A.indices,A.indptr

ix=np.array(range(1,A.nnz+1))

ia=sp.csr_matrix((ix,i,p))

iac=ia.tocsc()
iatc=sp.csr_matrix((iac.data,iac.indices,iac.indptr))
iat=sp.csc_matrix((iac.data,iac.indices,iac.indptr))

ixp=list(iac.data-1)


aa=sp.csc_matrix((ix[ixp],iac.indices,iac.indptr))
Ar=sp.csc_matrix((A.data[ixp],iac.indices,iac.indptr))

print(ia.todense())
print(aa.todense())

'''
At=A.T.tocsr()
it,pt=At.indices,At.indptr

Att=sp.csr_matrix((ix,it,pt))
iix=ix-1
AT=sp.csr_matrix((A.data[iix],it,pt))
'''