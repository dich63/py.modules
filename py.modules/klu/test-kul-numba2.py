# -*- coding: utf-8 -*-
"""
Created on Sun May 22 18:27:51 2022

@author: wwww
"""

%reset -sf


from klu.klu_numba import *

    
    
    
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
fn=r'O:\__ss\matrix\sfFEM128k.json'
#fn=r'O:\__ss\matrix\KGM.json'
d=decode(fn,1)
k=d.K
k=k.tocsr();
k=k.tocsr();
q=k.has_sorted_indices

#K,G,M=[m.tocsr() for m in [d.K,d.G,d.M]]

zz=[7+17j,2+1j,2+4j,2-4j,1+1j,2+1j,2+4j,2-4j]
zz=[7+17j,2+1j,2+4j,2-4j]
pp=1
mp=8
zz=[7+17j for k in range(mp)]
nt=0
#zz=[1+1j,2+1j]
#    zz=[7+17j,]

N=d.M.shape[0]
di=np.ones(N-1)

#G=sp.tril(G,format='csc')
#ab=[K+z*G+z**2*M for z in zz ];
K,G,M=[m.tocoo() for m in [d.K,d.G,d.M]]

daz=np.empty_like(K.data,dtype=np.complex128)
ab=[sp.coo_matrix((copy.copy(daz),(K.row,K.col)),shape=K.shape)  for z in zz ]

ds=[K.data+z*G.data+z**2*M.data for z in zz]

for k in range(len(ds)):
    ab[k].data[:]=ds[k];

#ab=[K+z*G+z**2*M for z in zz ];

#ab=[m.tocsr() for m in ab ];
ab=[m.tocsc() for m in ab ];

indptr,indices=ab[0].indptr,ab[0].indices
datas=[m.data for m in ab ];

common=common_update(ordering=0,scale=-1,btf=0,xtol=0.1);



hsymbolic=smart_ptr_t(h_context_t())
#hsymbolic=link_ptr_t(h_context_t())


tclass=np.int64(0x0109)
n=ab[0].shape[0]

'''
tic();hsymbolic=iklu_symbolic_batch(n,indptr,indices,tclass,common);toc('_klub_symbolic_batch:')


print('iklu_symbolic=',hsymbolic.state);

print('pp=',pp);
tic();factors=iklu_numeric_batch(n,indptr,indices,datas,tclass,hsymbolic,parallel=pp);toc('iklu_numeric_batch:')   

print('iklu_numeric=',factors.states);
'''

factors=to_handles(release=iklu_release,addref=iklu_addref);

for k in range(mp):
    tic();hsymbolic=iklu_symbolic_batch(n,indptr,indices,tclass,common);toc('_klub_symbolic_batch:')
    tic();fs=iklu_numeric_batch(n,indptr,indices,datas[0:1],tclass,hsymbolic,parallel=pp);toc('iklu_numeric_batch:')   
    factors=factors+fs


xx=np.random.randn(mp,n)+1j*np.random.randn(mp,n)

bb=1*xx;
pbb=to_handles(bb)

tic();st=iklu_solve_ts_batch(factors,pbb,parallel=pp) ;toc('iklu_solve_ts_batch:')
print('iklu_solve_ts_batch=',st)
err=0.0
for m in range(mp):
    err+=norm(ab[m].dot(bb[m]) -xx[m])/norm(xx[m]);
    
print('err=',err/mp)


