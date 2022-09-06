# -*- coding: utf-8 -*-
"""
Created on Sun May 22 18:27:51 2022

@author: wwww
"""

#%reset -sf


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
#fn=r'O:\__ss\matrix\sfFEM128k.json'

#fn=r'O:\__ss\matrix\sfFEM24k.json'
fn=r'O:\__ss\matrix\KGM.json'
d=decode(fn,1)
k=d.K
k=k.tocsr();
k=k.tocsr();
q=k.has_sorted_indices

print('decode ',fn)

#K,G,M=[m.tocsr() for m in [d.K,d.G,d.M]]

zz=[7+17j,2+1j,2+4j,2-4j,1+1j,2+1j,2+4j,2-4j]
zz=[7+17j,2+1j,2+4j,2-4j]
pp=1
mp=16
print('mp',mp)
zz=[7+17j for k in range(mp)]
nt=0
#zz=[1+1j,2+1j]
#    zz=[7+17j,]

N=d.M.shape[0]
di=np.ones(N-1)

#G=sp.tril(G,format='csc')
#ab=[K+z*G+z**2*M for z in zz ];
K,G,M=[m.tocoo() for m in [d.K,d.G,d.M]]

#k=K.todia()

daz=np.empty_like(K.data,dtype=np.complex128)
ab=[sp.coo_matrix((copy.copy(daz),(K.row,K.col)),shape=K.shape)  for z in zz ]

ds=[K.data+z*G.data+z**2*M.data for z in zz]

for k in range(len(ds)):
    ab[k].data[:]=ds[k];

#ab=[K+z*G+z**2*M for z in zz ];

#ab=[m.tocsr() for m in ab ];
ab[0]=ab[0]+sp.tril(ab[0]);


w=ab[0].tocoo()
w=w.tocsr()
ab=[m.tocsr() for m in ab ];
a0=ab[0]
a1=ab[1]
'''
a1=a1.tocoo()
a0=a0.tocoo()
'''
kk=a0+a1
'''
Ma=decode(r'O:\__ss\matrix\M.json').tocsr()
ab=[Ma for m in range(mp)];
'''

indptr,indices=ab[0].indptr,ab[0].indices
datas=[m.data for m in ab ];

cc=sp.csr_matrix((datas[0],indices,indptr));

N=ab[0].shape[0]

common=common_update(ordering=0,scale=-1,btf=0,tol=1e-15);


print('start klu ')

hsymbolic=smart_ptr_t(h_context_t())
#hsymbolic=link_ptr_t(h_context_t())


tclass=np.int64(0x0109)
n=ab[0].shape[0]

tic();hsymbolic=iklu_symbolic_batch(n,indptr,indices,tclass,common);toc('_klub_symbolic_batch:')


print('iklu_symbolic=',hsymbolic.state);

print('pp=',pp);
tic();factors=iklu_numeric_batch(n,indptr,indices,datas,tclass,hsymbolic,parallel=pp);toc('iklu_numeric_batch:')   

print('iklu_numeric=',factors.states);


xx=np.random.randn(mp,n)+1j*np.random.randn(mp,n)

bb=1*xx;
pbb=to_handles(bb)
print('pid=',os.getpid())
st=iklu_tsolve_ts_batch(factors,pbb,parallel=pp)
tic();st=iklu_tsolve_ts_batch(factors,pbb,parallel=pp) ;t0=toc('iklu_solve_ts_batch1:')
bb[:]=xx
common=klu_common_t(user_data=0xbabaeb)
#print("common.user_data=",common.user_data)
print('pid=',os.getpid())
tic();st=iklu_tsolve_ts_batch(factors,pbb,parallel=pp,common=common) ;t1=toc('iklu_solve_ts_batch2:')
#print("common.user_data=",common.user_data)
print('iklu_solve_ts_batch=',st)
print('perfQP=',t0/t1)


bb[:]=xx
bbuf=1*xx;
print('pid=',os.getpid())

(err,Q,P)=get_symbolic_permutes(hsymbolic.ptr)
(Ps,Rs)=iklu_numeric_permuts(factors)

tic();st=iklu_tsolve_ts_batch_pp(factors,bb,bbuf,(Q,Ps,Rs)) ;t1=toc('iklu_solve_ts_batchpp:')
#print("common.user_data=",common.user_data)
print('iklu_solve_ts_batch=',st)
print('perfQP1=',t0/t1)


err=0.0
for m in range(mp):
    err+=norm(ab[m].dot(bb[m]) -xx[m])/norm(xx[m]);
    
print('err=',err/mp)

fs=factors.subset(range(1))
tic()
for k in range(mp):
    st=iklu_tsolve_ts_batch(fs,pbb,parallel=0) ;
t2=toc('iklu_solve_ts_batch:')
print('perf0=',t2/t0)
print('perf=',t2/t1)
'''
(err,Q,P)=get_symbolic_permutes(hsymbolic.ptr)
(err,Pn,Rs)=get_number_PRS(factors[0])


i=np.random.permutation(np.arange(n));
q=bb[0]

tic();p=q[i];q[:]=p[i];toc('perm2')
'''

