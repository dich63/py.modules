
#
import sys;sys.path.append('v:/ipc/py.modules')

import env;

import  parallel.sparse as psp 
import numpy as np
from asyn.SharedWorker import Tic

import jsonrpc.sparse_marshal
import jsonrpc.jsonclass as jsc

from lipa.kernel import pade_exp_poles_res

import copy;

tic=Tic();

norm=lambda x : np.linalg.norm(x,np.inf) 
#norm=lambda x : np.linalg.norm(x) 
Ntic=40;
Ms=1e-3
t=0.01*Ms;
prz=pade_exp_poles_res.get((2,4))
zt=prz[0][0]/t;


A=[
    [0.,  2, 3,4 ],
    [5, 1., 7,8 ],
    [9, 10, 1,12],
    [13,14, 15,50]
    ]


nA=np.array(A)

b=np.array([1,2,3,4+1j]);

x=np.dot(nA,b);



#print(' x=',x)



jfile='c:/temp/fro_test_ec3.json'
#jfile='c:/temp/fro_test.json'
print('loading "',jfile,'"...');
tic.start()
s=jsc.decode_from_file(jfile,1)

nA=zt*s.C+s.A;
n=nA.shape[0]
b=(1.0+np.random.rand(n))+1j*(0.0-np.random.rand(n));
nA=nA.tocsc();
x=nA*b;
t=tic.sec();
print('loading end: ',round(t,5),' sec');
print(' |x|=',norm(x))

opts={"common":{"numerical_rank":-4,"scale":-1}}
print('opts:',opts)
sol=psp.sparse_solver_invoker(nA,opts=opts)
tic.start()
print('factorize start...');
sol.factorize()
t=tic.sec();
print('factorize end: ',round(t,5),' sec');

sol.y=x[:];
sol.solve()
y=copy.copy(sol.y);
tic.start()
print('solve  start...');

for k in range(1,Ntic):
    sol.solve()
t=tic.sec();
print('solve end: ',round(t,5),' sec ',Ntic,' times');



print(' err=|y-b|/|b|=',norm(y-b)/norm(b));

from utils.derr2m import derr2m
derr=derr2m(y,b);

#print(' derr2m(y,b)=',100*np.mean(np.abs(derr)),'%');
print(' derr2m(y,b)=',100*np.max(np.abs(derr)),'%');
cn=sol.cond;
print('estimate condnum={:e}'.format(cn));
import ctypes as cp
common=psp.klu_common_t()
err=sol.context('i',cp.byref(common))
mb=1.*2**20
print('mempeak=',common.mempeak/mb,' mb')
print('memusage=',common.memusage/mb,' mb')
t=11;

#print(' derr2m(y,b)=',100*(np.abs(derr)));