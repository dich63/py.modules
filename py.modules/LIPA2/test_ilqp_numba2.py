# -*- coding: utf-8 -*-
"""
Created on Fri Apr  1 05:24:30 2022

@author: wwww
"""

from LIPA2.ilaplace import *
from LIPA2.qp_solver import *

from jet.jet_csr_iklu import *

from utils import *
from utils.interpolate_1D import *

import matplotlib.pyplot as plt
fig=plt.figure(figsize=(18,18))
#fig.suptitle('sensors data + pure response ', fontsize=16)
plt.grid(True,which='major')
plt.minorticks_on()
plt.grid(True,which='minor',alpha=0.2)
'''



'''
import numpy as np
norm=np.linalg.norm

syms_def(locals(),' c0, c1 ,c0,c1,c2,c3,c4,b1, b2, b3, b4, g, a1, a2, a3, a4, g1')

tic();

import dill
dill.settings['recurse'] = True

if 0:
    #[p,d,q]=sym_reduce([g],[( [b1,b2,b3,b4] ,g  )]);
    [p,d,q]=sym_reduce([c0,c1,c2,c3],[( [b1,b2,b3,b4] ,g  )]);
    dill.dump([p,d,q],open(__file__+'.picle','wb'))
else:
    [p,d,q]=dill.load(open(__file__+'.picle','rb'))
    
#[p,d,q]=sym_reduce([c0 ,c1,c0,c1,c2,c3],[( [0,b1,b2,b3] ,g  )]);

toc('reduce:')
#jo=arg2jso(c1=1.1,c2=2.3,c3=1.11,c4=21,c0=10j,g=3.11,g1=2j,b1=2.1,b2=1.2,b3=1.1j,a4=2)
jo={
    
    c0:-5.5+250j,
    c1:-5.5-250.5j,
    c2:-1.5+53j,
    c3:-1.5-53j,
    c4:-21,
    
    g:-57.5+260j,
    #g:0,                
    b1:0,
    b2:-1j,
    b3:0,
    b4:1,
    }


jo={
    
    c0:-.5+250j,
    c1:-.5-265.5j,
    c2:-1.5+258j,
    c3:-1.5-268j,
    c4:-21,
    
    g:-57.5+260j,
    #g:0,                
    b1:10,
    b2:-1j,
    b3:3,
    b4:11,
    }




tic();
nd=0

[GF,exg]=ilaplace_functor(p,jo,nd)
[F,exf]=ilaplace_functor(q,jo)
toc('functors:')


tR=10
tt=np.linspace(0,tR,1200)
rn=lambda x: x.real/np.max(np.abs(x.real))
#plt.plot(tt,rn(F(tt)),tt,rn(GF(tt)))
plt.plot(tt,rn(F(tt)),tt,rn(GF(tt)))
plt.show()


[DC,QP,G,FF]=make_lipa_data(d,q,jo);




LM=[4,12]
LM=[2,5]
dt=0.01
LM=[8,18];
dt=0.0225*1.8;
#dt=1
'''
LM=[4,8];
dt=0.00225*1;

'''
print(QP.shape)
print(G.shape)
print(FF.shape)
D=DC.size-1


#poles,res=pade_exp_poles_res(LM,t=dt);

#QP=np.array([[0j,0]]);
#QP=[[0]]
cQP=QP.copy()

lqp=lipa_qp_number(DC,LM=LM,FF=FF,qp=cQP,g=G,nd=D+1).reset(dt);
jet=jet_csr_klu_number(DC,D=D,FF=FF,qp=QP,g=G).reset(dt=dt,LM=LM);

#lqp.x=1+0j
#jet.x=1+0j
#imzq=np.array([p.qp_laplace.imz for p in reversed(lqp.psolvers)])
imzq=np.array([p.qp_laplace.imz for p in lqp.psolvers])
imzj=jet.source.qpl.imz



'''
yl=lqp()
yj=jet()
print('norm(yj-yl)/norm(yj)=',norm(yj-yl)/norm(yj))
raise SystemExit();
'''
nt=int(tR/dt)
#nt=400
tt=dt*np.arange(1,nt+1)
tm=np.max(tt)
ttx=np.linspace(dt,tm,10000)
x=GF(tt)
ynn=lqp.dump(nt)
yjnn=jet.dump(nt)

ynn=ynn.reshape(ynn.shape[0],-1).T
yjnn=yjnn.reshape(yjnn.shape[0],-1).T

y=ynn[nd]
yj=yjnn[nd]
jG=jet_spline(tt,yjnn,D)
fu= lambda x: np.abs(x)
fu= lambda x: x
err=norm(x-y)/norm(x);
plt.plot(ttx,fu(GF(ttx)),ttx,fu(jG(ttx)))
plt.plot(tt,fu(y), marker='o',linewidth=0)
plt.plot(tt,fu(yj), marker='x',linewidth=0)
plt.show()
#print(x)    
#print(y[0])    
printf('error=%3.4e%%\n',100*err)
errjq=norm(ynn-yjnn)/norm(ynn)
printf('err-jet-qp=%3.4e%%\n',errjq*100)
'''
tic()
[p,d,q]=sym_reduce([c0,c1,c2,c3],[( [b1,b2,b3,b4] ,g  )]);
toc('[p,d,q]=sym_reduce([c0,c1,c2,c3],[( [b1,b2,b3,b4] ,g  )])->')
'''
