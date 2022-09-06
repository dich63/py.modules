# -*- coding: utf-8 -*-
"""
Created on Fri Apr  1 05:24:30 2022

@author: wwww
"""

from LIPA2.ilaplace import *
from LIPA2.qp_solver import *

from utils import *
from utils.interpolate_1D import *

import matplotlib.pyplot as plt
fig=plt.figure(figsize=(12,8))
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
    
    g:-17.5+250j,
    #g:0,                
    b1:0,
    b2:-1j,
    b3:0,
    b4:1,
    }



tic();
nd=0
[GF,exg]=ilaplace_functor(p,jo,nd)
[F,exf]=ilaplace_functor(q,jo)
toc('functors:')



tt=np.linspace(0,10,1200)
rn=lambda x: x.real/np.max(np.abs(x.real))
#plt.plot(tt,rn(F(tt)),tt,rn(GF(tt)))
plt.plot(tt,rn(F(tt)),tt,rn(GF(tt)))
plt.show()


[DC,QP,G,FF]=make_lipa_data(d,q,jo);

LM=[4,8]
dt=0.02
print(QP.shape)

print(G.shape)
print(FF.shape)
lqp=lipa_qp_number(DC,LM=LM,FF=FF,qp=QP,g=G,nd=nd+1).reset(dt);

nt=200
tt=dt*np.arange(1,nt+1)
tm=np.max(tt)
ttx=np.linspace(0,tm,10000)
x=GF(tt)
ynn=lqp.dump(nt)
ynn=ynn.reshape(ynn.shape[0],-1).T
y=ynn[nd]

jG=jet_spline(tt,ynn,3)

err=norm(x-y)/norm(x);
plt.plot(ttx,GF(ttx),ttx,jG(ttx))
plt.plot(tt,y, marker='o',linewidth=0)
plt.show()
#print(x)    
#print(y[0])    
printf('error=%3.4g%%\n',100*err)

tic()
[p,d,q]=sym_reduce([c0,c1,c2,c3],[( [b1,b2,b3,b4] ,g  )]);
toc('[p,d,q]=sym_reduce([c0,c1,c2,c3],[( [b1,b2,b3,b4] ,g  )])->')

