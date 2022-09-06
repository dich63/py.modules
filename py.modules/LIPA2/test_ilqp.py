# -*- coding: utf-8 -*-
"""
Created on Fri Apr  1 05:24:30 2022

@author: wwww
"""

from LIPA2.ilaplace import *
from LIPA2.qp_solver import *

from utils import *

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

syms_def(locals(),' c0, c1 ,c0,c1,c2,c3,c4,b1, b2, b3, g, a1, a2, a3, a4, g1')

tic();

import dill
dill.settings['recurse'] = True

if 0:
    [p,d,q]=sym_reduce([(c0,2),(c1,2),(c3,2)],[( [0,b1,b2,b3] ,g  ),([a1,0,0,a4],g1)]);
    dill.dump([p,d,q],open(__file__+'.picle','wb'))
else:
    [p,d,q]=dill.load(open(__file__+'.picle','rb'))
    
#[p,d,q]=sym_reduce([c0 ,c1,c0,c1,c2,c3],[( [0,b1,b2,b3] ,g  )]);

toc('reduce:')
#jo=arg2jso(c1=1.1,c2=2.3,c3=1.11,c4=21,c0=10j,g=3.11,g1=2j,b1=2.1,b2=1.2,b3=1.1j,a4=2)
jo={
    c1:-0.01+3j,
    c2:-3.2,
    c3:-0.13+2j,
    c4:-21,
    c0:-0.01-3j,#-0.01+2.1j,
    g:-1.051+3j,        
    b1:2.1e-3,
    b2:1.2,
    b3:10+1.1j,
    g1:-1+2j,
    a1:0.0,
    a4:2.0e-2}



tic();

[GF,exg]=ilaplace_functor(p,jo)
[F,exf]=ilaplace_functor(q,jo)
toc('functors:')



tt=np.linspace(0,100,1200)
rn=lambda x: x.real/np.max(np.abs(x))
#plt.plot(tt,rn(F(tt)),tt,rn(GF(tt)))
plt.plot(tt,rn(F(tt)),tt,rn(GF(tt)))


[DC,QP,G,FF]=make_lipa_data(d,q,jo);

LM=[4,8]
dt=0.2
print(QP.shape)
print(G.shape)
print(FF.shape)
lqp=lipa_qp_number(DC,LM=LM,FF=FF,qp=QP,g=G).reset(dt);
x=GF(dt)
y=lqp()

print(x)    
print(y[0])    




