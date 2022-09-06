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

%matplotlib auto
fig=plt.figure(figsize=(18,18))
#fig.suptitle('sensors data + pure response ', fontsize=16)
plt.grid(True,which='major')
plt.minorticks_on()
plt.grid(True,which='minor',alpha=0.2)
'''



'''
import numpy as np
norm=np.linalg.norm

syms_def(locals(),' c0, c1 ,c0,c1,c2,c3,c4,b5,b0,b1, b2, b3, b4, g, a1, a2, a3, a4, g1')

tic();

import dill
dill.settings['recurse'] = True

if 0:
    #[p,d,q]=sym_reduce([g],[( [b1,b2,b3,b4] ,g  )]);
    [p,d,q]=sym_reduce([c0,c1,c2,c3],[( [b0,b1,b2,b3,b4] ,g  )]);
    dill.dump([p,d,q],open(__file__+'.picle','wb'))
else:
    [p,d,q]=dill.load(open(__file__+'.picle','rb'))
    
#[p,d,q]=sym_reduce([c0 ,c1,c0,c1,c2,c3],[( [0,b1,b2,b3] ,g  )]);

toc('reduce:')
#jo=arg2jso(c1=1.1,c2=2.3,c3=1.11,c4=21,c0=10j,g=3.11,g1=2j,b1=2.1,b2=1.2,b3=1.1j,a4=2)

jo={
    
    c0:-.5+250j,
    c1:-.5-265.5j,
    c2:-1.5+258j,
    c3:-1.5-268j,
    c4:-21,
    
    g:-3777.5+8110j,
    #g:0,                
    b0:0,
    b1:0,
    b2:0,
    b3:1j,
    b4:1,
    }



tic();
nd=0

t0=0.001
sc=1

LM=[4,12]
LM=[2,5]
dt=0.01
LM=[8,18];
dt=0.0225*1.8*1;

#t0=0

#dt=0.0225*1.0/2;
#dt=1

'''
LM=[4,8];
LM=[4,12];
dt=0.00225*1;
'''
t0=dt/8

[GF,exg]=ilaplace_functor(p,jo,nd)
[F,exf]=ilaplace_functor(q,jo)
toc('functors:')


tR=1
tt=np.linspace(0,tR,120000)
rn=lambda x: x.real/np.max(np.abs(x.real))
rn=lambda x: x.real
#plt.plot(tt,rn(F(tt)),tt,rn(GF(tt)))
#plt.plot(tt,rn(F(tt)-F(tt-t0)),tt,rn(GF(tt)-GF(tt-t0)))
#plt.plot(tt,(F(tt)-F(tt-t0)),tt,F(tt))



GFt0=lambda t: GF(t)-sc*GF(t-t0)
Ft0=lambda t: F(t)-sc*F(t-t0)

ma=max(np.abs(F(tt).real))
plt.plot([t0,t0],[-ma,ma],color='#AAAAAA')
plt.plot([dt,dt],[-ma,ma],color='#AAAAAA')
plt.plot(tt,rn(F(tt)),color='#AAFFAA')
plt.plot(tt,rn(F(tt-t0)),color='#FFAAFF')
plt.plot(tt,rn(Ft0(tt)))


plt.show()


[DC,QP,G,FF]=make_lipa_data(d,q,jo);





print(QP.shape)
print(G.shape)
print(FF.shape)
D=DC.size-1
D=12

#poles,res=pade_exp_poles_res(LM,t=dt);

#QP=np.array([[0j,0]]);
#QP=[[0]]
cQP=QP.copy()
QP2=QP.copy()



#lqp=lipa_qp_number(DC,LM=LM,FF=FF,qp=cQP,g=G,nd=D+1).reset(dt);
FF2=np.array([1+0j,1]).T
lqp=jet_csr_klu_number(DC,D=D,FF=FF,qp=cQP,g=G).reset(dt=dt,LM=LM);
jet=jet_csr_klu_number(DC,D=D,FF=[FF[0],FF[0]],qp=[1*QP[0],-sc*QP[0]],g=[G[0],G[0]],t0=[0,t0]).reset(dt=dt,LM=LM);
#jet=jet_csr_klu_number(DC,D=D,FF=FF,qp=[QP[0]],g=[G[0]],t0=[0]).reset(dt=dt,LM=LM);

#lqp.x=1+0j
#jet.x=1+0j
#imzq=np.array([p.qp_laplace.imz for p in reversed(lqp.psolvers)])
#imzq=np.array([p.qp_laplace.imz for p in lqp.psolvers])
#imzj=jet.source.qpl.imz



'''
yl=lqp()
yj=jet()
print('norm(yj-yl)/norm(yj)=',norm(yj-yl)/norm(yj))
raise SystemExit();
'''
nt=int(tR/dt)
#nt=400
tt=dt*np.arange(0,nt+1)
tm=np.max(tt)
ttx=np.linspace(dt,tm,10000)
ttx=np.linspace(0,tm,10000)
x=GF(tt)
x0=GFt0(tt)
#x=GF(tt)-GF(tt-t0)
tic('start ynn=lqp.dump(nt)')
ynn=lqp.dump(nt)
toc('ynn=lqp.dump(nt)')
tic('start yjnn=jet.dump(nt)')
yjnn=jet.dump(nt)
toc(' yjnn=jet.dump(nt)')
ynn=ynn.reshape(ynn.shape[0],-1).T
#yjnn=yjnn.reshape(yjnn.shape[0],-1).T
yjnn0=yjnn.reshape(yjnn.shape[0],-1).T
d1,d2=yjnn0.shape;
yjnn=np.zeros((d1,d2+1),dtype=complex)
yjnn[:,1:]=yjnn0



y=ynn[nd]
yj=yjnn[nd]
jG=jet_spline(tt,yjnn,D)
fu= lambda x: np.abs(x)
fu= lambda x: x



fig2=plt.figure(figsize=(18,18))
#fig.suptitle('sensors data + pure response ', fontsize=16)
plt.grid(True,which='major')
plt.minorticks_on()
plt.grid(True,which='minor',alpha=0.2)

ma=max(np.abs(fu(GFt0(ttx))))
#plt.plot(ttx,fu(GF(ttx)),label='$G^{ex}$',color='#AAAAAA')
#plt.plot(ttx,fu(GF(ttx-t0)),label='$G^{t0}$')
plt.plot(ttx,fu(GFt0(ttx)),label='$G^{ex}_{full}$')
plt.plot(ttx,fu(jG(ttx)),label='$Spline_{t}$')
#plt.plot(tt,fu(y), marker='o',linewidth=0)
plt.plot(tt,fu(yj), marker='o',linewidth=0,label='$Jet_{t}$')
plt.plot([t0,t0],[-ma,ma],label='$t_{0}$',color='#AAAAAA')
plt.plot([dt,dt],[-ma,ma],label='$\delta t$',color='#AAAAAA')
plt.title('aa')
legend = plt.legend(loc='upper right', shadow=True, fontsize='x-large')


#print(x)    
#print(y[0])    
err=norm(x[1:]-y)/norm(x[1:]);
printf('error=%3.2e%%\n',100*err)

errj=norm(x0-yj)/norm(x0);
printf('errorj=%3.2e%%\n',100*errj)

plt.title(sprintf('error=%3.1e%% $error_{spline}$=%3.1e%%\n',100*err,100*errj))
plt.show()
'''
errjq=norm(ynn-yjnn)/norm(ynn)
printf('err-jet-qp=%3.2e%%\n',errjq*100)

tic()
[p,d,q]=sym_reduce([c0,c1,c2,c3],[( [b1,b2,b3,b4] ,g  )]);
toc('[p,d,q]=sym_reduce([c0,c1,c2,c3],[( [b1,b2,b3,b4] ,g  )])->')
'''
