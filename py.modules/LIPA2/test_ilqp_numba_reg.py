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





import numpy as np
norm=np.linalg.norm

syms_def(locals(),' c0, c1 ,c0,c1,c2,c3,c4,b1,b0, b2, b3, b4, g, a1, a2, a3, a4, g1')

tic();


#[p,d,q]=sym_reduce([c0 ,c1,c0,c1,c2,c3],[( [0,b1,b2,b3] ,g  )]);
[p,d,q]=sym_reduce([c0],[( [b0,b1,b2,b3,b4] ,g  )]);
toc('reduce:')
#jo=arg2jso(c1=1.1,c2=2.3,c3=1.11,c4=21,c0=10j,g=3.11,g1=2j,b1=2.1,b2=1.2,b3=1.1j,a4=2)

jo={
    
    c0:0*(-5.5+250j),
    
    g:-7.5+0*1150j,
    #g:0,                
    b0:.0,
    b1:0,
    b2:0,
    b3:1,
    b4:0,
    }




tic();
nd=0


sc=1

LM=[4,12]
LM=[2,5]
dt=0.01
LM=[8,18];
dt=0.0225*1.8*1;
#dt=0.0225*1.0/2;
#dt=1


LM=[2,4];
dt=0.00225*3;
'''
'''

t0=dt/2
t0=0

[GF,exg]=ilaplace_functor(p,jo,nd)
[F,exf]=ilaplace_functor(q,jo)
toc('functors:')


tR=1
tt=np.linspace(0,tR,12000)
rn=lambda x: x.real/np.max(np.abs(x.real))
rn=lambda x: x.real
#plt.plot(tt,rn(F(tt)),tt,rn(GF(tt)))
#plt.plot(tt,rn(F(tt)-F(tt-t0)),tt,rn(GF(tt)-GF(tt-t0)))
#plt.plot(tt,(F(tt)-F(tt-t0)),tt,F(tt))





[DC,QP,G,FF]=make_lipa_data(d,q,jo);





print(QP.shape)
print(G.shape)
print(FF.shape)
D=DC.size-1
D=12

#jet=jet_csr_klu_number(DC,D=D,FF=FF,qp=QP,g=G,t0=[t0]).reset(dt=dt,LM=LM);
jet=jet_csr_klu_number(DC,D=D,FF=FF,qp=QP,g=G).reset(dt=dt,LM=LM);

y1=jet()[0][0]
x1=GF(dt)
print('x1=',x1)
print('y1=',y1)
print('err=',norm(y1-x1)/norm(x1))

fy=jet.source.qpm.qp[0][0]
fx=F(dt)
print('fx=',fx)
print('fy=',fy)
print('errf=',norm(fy-fx)/norm(fx))


raise SystemExit()

GFt0=lambda t: GF(t-t0)
Ft0=lambda t: F(t-t0)

ma=max(np.abs(F(tt)))
plt.plot([t0,t0],[-ma,ma],color='#AAAAAA')
plt.plot([dt,dt],[-ma,ma],color='#AAAAFF')
plt.plot(tt,rn(F(tt)),color='#AAAAAA')
plt.plot(tt,rn(Ft0(tt)),color='#FF0000')


plt.show()



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

import matplotlib.pyplot as plt
%matplotlib auto


fig=plt.figure(figsize=(18,18))
#fig.suptitle('sensors data + pure response ', fontsize=16)
plt.grid(True,which='major')
plt.minorticks_on()
plt.grid(True,which='minor',alpha=0.2)
'''




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
tic('start yjnn=jet.dump(nt)')
yjnn=jet.dump(nt)
toc(' yjnn=jet.dump(nt)')

#yjnn=yjnn.reshape(yjnn.shape[0],-1).T
yjnn0=yjnn.reshape(yjnn.shape[0],-1).T
d1,d2=yjnn0.shape;
yjnn=np.zeros((d1,d2+1),dtype=complex)
yjnn[:,1:]=yjnn0

yj=yjnn[nd]
jG=jet_spline(tt,yjnn,D)
fu= lambda x: np.abs(x)
fu= lambda x: x
#fu= lambda x: x/np.mean(np.abs(x))


fig2=plt.figure(figsize=(18,18))
#fig.suptitle('sensors data + pure response ', fontsize=16)
plt.grid(True,which='major')
plt.minorticks_on()
plt.grid(True,which='minor',alpha=0.2)

ma=max(np.abs(fu(GF(ttx)).real))
plt.plot(ttx,fu(GF(ttx)),label='$G^{ex}$',color='#AAAAAA')
#plt.plot(ttx,fu(GF(ttx-t0)),label='$G^{t0}$')
plt.plot(ttx,fu(GFt0(ttx)),label='$G^{ex}_{full}$')
#plt.plot(ttx,fu(jG(ttx)),label='$Spline_{t}$')
#plt.plot(tt,fu(y), marker='o',linewidth=0)
plt.plot(tt,fu(yj), marker='o',linewidth=0,label='$Jet_{t}$')
'''
plt.plot([t0,t0],[-ma,ma],label='$t_{0}$',color='#AAAAAA')
plt.plot([dt,dt],[-ma,ma],label='$\delta t$',color='#AAAAAA')
'''
plt.title('aa')
legend = plt.legend(loc='upper right', shadow=True, fontsize='x-large')


#print(x)    
#print(y[0])    

err=norm(x0-yj)/norm(x0);

tts=ttx[ttx>1*dt]
ys=jG(tts)
xs=GFt0(tts)
errs=norm(xs-ys)/norm(xs);

plt.title(sprintf('error=%3.1e%% $error_{spline}$=%3.1e%%\n',100*err,100*errs))
plt.show()
'''
errjq=norm(ynn-yjnn)/norm(ynn)
printf('err-jet-qp=%3.2e%%\n',errjq*100)

tic()
[p,d,q]=sym_reduce([c0,c1,c2,c3],[( [b1,b2,b3,b4] ,g  )]);
toc('[p,d,q]=sym_reduce([c0,c1,c2,c3],[( [b1,b2,b3,b4] ,g  )])->')
'''
sc=np.mean(np.abs(x0))/np.mean(np.abs(yj))
printf('error=%3.2e%%\n',100*err)
printf('error_s=%3.2e%%\n',100*errs)

printf('sc=%g\n 1/sc=%g \n',sc,1./sc)
