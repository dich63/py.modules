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
plt.grid(True,which='minor',alpha=0.8)
'''



'''
import numpy as np
norm=np.linalg.norm

syms_def(locals(),' c0, c1 ,c0,c1,c2,c3,c4,b0,b1, b2, b3, b4,a0, g, a1, a2, a3, a4, g1')

tic();

[p,d,q,sq]=sym_reduce_ex([c0,c1],[],[a0,a1,a2]);
toc('reduce:')
#jo=arg2jso(c1=1.1,c2=2.3,c3=1.11,c4=21,c0=10j,g=3.11,g1=2j,b1=2.1,b2=1.2,b3=1.1j,a4=2)



jo={
    
    c0:-11.5+125j,    
    c1:-11.5-225j,    
    a0:1,
    a1:0,
    a2:0
    }






tic();
nd=0

t0=0.02
sc=0

LM=[4,12]
LM=[2,5]
dt=0.01
LM=[8,18];
dt=0.0225*3.8/2;

#dt=0.0225*1.0/2;
#dt=1
'''
LM=[4,8];
dt=0.00225*1;


'''

t0=2*dt/4
#t0=0

[GF0,exg]=ilaplace_functor(p,jo,nd=nd)

GF =lambda t: GF0(t-t0)
#[GF1,exg1]=ilaplace_functor(p,jo,nd=nd+1)
#[GF2,exg2]=ilaplace_functor(p,jo,nd=nd+2)

toc('functors:')


tR=1
tt=np.linspace(0,tR,12000)
rn=lambda x: x.real/np.max(np.abs(x.real))
#rn=lambda x: x.real
#plt.plot(tt,rn(F(tt)),tt,rn(GF(tt)))
#plt.plot(tt,rn(F(tt)-F(tt-t0)),tt,rn(GF(tt)-GF(tt-t0)))
#plt.plot(tt,(F(tt)-F(tt-t0)),tt,F(tt))



%matplotlib auto
#%matplotlib inline

#plt.plot(tt,rn(F(tt)),tt,rn(Ft0(tt)))
'''
plt.plot(tt,rn(GF(tt)))
plt.plot([t0],[0], marker='o',linewidth=0)
plt.plot([dt],[0], marker='x',linewidth=0)
plt.show()
'''


[DC,QP,G,sing,FF,FFS]=make_lipa_data_ex(d,q,[sq],jo);



'''
print(QP.shape)
print(G.shape)
print(FF.shape)
'''
D=DC.size-1
D=10

#poles,res=pade_exp_poles_res(LM,t=dt);

#QP=np.array([[0j,0]]);
#QP=[[0]]
cQP=QP.copy()
QP2=QP.copy()



#lqp=lipa_qp_number(DC,LM=LM,FF=FF,qp=cQP,g=G,nd=D+1).reset(dt);
FF2=np.array([1+0j,1]).T

jet=jet_csr_klu_number(DC,D=D).reset(dt=dt,LM=LM);
jet0=jet_csr_klu_number(DC).reset(dt=dt,LM=LM);

ssc=source_singular_csr(sing,t0=[t0])

#
jet.reset_singular_source(ssc,FFS)
#jet.reset_singular_source(ssc,FFS[:,0:1])

'''
jet0.xx[0]=-1
jet0.xx[1]=1


jet()
jet0()
jet()
jet0()
'''

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
tt=dt*np.arange(1,nt+1)
tm=np.max(tt)
ttx=np.linspace(dt,tm,10000)
ttx0=np.linspace(0,tm,10000)
x=GF(tt)
#x1=GF1(tt)
#x2=GF2(tt)

#x=GF(tt)-GF(tt-t0)
tic('start yjnn=jet.dump(nt)')
yjnn=jet.dump(nt)
toc(' yjnn=jet.dump(nt)')

yjnn=yjnn.reshape(yjnn.shape[0],-1).T


yj=yjnn[nd]

jG=jet_spline(tt,yjnn,D)
fu= lambda x: np.abs(x)
fu=lambda x: x.real/np.max(np.abs(x))
fu= lambda x: x
#fu=lambda x: np.abs(x)/np.max(np.abs(x))
#fu= lambda x: x


#plt.plot(ttx,fu(GF(ttx)),ttx,fu(jG(ttx)))
#plt.plot(ttx,fu(GF(ttx).imag),ttx,fu(jG(ttx).imag))
#plt.plot(ttx,fu(GF0(ttx))); 

plt.plot(ttx0,fu(GF(ttx0)),label='$G^{ex}_{real}$'); 
plt.plot(ttx,fu(jG(ttx)),label='$spline_{real}$')

plt.plot(ttx0,fu(GF(ttx0).imag),label='$G^{ex}_{imag}$')
plt.plot(ttx,fu(jG(ttx).imag),label='$spline_{imag}$')
#plt.plot(tt,fu(y), marker='o',linewidth=0)
plt.plot(tt,fu(yj), marker='o',linewidth=0,label='$jet_{real}$')
plt.plot(tt,fu(yj.imag), marker='x',linewidth=0,label='$jet_{imag}$')
plt.plot([t0,dt],[0,0], marker='+',linewidth=0)
'''
'''
legend = plt.legend(loc='upper right', shadow=True, fontsize='x-large')
plt.grid(True,which='major')
plt.minorticks_on()
plt.grid(True,which='minor',alpha=0.8)

plt.show()
#print(x)    
#print(y[0])    

errj=norm(x-yj)/norm(x);
printf('errorj=%3.2e%%\n',100*errj)
gfj=GF(ttx);
spj=jG(ttx);
errsp=norm(gfj-spj)/norm(gfj);
printf('error_spline=%3.2e%%\n',100*errsp)

