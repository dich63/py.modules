import sys,os
sys.path.append('v:/ipc/py.modules')


os.environ['parallel_wrapper_lib']="V:/Projects/tbb_wrapper/x64/Release/tbb_wrapper.dll"
os.environ['sparse_context_lib']='O:\\__ss\\suitesparse-metis-for-windows\\bin\\x64\\release\\tpt_klu.dll'


from sympy import *
from sympy.integrals.transforms import inverse_laplace_transform
import matplotlib.pyplot as plt

import numpy as np

from lipa.solvers import *
#from parallel.LIPA import *
from lipa_ctx.LIPA_ctx import *



#AC=[np.array([1]),np.array([1])]
#sol=LIPA_solver_st(AC,dt=1);


p=Symbol('p')
x=Symbol('x')
#g=Symbol('g')


vc=4.;
q=(p**2+0.003*p+1./vc)
sAC=q*(q+0.1)*(q-0.1);


q1=(p**2+0.025*p+1./vc)
q2=(p**2+0.035*p+1./vc)


q1=(p**2+0.0025*p+1.1/vc)
q2=(p**2+0.0035*p+1.2/vc)
sAC=q*q1*q2;
sAC=(p**2+0.003*p+9);
#sAC=(p**2+0.3e-4*p+9);

g=.01+1j;

Rz=1/(sAC*p)
#Rz=1/(p**2+1e-6*p+1)
#'''
print('sG ...')
sG=inverse_laplace_transform(Rz,p,x)

sG=sG/Heaviside(x);


print(sG)
pa=Poly(sAC,p)
print(pa)
c=pa.coeffs()

fG=lambdify(x,sG);

fG1=lambdify(x,diff(sG,x,1));
fG2=lambdify(x,diff(sG,x,2));

tt=[ 0.5*k for k in  range(4000)];

yy=[fG(t) for t in tt]
yy1=[fG1(t) for t in tt]
yy2=[fG2(t) for t in tt]
#'''


'''
c=[1.10000000000000,
 0.0105000000000000,
 2.10002500000000,
 0.0100000000000000,
 1.00000000000000]
'''

c=[k for k in reversed(c)]
print ('c=',c)
#c=[1,0,1]

floatT=np.float64;
#floatT=np.complex128;
AC=[np.array([k],dtype=floatT) for k in c]


print(AC)

'''
from scipy import sparse as sp

AC[0]=sp.csc_matrix([1.0],dtype=floatT)
AC[1]=sp.csc_matrix([1.0],dtype=floatT)
'''
print('LIPA_solver_st start...')
dt=24
dt=32
dt=48
dt=43
dt=0.1
#sol=LIPA_solver_st(AC,dt=dt,fcomplex=False,pade_nm=(17,18),nd=len(AC));

sol=LIPA_solver_ctx(AC,dt=dt,fcomplex=False,pade_nm=(8,8),nd=len(AC));
sol.reset_J(1*np.array([1],dtype=floatT))

#sol.x=100;
yyn=[];
yy1n=[];
yy2n=[];

ttn=[dt*k  for k in range(int(2000./dt))];
ttn=[dt*k  for k in range(int(200./dt))];
yy=[fG(t) for t in ttn]

for n in ttn:
    y=sol.x[0];
    y1=sol.xn[1][0];
    y2=sol.xn[2][0];
    yyn.append(y)
    yy1n.append(y1)
    yy2n.append(y2)
    sol.step();


yy=np.array(yy)
yyn=np.array(yyn)
yy1=np.array(yy1)
yy1n=np.array(yy1n)
yy2=np.array(yy2)
yy2n=np.array(yy2n)


from matplotlib import rc
#rc('text', usetex=True)

fig = plt.figure();

#fig.suptitle(str(expand(sAC))+'=1',fontsize=12)
fig.suptitle(str(sAC)+'=1',fontsize=12)

fig=plt.subplot(411)

clra,clrn='#555555','#000088'
#
plt.plot(tt,yy, color=clra,label='anal')
#
plt.plot(ttn,yyn, color=clrn,label='lipa',marker='.',ls='',markersize=8)
#plt.legend()
fig=plt.subplot(412)
#
plt.plot(tt,yy1, color=clra,label='anal')
#
plt.plot(ttn,yy1n, color=clrn,label='lipa',marker='.',ls='',markersize=8)

fig=plt.subplot(413)
#
plt.plot(tt,yy2, color=clra,label='anal')
#
plt.plot(ttn,yy2n, color=clrn,label='lipa',marker='.',ls='',markersize=8)

fig=plt.subplot(414)

w0=np.sqrt(vc)*1j
ay=np.abs(yy1+w0*yy2)
ayn=np.abs(yy1n+w0*yy2n)
#

plt.plot(tt,ay, color=clra,label='anal')

#
plt.plot(ttn,ayn, color=clrn,label='lipa',marker='.',ls='',markersize=8)

#plt.legend()
print('save...')
plt.savefig('zzz.svg')
print('ok')
os.system('start chrome zzz.svg')
print(LIPA_solver_ctx)
#