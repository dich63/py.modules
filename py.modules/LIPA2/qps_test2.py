# -*- coding: utf-8 -*-
"""
Created on Mon Mar 28 21:11:43 2022

@author: wwww
"""


from sympy import *
from sympy.matrices import Matrix, eye, zeros, ones, diag, GramSchmidt
from sympy.abc import x
#from sympy.functions.special.hyper import *
from sympy.functions.special.hyper import *
from mpmath import *
from sympy.integrals.transforms import inverse_laplace_transform

mp.dps = 55; 

x0=Symbol('x0');x1=Symbol('x1');x2=Symbol('x2');x3=Symbol('x3');#f=Symbol('f');
p=Symbol('p');p1=Symbol('p1');p2=Symbol('p2');
p3=Symbol('p3');p4=Symbol('p4');p0=Symbol('p0');
y=Symbol('y');
d=Symbol('d');
z=Symbol('z');

c1=Symbol('c1');c2=Symbol('c2');c3=Symbol('c3');c4=Symbol('c4');
g1=Symbol('g1');g2=Symbol('g2');
G=Symbol('G');w=Symbol('w');
J1=Symbol('J1');J3=Symbol('J3');J2=Symbol('J2');J0=Symbol('J0')

t=Symbol('t')
q=series((x**2+g1*x*y+c1**2*y**2)*(x**2+g2*y*x+c2**2*y**2),x)
#q=series(x**2+g2*y*x+c2**2*y**2,x)
q1=q.subs(y,1)
r1=1/q.subs(y,1)
rr=r1.subs([[c1,2],[c2,2+0.1],[g1,1/30],[g2,1/29]])

rr=r1.subs([[c1,1],[c2,1+0.1],[g1,1/30],[g2,1/29]])

qr=q1.subs([[c1,1],[c2,1],[g1,1/3],[g2,1]])

#e=inverse_laplace_transform(rr/(x+8)**4,x,z)

from  LIPA2.ilaplace import ilaplace,fast_inverse_laplace
from utils import *



qq=rr/(x+3+7j)**3
pp=(0/(x+3+7j)**2+1/(x+3+7j)**3)
pp=simplify(70/(x+3-7j)**2+1/(x+3-7j)**3)
pp=simplify(70/(x+1e-3)**3)


#pp=simplify(70/(x+3-7j)+1/(x+5-17j))
#pp2=simplify(7/(x+4+3j)**2)
#pp=7/(x+3+7j)
#pp=7/x**1
qqp=simplify(qq*pp)
tic()

#w=fast_inverse_laplace(qq,x,t)
#d=InverseLaplaceTransform(exp=qq,s=x,plane=1)
ju,jif=ilaplace(pp ,x,5)
fu,iif=ilaplace(qqp ,x,5)

#fu,iif=ilaplace(rr/((x+4)*(x+8)**2) ,x)
toc('create ilaplace')

tic()
v=fu([1,2,3,4,5,6,7,8,9,10])
toc('calc')
print(v)

import matplotlib.pyplot as plt
#%matplotlib widget
fig=plt.figure(figsize=(12,8))
t=arange(0,300,0.5)
tic()
jt=ju(t);
toc('plot j:')
ft=fu(t);
toc('plot x:')
import numpy as np
nr= lambda f: f/np.max(np.abs(np.array(f).real))
#plt.plot(t,nr(np.abs(jt)),t,nr(np.abs(ft)))
plt.plot(t,nr(jt),t,-nr((ft)))

raise SystemExit()

#fua=lambdify(x,e,modules = ["numpy", "scipy","mpmath"])
fua=lambdify(x,e,modules = ["sympy","mpmath"])
fua(11.1)

from utils import *
q1=q.subs([[x**4,p4]]).subs([[x**3,p3]]).subs([[x**3,p3]]).subs([[x**2,p2]]).subs([[x,p1]])


q1=q1.subs([[y**4,p0]])
q1=q1.subs([[y,1]])
#q1=q1.subs([[p4,d*p4]])

coe=[ q1.coeff('p'+str(k)) for k in range(4,-1,-1)]

f=Function('f')


sr=x#*exp(-16*x)#*sin(w*x)♥
sr=exp(-G*x)

sr=J2*x**2+J1*x+J0

qq=q1.subs(p0,f(x).diff(x,0)).subs(p1,f(x).diff(x,1)).subs(p2,f(x).diff(x,2)).subs(p3,f(x).diff(x,3)).subs(p4,f(x).diff(x,4))

tic()


res=dsolve(qq+sr,f(x),ics={
    f(0):x0,
    f(x).diff(x,1).subs(x,0):x1,
    f(x).diff(x,2).subs(x,0):x2,
    f(x).diff(x,3).subs(x,0):x3
    },simplify=False);

'''

res=dsolve(qq+sr,f(x),ics={
    f(0):0,
    f(x).diff(x,1).subs(x,0):0,
    f(x).diff(x,2).subs(x,0):0,
    f(x).diff(x,3).subs(x,0):0
    },simplify=False);
'''


toc('dsolve end:')


'''
g=Symbol('g')
c=Symbol('c')

f=Function('f')
s=(x**1)*exp(-x)*sin(16*x)
#rq=dsolve(Derivative(f(x),x,x,x)10+Derivative(f(x),x,x)+9*f(x),f(x),\
#res=dsolve((Derivative(f(x),x,x)+g*Derivative(f(x),x)+(c**2)*f(x)+s),f(x),\
res=dsolve((f(x).diff(x,2)+g*f(x).diff(x,1)+(c**2)*f(x)+s),f(x),\
ics={f(0):x0,\
     f(x).diff(x,1).subs(x,0):x1\
         })

r=res.rhs




rgc=r.subs([[g,0.01],[c,3]])
j=lambdify(x,s)
fua=lambdify((x,x0,x1),rgc)    
print(fua(1,2,3))    

import numpy as np
#%matplotlib inline
import matplotlib.pyplot as plt
#%matplotlib widget
fig=plt.figure(figsize=(12,8))
#fig.suptitle('sensors data + pure response ', fontsize=16)
plt.grid(True,which='major')
plt.minorticks_on()
plt.grid(True,which='minor',alpha=0.2)
tt=np.linspace(0,20,2000)
plt.plot(tt,fua(tt,0,0),tt,j(tt)/400)

import dill
dill.settings['recurse']=True
dill.dump(r, open("res_rhs.dill", "wb"))
'''
