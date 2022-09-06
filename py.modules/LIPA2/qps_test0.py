# -*- coding: utf-8 -*-
"""
Created on Mon Mar 28 21:11:43 2022

@author: wwww
"""
from sympy import *
from sympy.matrices import Matrix, eye, zeros, ones, diag, GramSchmidt
from sympy.abc import x

x0=Symbol('x0');x1=Symbol('x1');x2=Symbol('x3');x3=Symbol('x3');#f=Symbol('f');

g=Symbol('g')
c=Symbol('c')

G=Symbol('G')

f=Function('f')
s=(x**2)*exp(-G*x)
#rq=dsolve(Derivative(f(x),x,x,x)10+Derivative(f(x),x,x)+9*f(x),f(x),\
#res=dsolve((Derivative(f(x),x,x)+g*Derivative(f(x),x)+(c**2)*f(x)+s),f(x),\
res=dsolve((f(x).diff(x,2)+g*Derivative(f(x),x)+(c**2)*f(x)+s),f(x),\
ics={f(0):x0,\
     f(x).diff(x,1).subs(x,0):x1\
         })

r=res.rhs



Gn=1+16j;
rgc=r.subs([[g,0.01],[c,3]])
rgc=r.subs([[g,0.01],[c,3],[G,Gn]])
j=lambdify(x,s.subs(G,Gn))
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

