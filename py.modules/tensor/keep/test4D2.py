# -*- coding: utf-8 -*-
"""
Created on Thu Feb 18 18:30:27 2021

@author: wwww
"""

from numpy import *
from utils import *
import numpy as np


from mpl_toolkits import mplot3d
#%matplotlib inline
import numpy as np
import matplotlib.pyplot as plt


from tensor.spline4D import *




def f(x, y):
    return np.sin(np.sqrt(x ** 2 + y ** 2))

Na=100
Np=14

nc=int(Np/2)

X = linspace(-2, 2, nc)
Y = linspace(-2, 2, Np)
Z = linspace(-2, 2, Np)
T = linspace(-2, 2, Np)



XX,YY,ZZ,TT = meshgrid(X,Y,Z,T)
FF=random.normal(size=[Np,Np,Np,nc]);

'''♀
plt.rcParams["figure.figsize"]=15,15
ax = plt.axes(projection='3d')
ax.plot_surface(XX[:,:,nc,nc], YY[:,:,nc,nc], FF[:,:,nc,nc], rstride=1, cstride=1,cmap='viridis', edgecolor='none');
plt.pause(0.001);
'''
tic('tensor_spline_4D init')
ts=tensor_spline_4D(FF,X,Y,Z,T,full=0)
toc(':')

X2 = linspace(-0.5, 0.5, Na)
Y2 = linspace(-0.5, 0.5, Na)

XX2,YY2 = meshgrid(X2,Y2);

tic('tensor_spline_4D make')
ff=ts(XX2,YY2,0.0,0.0);
toc(':')

ax = plt.axes(projection='3d')
ax.plot_surface(XX2, YY2, ff, rstride=1, cstride=1,cmap='viridis', edgecolor='none');

