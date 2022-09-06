# -*- coding: utf-8 -*-
"""
Created on Tue Jun  1 14:52:22 2021

@author: wwww
"""

from numpy import *
from utils import *

import numpy as np
from mpl_toolkits import mplot3d
import matplotlib.pyplot as plt


from tensor.spline4D import *

import tracemalloc

tracemalloc.start()

def tracedump():
    current, peak = tracemalloc.get_traced_memory()
    
    print(f"Current memory usage is {round(current / 2**30,3)} GB; Peak was {round(peak / 2**30,3)} GB")


Na=100
Np=50
Np2=int(Np/2)
X = Y = Z = T = linspace(-2, 2, Np)
Npx,Npy,Npz,Npt=15,25,23,35

X , Y , Z , T =[ linspace(-2, 2, Nk) for Nk in (Npx,Npy,Npz,Npt)  ]

#X= linspace(-2, 2, Np2)
#nc=int(Np/2)

nc=20

XX,YY,ZZ,TT = meshgrid(X,Y,Z,T)
#FF=random.normal(size=[Np2,Np,Np,Np]);

FF=random.normal(size=[Npx,Npy,Npz,Npt]);


tracedump()

tic('tensor_spline_4D init')
ts=tensor_spline_4D(FF,X,Y,Z,T,full=0)
toc(':')
tracedump()
d=2-1e-5
X2 = Y2 = linspace(-d,d,Na)
XX0,YY0 = meshgrid(X2,Y2);

tic('tensor_spline_4D  [pass 1] make  ')
ff0=ts(XX0,YY0,0.01,0.01);
toc(':')
tracedump()
print('cache_count=',ts.cache_count)

X2 = Y2 = linspace(-0.5, 0.5, Na)


XX2,YY2 = meshgrid(X2,Y2);


tic('tensor_spline_4D [pass 2] make')
ff=ts(XX2,YY2,0.,0.);
toc(':')
print('cache_count=',ts.cache_count)
tracedump()
tracemalloc.stop()

# ---- Graphics.....
# ---- Graphics.....
# ---- Graphics.....
plt.rcParams["figure.figsize"]=15,15
ax = plt.axes(projection='3d')
ax.plot_surface(XX[:,:,nc,nc], YY[:,:,nc,nc], FF[:,:,nc,nc].T, rstride=1, cstride=1,cmap='viridis', edgecolor='none');

plt.pause(0.001);
ax = plt.axes(projection='3d')
ax.plot_surface(XX0, YY0, ff0, rstride=1, cstride=1,cmap='viridis', edgecolor='none');

plt.pause(0.001);
ax = plt.axes(projection='3d')
ax.plot_surface(XX2, YY2, ff, rstride=1, cstride=1,cmap='viridis', edgecolor='none');