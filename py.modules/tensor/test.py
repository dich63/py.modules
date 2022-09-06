# -*- coding: utf-8 -*-
import numpy as np
import h5py
import matplotlib.pyplot as plt
from utils import *
from tensor.spline4D import *
import tracemalloc

ff = h5py.File('FF.h5', 'r+')
ff = np.asarray(ff['FF'])
T = h5py.File('T.h5', 'r+')
T = np.asarray(T['T'])
X = h5py.File('X.h5', 'r+')
X = np.asarray(X['X'])
Y = h5py.File('Y.h5', 'r+')
Y = np.asarray(Y['Y'])
Z = h5py.File('Z.h5', 'r+')
Z = np.asarray(Z['Z'])

ts = tensor_spline_4D(ff, X, Y, Z, T, full = 0)

Nxy = 60
Nz = 10
Nt = 10

Z2 = np.linspace(Z[0], Z[-1], Nz)
X2 = np.linspace(X[0], X[-1], Nxy)
Y2 = np.linspace(Y[0], Y[-1], Nxy)
T2 = np.linspace(T[0], T[-1], Nt)
XX0, YY0, ZZ0, TT0 = np.meshgrid(X2, Y2, Z2, T2)

ff2 = ts(XX0, YY0, ZZ0, TT0)

fig, ax = plt.subplots(figsize=(10,10))
im = ax.imshow(ff[:,:,1,1])
plt.colorbar(im)

fig, ax = plt.subplots(figsize=(10,10))
im = ax.imshow(ff2[:,:,5,5])
plt.colorbar(im)
