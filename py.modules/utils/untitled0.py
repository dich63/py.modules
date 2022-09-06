# -*- coding: utf-8 -*-
"""
Created on Tue Dec  7 16:43:25 2021

@author: wwww
"""

import k3d
import numpy as np

Nx, Ny  = 34, 33
xmin, xmax = -3, 4
ymin, ymax = -0, 3

x = np.linspace(xmin, xmax, Nx)
y = np.linspace(ymin, ymax, Ny)
x, y = np.meshgrid(x, y)
f = np.sin(x**2 + y**2)

plot = k3d.plot()
plt_surface = k3d.surface(f.astype(np.float32), bounds=[xmin,xmax,ymin,ymax])
plot += plt_surface
plot.display()