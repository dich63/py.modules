# -*- coding: utf-8 -*-
"""
Created on Tue Mar 16 20:54:06 2021

@author: wwww
"""

import numpy as np

x=np.arange(12).reshape(3,4)
y=x.copy()
for i, xi in enumerate(x):
    print(i,xi);
    y[:,i]=xi