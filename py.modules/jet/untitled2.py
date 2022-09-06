# -*- coding: utf-8 -*-
"""
Created on Tue Jul 19 00:46:30 2022

@author: DICH
"""

import numpy as np
import sympy 
from sympy import *
from mpmath import *
from utils import *
import jsobj 
mp.dps = 36; 
mp.pretty = True

s=np.random.rand(512,512)
s=np.random.rand(256,256)
print('matrix:',s.shape)
m=matrix(s);
e=matrix(np.eye(s.shape[0]));
tic();mz=m-1j*e;toc('mz=m-1j:');
tic();im=mz**-1;toc('im=mz**-1;');