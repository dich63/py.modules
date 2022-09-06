# -*- coding: utf-8 -*-
"""
Created on Tue Jun 28 14:45:20 2022

@author: DICH
"""
from fs1D_model import *




model=':file:'+__file__+'/../mod0.json'
model=':file:'+__file__+'/../model1D.json'

fs=fs1D_t(model);
fs.solve(dt=1/100,dz=0.1,Tm=1)
    