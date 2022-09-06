# -*- coding: utf-8 -*-
"""
Created on Tue Apr 05 17:24:55 2016

@author: dich6_000
"""
import numpy as np


def derr2m(x,y,epsa=1e-18):
    
    x,y=np.array(x),np.array(y)
    axx = np.abs(y) + np.abs(x)
    axx[axx<epsa]=1
    xy=y-x
    xy[abs(xy)<epsa]=0
    err=xy/(axx)
    err[abs(err)>1] = 1
    return err

