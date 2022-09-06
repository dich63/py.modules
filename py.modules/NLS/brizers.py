# -*- coding: utf-8 -*-
"""
Created on Thu Mar  4 17:09:14 2021

@author: wwww
"""

import numpy as np


def brizerKM_periods(a1):    
    a1=a1+0j;
    p=np.sqrt(2-4*a1);
    pi2=(2.0*np.pi);
    bb=np.sqrt(8*a1*(1-2*a1));
    return (pi2/p,1j*pi2/bb)

def brizerKM(t,z,a1=1/4):
    
    a1=a1+0j;
    bb=np.sqrt(8*a1*(1-2*a1));
    p=np.sqrt(2-4*a1);
    cospt=np.cos(p*t);
    chbz=np.cosh(bb*z);
    shbz=np.sinh(bb*z);
    
    sqrt2a1=np.sqrt(2*a1);
    
    N=(1-4*a1)*chbz+sqrt2a1*cospt+1j*bb*shbz;
    D=(chbz-sqrt2a1*cospt);
    
    F=(N/D)*np.exp(1j*z);
    return F;


def rbrizerKM(N,):
    pass

if __name__=='__main__':
    brizerKM_periods(0.71)
