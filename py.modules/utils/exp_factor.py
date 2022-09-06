# -*- coding: utf-8 -*-
"""
Created on Fri Sep  2 02:40:09 2022

@author: DICH
"""

import numpy as np
# x_n ~ A*exp(G*n)
def exp_factor(x,dt=1.0):        
    x=np.array(x);
    x0=np.mean(x)    
    x/=x0;
    x=x[x>0];
    N=len(x)-1;      
    x=np.log(x);
    a=N*(N+1)/2.0
    b=N+1
    c=N*(N+1)*(2*N+1)/6.0
    d=a;
    det=a*d-b*c;
    sx0=1*np.sum(x)
    nn=np.arange(N+1)
    sx1=1*np.sum(x*nn);
    r=[[d,-b],[-c,a]]@np.array([sx0,sx1]);
    r/=det;
    g=r[0]/dt;
    A=np.exp(r[1]);
    xm=np.mean(x);
    d=np.linalg.norm(x-r[0]*nn-r[1])/np.abs(xm)    
    return g,A*x0,d;

if __name__=='__main__':  
    
    from utils import *
    
    x=-10.001*np.exp(-1.3334*np.arange(111));
    print('exp_factor(x)=',exp_factor(x))
    r=np.random.randn(len(x))
    x+=1.4*r*x;
    print('exp_factor(x+r)=',exp_factor(x))
    
    