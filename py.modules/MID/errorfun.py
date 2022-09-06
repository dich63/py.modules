# -*- coding: utf-8 -*-
"""
Created on Tue Apr 05 17:24:55 2016

@author: dich6_000
"""
import numpy as np

def errfun(uet,u, sig, t):
    uet = copy.deepcopy(uet)
    u = copy.deepcopy(u)
    me = np.mean(uet)
    m = np.mean(u)
    uet = uet/me
    u=u/m
    #e=norm(derr2m(uet,u),inf);
    #e=norm(derr2m(uet,u),1);
    e=np.linalg.norm(derr2m(uet,u))
    a=me/m

    # dt = np.diff(t, n=1, axis=0)
    # dt = np.append( dt, dt[-1] )
    # w = np.sqrt(dt)/np.sqrt( np.power(uet, 2) + sig**2 )
    # v  = w*u
    # a = np.dot(v,uet)/np.dot(v,u)
    # eps = np.abs(uet - a*u) / np.sqrt(np.power(uet, 2) + sig**2)
    # e = np.dot(eps, eps)
    return e, a

def derr2m(x,y):
    epsa=2^-18
    x,y=np.array(x).flatten(),np.array(y).flatten()
    axx = np.abs(y) + np.abs(x)
    axx[np.where(axx<epsa)[0]]=1
    xy=y-x
    xy[np.where(abs(xy)<epsa)[0]]=0
    err=xy/(axx+epsa)
    err[np.where(abs(err)>1)[0]] = 1
    return err