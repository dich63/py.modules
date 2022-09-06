# -*- coding: utf-8 -*-
"""
Created on Wed Jan 27 00:12:00 2016

@author: dich6_000
"""
#import copy
import numpy as np
from numpy import linalg as la

"""
def rect2quad(r):
    xb,yb,xe,ye=r
    q=[[xb,yb],[xb,ye],[xe,ye],[xe,yb]]
    return q
"""

def rect_xxyy2quad(r):
    xb,xe,yb,ye=r
    q=[[xb,yb],[xb,ye],[xe,ye],[xe,yb]]
    return q

def rect_xyxy2quad(r):
    xb,yb,xe,ye=r
    q=[[xb,yb],[xb,ye],[xe,ye],[xe,yb]]
    return q

rect2quad=rect_xxyy2quad


def quad2quad_projective(x,y):
    """
       | a b c | xn[0] = yn[0]
G--->  | d e f | xn[1] = yn[1]
       | g h 1 |   1   = Qn
     
     x[n][0]*a+x[n][1]*b+c=(x[n][0]*g + x[n][1]*h+1)*y[n][0]
     x[n][0]*d+x[n][1]*e+f=(x[n][0]*g + x[n][1]*h+1)*y[n][1]
     
     x[0][0] | x[0][1] | 1 |      0  |      0  | 0 | -x[0][0]*y[0][0]  -x[0][1]*y[0][0] =y[0][0]
          0  |      0  | 0 | x[0][0] | x[0][1] | 1 | -x[0][0]*y[0][1]  -x[0][1]*y[0][1] =y[0][1]
          
     x[n][0] | x[n][1] | 1 |      0  |      0  | 0 | -x[n][0]*y[n][0]  -x[n][1]*y[n][0] =y[n][0]
          0  |      0  | 1 | x[n][0] | x[n][1] | 1 | -x[n][0]*y[n][1]  -x[n][1]*y[n][1] =y[n][1]
     """
    GG=()
    B=()
    for n in range(4):
        GG+=((x[n][0],x[n][1],1,0,0,0,-x[n][0]*y[n][0],-x[n][1]*y[n][0])
            ,(0,0,0,x[n][0],x[n][1],1,-x[n][0]*y[n][1],-x[n][1]*y[n][1]))
        B+=(y[n][0],y[n][1])
    B=np.array(B,dtype='d')
    GG=np.array(GG,dtype='d')
    G=la.solve(GG,B);
    
    TT=((G[0],G[1],G[2]),
        (G[3],G[4],G[5]),
        (G[6],G[7],1.0))
    
    return TT;

def projmul(TT,x):
    if len(x)==4:
        return [projmul(TT,x[0]),projmul(TT,x[1]),projmul(TT,x[2]),projmul(TT,x[3])]
    
    x=(x[0],x[1],1.0)
    y=np.matmul(TT,x);
    by=1./y[2];
    return [by*y[0],by*y[1]]


