# -*- coding: utf-8 -*-
"""
Created on Tue Dec 22 02:09:18 2015

@author: dich6_000
"""
import numpy as np
import numbers



def extend_exists(d1,d2=None):
    d={}
    for k in d1:
        d[k]=d2.get(k,d1[k])
    return d

def extend(d1,d2=None):
    d={}
    if d2:
        d.update(d2);
    if d1:
        d.update(d1);
    return d

def tolist(l):
    if (type(l)==tuple) or (type(l)==list):
        return l;
    else:
        return (l,);

def calc_Cx1(xc,z,CC,xx,x_tmp):
    l=len(CC);
    if l>1:
        if l>len(xx):
            raise Exception('CC>xx')
        x_tmp[:]=xx[0];
        xf=x_tmp;
        for k in range(1,l):
            xf*=z;
            xf+=xx[k];
            xc+=CC[k]*xf;
        
    return xc;

def calc_Cx0(xc,CC,xx):
    if CC:
        xc[:]=CC[0]*xx[0]
    else:
        xc[:]=xx[0]
    return xc;

def calc_Cx0_fake(xc,CC,xx):
    xc[:]=0;
    return xc;
"""
def calc_Cx(xc,z,CC,xx,x_tmp):
    return calc_Cx1(calc_Cx0(xc,CC,xx),z,CC,xx,x_tmp)
"""
def calc_Cx(xc,z,CC,xx,x_tmp):
    l=len(CC);
    if l>0:
        if l>len(xx):
            raise Exception('CC>xx')

        x_tmp[:]=xx[0];
        xf=x_tmp;
        xc+=CC[0]*xf;
        for k in range(1,l):
            xf*=z;                        
            xf+=xx[k];
            xc+=CC[k]*xf;
    """            
        x_tmp[:]=0;
        xf=x_tmp;
        for k in range(0,l):            
            xf+=xx[k];
            xc+=CC[k]*xf;
            xf*=z;                        
    """      
            
        
    return xc;


class poly_J(object):
    def __init__(self,ndeg,dt=1,dec=0):
        n=ndeg;
        nn=np.zeros((n,1),dtype=np.double);
        p=1;
        nn[0]=1;
        for k in range(1,n):
            p*=dt/np.float64(k);
            nn[k]=p;

        #nn=nn.reshape((1,n));
        self.dec=dec;

        Udt=np.matrix(np.zeros((n,n),dtype=np.complex128))
        for k in range(n):
            #Udt[k,k:n]=nn[0,0:n-k];
            Udt[k:n,k]=nn[0:n-k];
        Udt*=np.exp(-dec*dt);
        self.Udt=Udt;

    def renorm(self,c):
        r=c*self.Udt;
        return r.tolist()[0];

class poly_Jc(poly_J):
    def __init__(self,cJ,dt=1,dec=0):
        ndeg=len(cJ);
        self.cJ=cJ;
        poly_J.__init__(self,ndeg,dt,dec);
    def renorm(self):
        cJ = super(poly_Jc,self).renorm(self.cJ);
        return cJ;
    def update(self):
        self.cJ=self.renorm();
        return self.cJ
    def laplace(self,z):
        bz=np.complex128(1)/(z+self.dec);
        c=np.flipud(self.cJ)
        return bz*np.polyval(c,bz);


class poly_Jcs:
    def __init__(self,Js,dt=1):
        pjs=[];
        for j in Js:
            d={'c':(1,),'dec':0};
            d.update(j);
            pjs.append(poly_Jc(d['c'],dt,d['dec']));
        self.pjs=pjs;

    def update(self):
        for p in self.pjs:
            p.update();
    def laplace(self,z):
        v=np.complex128(0);
        for p in self.pjs:
            v+=p.laplace(z);
        return v;
    def j0(self):
        v=np.complex128(0);
        for p in self.pjs:
            v+=p.cJ[0];
        return v;

class qpJ:
    def __init__(self,qp=None,dt=1):
        jl=[];
        xx=[];
        if qp:
            if type(qp)==dict:
                qp=(qp,);
            for j in qp:
                q=j.get('qp',{})
                #x0=j.get('x',1);
                x0=j['x'];
                jl.append(poly_Jcs(q,dt));
                xx.append(x0.flatten());
        self.jl=jl;
        self.xx=xx;
    def set_J(self,z,xb=0):
        #if not isinstance(xb,numbers.Number):
        #    xb[:]=0;
        ix=iter(self.xx)
        for j in self.jl:
            x=next(ix);
            xb+=j.laplace(z)*x;
            j.update();
        return xb;


