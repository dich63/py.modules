#
import os
import numpy as np
import numbers
import math
from scipy.interpolate import interp1d
from utils.derr2m import derr2m
from jsonrpc.jsonclass import *  
  
def interpolate_fun(ty1,ty2,fbf=(lambda x:x,lambda x:x ),   kind='linear'):
    def prj(x):
        if type(x) is numbers.Number:
            return x;
        else:
            return x[0];
            
    tt=[ prj(t)  for t in ty1 ];
    tt=np.array(tt,dtype=np.float64);
    
    tti=[ t[0]  for t in ty2 ];
    yyi=[ t[1]  for t in ty2 ];
    
    fi = interp1d(tti, fbf[0](yyi),kind=kind);
    
    f= lambda t : fbf[1](fi(t))
    
    
    ty=[]
    for t in tt:
        ty.append((t,f(t)));
        
    return ty;
     

class jet_scheme(object):
    def __init__(self,n):
        n2=2*n;
        nn=np.empty(n2,dtype=np.double);
        r=1.0;
        nn[0]=1;
        for k in range(1,n2):
            r*=k;
            nn[k]=r;
        self.nn=nn;
        self.bnn=1./nn;
        pv=np.poly1d(np.ones(n2,dtype=np.double));
        vv=np.zeros((n,n2),dtype=np.double);
        for k in range(0,n):
            #print k;            
            vv[k][k:n2]=np.flipud(pv.coeffs);
            pv=pv.deriv();
        self.n=n;
        self.vv=vv;    
        #self.cc=vv[:,0:n];
        #self.cb=vv[:,n:n2];
        
    def sort_pair(self,tt,ff):
        tt=np.array(tt,dtype=np.double);
        ii=np.argsort(tt);
        return (tt[ii],ff[:,ii]);
        
        
    def make_coeffs_b(self,ff):
        n,bnn=self.n,self.bnn;
        f0=ff[0:n];
        c0=f0*bnn[0:n];
        c0=np.flipud(c0);
        return c0;
        
    def make_coeffs(self,dt,ff):       
        
        
        n,nn,bnn=self.n,self.nn,self.bnn;
        
        n2=2*n;
        f0,f1=ff[0:n,0],ff[0:n,1];
        d=1;
        dtt=np.empty(n2,dtype=np.double);
        
        for k in range(n2):            
            dtt[k]=d;
            d*=dt;
            
        c0=f0*bnn[0:n];
        vv=self.vv;
        vvt=np.zeros_like(vv);

        vvt[0]=vv[0]*dtt;
        for k in range(1,n):
            vvt[k,k:n2]=vv[k,k:n2]*dtt[0:-k];

        #vvt=self.vv*dtt;
        vb=vvt[:,0:n];
        bt=f1-np.matmul(vb,c0);
        At=vvt[:,n:n2];
        c1=np.linalg.solve(At,bt);
        q=np.empty(n2,dtype=np.complex);
        q[0:n]=c0;
        q[n:n2]=c1;
        q=np.flipud(q);
        return q;
        
def isiterable(p_object):
     try:
         it = iter(p_object)
     except TypeError: 
         return False
     return True

class jet_spline(object):

     def __init__(self,tt,ff,deg=-1):
         
         n= ff.shape[0] if deg<0 else deg+1;
         
         scheme=self.scheme=jet_scheme(n);
              
         tt=np.array(tt).flatten();
         
         tt,ff=scheme.sort_pair(tt,ff);
         
         self.nt=nt=np.size(tt);
         self.tt=tt;
         self.ff=ff;         
         
         #self.vv=[None for k in range(nt+2)];
         self.vv=np.full(nt+2,None);
         
         
     def interpolate(self,t):
         #print(t)
         if isiterable(t):
             return [self.interpolate(i) for i in t ];
         
         #print(t)
         tt,ff,vv=self.tt,self.ff,self.vv;
         nt=self.nt;
         #print tt,t
         i=np.searchsorted(tt,t);
         #print('sorted i=',i)
         q=vv[i];
         if q is None:             
             if i==0:
                 q=self.scheme.make_coeffs_b(ff[:,0]);
             elif i==nt:
                 q=self.scheme.make_coeffs_b(ff[:,nt-1]);
             else:
                 dt=tt[i]-tt[i-1];
                 q=self.scheme.make_coeffs(dt,ff[:,(i-1):(i+1)]);
                 
             vv[i]=q;
         
         if  i!=0:
             dt=t-tt[i-1];
         else:
             dt=t-tt[0];
             
         return np.polyval(q,dt);
     
     def __call__(self,t):
         if isiterable(t):
             return np.array(self.interpolate(t));
         else:
             return self.interpolate(t);
             
         
def loglog_jet_spline(tt,ff,deg=3):
    nd=np.min([ff.shape[0]-1,deg,3]);
    tt=np.array(tt).flatten();
    #nt=tt.shape[0];
    s=np.log(tt);
    ff=np.array(ff,dtype=np.complex);
    lff=np.empty_like(ff,dtype=np.complex);
    z=lff[0];
    
    
    z[:]=np.log(ff[0]);

    if nd>0:
        z1=lff[1];
        z1[:]=tt*ff[1]/ff[0];    
        
        if nd>1:
            z2=lff[2];
            t2t=tt*tt;
            f2f=t2t*ff[2]/ff[0];
            z2[:]=z1*(1.0-z1)+f2f;
            
            if nd>2:
                z3=lff[3];
                z3[:]=z2-2*z1*z2+2*f2f+tt*t2t*ff[3]/ff[0] -  f2f*z1;  
           
            #z2=lff[1][:]+tt*tt*(ff[2] - ff[1]*ff[1]/ff[0])/ff[0];
            

    jsp=jet_spline(s,lff,nd);
    fu= lambda t:  np.exp(jsp(np.log(t)));
    return fu;

def log_jet_spline(tt,ff,deg=3):
    nd=np.min([ff.shape[0]-1,deg,3]);
    tt=np.array(tt).flatten();
    #nt=tt.shape[0];    
    ff=np.array(ff,dtype=np.complex);
    lff=np.empty_like(ff,dtype=np.complex);
    z=lff[0];
    
    
    z[:]=np.log(ff[0]);

    if nd>0:
        z1=lff[1];
        z1[:]=ff[1]/ff[0];    
        
        if nd>1:
            z2=lff[2];
            f2f=ff[2]/ff[0];
            zq2=z1*z1;
            z2[:]=f2f-zq2;            
            if nd>2:
                z3=lff[3];
                z3[:]=z2-2*z1*z2+2*f2f+tt*t2t*ff[3]/ff[0] -  f2f*z1;  
           
            #z2=lff[1][:]+tt*tt*(ff[2] - ff[1]*ff[1]/ff[0])/ff[0];
            

    jsp=jet_spline(tt,lff,nd);
    fu= lambda t:  np.exp(jsp(t));
    return fu;

def t_log_jet_spline(tt,ff,deg=3):
    nd=np.min([ff.shape[0]-1,deg,3]);
    tt=np.array(tt).flatten();
    #nt=tt.shape[0];    
    ff=np.array(ff,dtype=np.complex);
    lff=np.empty_like(ff,dtype=np.complex);
    z=lff[0];
    
    
    z[:]=ff[0];

    if nd>0:
        z1=lff[1];
        z1[:]=ff[1]*tt;    
        
        if nd>1:
            z2=lff[2];
            tt2=tt*tt;
            z2[:]=z1+tt2*ff[2];            
            if nd>2:
                z3=lff[3];
                tt3=tt*tt2;
                z3[:]=z1+3*tt2*ff[2]+tt3*ff[3];            
            
    lt=np.log(tt);
    jsp=jet_spline(lt,lff,nd);
    fu= lambda t:  jsp(np.log(t));
    return fu;
                             
ll_jet_spline=loglog_jet_spline;
l_jet_spline=log_jet_spline;
tl_jet_spline=t_log_jet_spline;

