# -*- coding: utf-8 -*-
"""
Created on Wed Mar 01 05:51:38 2017

@author: dich
"""
from parallel.sparse import *
from lipa.qpl import * 
import numbers
from jsonrpc.jsonclass import  jsobject,ext_def,to_dict,to_obj
from  p23 import *
import numpy as np
from scipy.linalg import expm
#import scipy as sc;
import collections

def shiftQ(QQ,t):
    pass

class qp_source_params_t(Structure):
    _pack_=8
    _fields_ = [
        ('id',c_char*16),
        ('is_active',c_longlong),
        ('n',c_longlong),
        ('b',c_void_p),
        ("z",c_double*2),
        ("degQ",c_longlong),
        ('Q',c_void_p),
        ('transform_step',c_void_p)
        ];

class qp_source_list_t(Structure):
    _pack_=8
    _fields_ = [
        ('count',c_longlong),
        ('pp_qp_source_params',c_void_p)
        ];
    
class qp_source_apply_t(Structure):
    _pack_=8
    _fields_ = [
        ('mode',c_longlong),
        ('n',c_longlong),
        ("z",c_double*2),
        ('y',c_void_p)
        ];
        
def qp_to_array(a,n,dtype):
    ta=type(a)
    if ta is numbers.Number:
        a=np.full(shape=(n,),fill_value=a,dtype=dtype);
    else:
        a=np.array(a,dtype=dtype,copy=False) 

    a=a.reshape((n,));       
    return np.ascontiguousarray(a);
    
def qp_to_dict(s,n,dt,dtype=np.complex):
    opts={
        'id':None,
        'is_active':1,
        'n':n,
        'z':0j,
        'degQ':0,
        'Q':(1,), 
        'b':1,
        'x':None,       
        'transform_step':None,
        'b_addr':0,
        'transform_step_addr':0,
        'Q_addr':0
          };
        
        
    
    ts=type(s);
    
    if not ts is dict:
        s={'b':s};

    s=to_dict(s);


    opts=to_obj(ext_def(s,opts))

    if (not opts.x is None):
        opts.b=opts.x;

    q=opts.Q=np.array(opts.Q,dtype=dtype,copy=True);
    opts.degQ=ndeg=q.size-1;

    opts.Q_addr=q.ctypes.data

    opts.b=b=qp_to_array(opts.b,n,dtype=dtype);
    opts.b_addr=b.ctypes.data;
    z=opts.z;
    
    f=(ndeg > 0) or bool(z)
    t=opts.transform_step;
    if f and (t is None):
        #p=poly_J(ndeg+1,dt=dt,dec=-opts.z)
        #opts.Udt=Udt=p.Udt;
        
        d=dt*np.ones(ndeg);
        dP=np.diag(d,1);
        Udt=expm(dP);
        Udt=np.exp(opts.z*dt)*Udt;
        
        opts.Udt=Udt
        #
        t=Udt;        
        #        t=Udt.T;
    if not t is None:
        nn=(ndeg+1)**2;        
        t=qp_to_array(t,nn,dtype=dtype);
        opts.transform_step=t
        opts.transform_step_addr=t.ctypes.data;
    return opts;

    
def create_qp_source_params(sources,n,dt):
    lo,lqp=[],[]

    sources=tolist(sources);
    qsl=qp_source_list_t();
    n_id=0;
    for s in sources:
        n_id+=1;
        o=qp_to_dict(s,n,dt);
        sp=qp_source_params_t()
        lo.append(o);        
        lqp.append(sp);
        
        sp.id=  to_bytes(n_id) if o.id is None else to_bytes(o.id);
        sp.is_active=o.is_active;
        sp.n=n;
        sp.b=o.b_addr;
        z=o.z;
        sp.z[0],sp.z[1]=z.real,z.imag;
        sp.degQ=o.degQ
        sp.Q=o.Q_addr;
        sp.transform_step=o.transform_step_addr

    ll=[addressof(p) for p in lqp]
    laddr=np.array(ll,dtype=np.uint64);
    qsl.count=len(laddr);
    qsl.pp_qp_source_params=laddr.ctypes.data;

    return (qsl,(laddr,lqp,lo,ll));
    
    
class qp_context(hcontext):
    def __init__(self,n,dt):
        hcontext.__init__(self,ptype=p_base_context_t,constructor =lambda p : create_qp_context(p))
        self.n=n;
        self.dt=dt;
    
    
    def set_sources(self,sources,dt,n,dtype=np.complex):
        dtype=np.dtype(dtype);
        p=create_qp_source_params(sources,n,dt);
        self.sourses=p[0];
        return self.create_invoke_context('set',p[0],byref=True,links=p[1]);
        
    def transform_step(self):
       return self.create_invoke_context('t')
       
    def apply_sources(self,z,y,mode=0):        
                                
        app=qp_source_apply_t();
        app.mode=mode;
        app.n=y.size
        app.z[0],app.z[1]=z.real,z.imag;        
        app.y=y.ctypes.data;
        return self.create_invoke_context('apply',app,byref=True,links=y);
 
    def get_current(self,y=None):

        if y is None:
            y=np.zeros(self.n,dtype=np.complex);
            i=self.get_current(y);
            i();
            return y;

        app=qp_source_apply_t();
        app.n=y.size;
        app.y=y.ctypes.data;
        return self.create_invoke_context('get_jc',app,byref=True,links=y);