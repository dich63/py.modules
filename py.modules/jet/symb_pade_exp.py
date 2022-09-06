# -*- coding: utf-8 -*-
"""
Created on Fri Jul 15 20:05:54 2022

@author: DICH
"""
import numpy as np
import sympy 
from sympy import *
from mpmath import *
import jsobj 

mp.dps = 48; 
mp.pretty = True


def _reverse(L):
    return [L[k] for k in range(len(L)-1,-1,-1)]

def to_number(L,dtype=complex):
    if type(L) in (tuple,list):
        return [to_number(x,dtype=dtype) for x in L ]
    else:
        return dtype(L)


def _nn_exp(L,M):
    frl=sympy.factorial
    return frl(L+M)/frl(L)
    
def _P_exp(k,L,M):
    frl=sympy.factorial
    return frl(L+M-k)*frl(L)/(frl(L+M)*frl(L-k)*frl(k))
    
def _Q_exp(k,L,M):    
    return _P_exp(k,M,L)*(-1)**k

def P_exp(L,M,norm=False,reverse=False):    
    nn=_nn_exp(L,M) if norm else 1
    r=[ nn*_P_exp(k,L,M) for k in range(L+1)]
    return _reverse(r) if reverse else r

def Q_exp(L,M,norm=False,reverse=False):    
    nn=_nn_exp(L,M) if norm else 1
    r=[ nn*_Q_exp(k,L,M) for k in range(M+1)]
    return _reverse(r) if reverse else r

def Q_roots(L,M,error=False,extraprec=55,maxsteps=150):
    q=Q_exp(L,M,norm=True,reverse=True);
    q=[np.double(k) for k in q]    
    return polyroots(q,error=error,extraprec=extraprec,maxsteps=maxsteps);
    
def exp_poles_res(LM,extraprec=55,maxsteps=150,dtype=None,scale=1):    
    L,M=LM;
    z=Symbol('z');
    
    p=P_exp(L,M,norm=True,reverse=True);
    Q=Q_exp(L,M,norm=True,reverse=True);
    pz=polyval(p,z);
    [roots,err]=Q_roots(L,M,error=True,extraprec=55,maxsteps=150)
    
    def q1(k):
        q=Q[0];
        for i in range(len(roots)):
            if i==k:
                continue
            q=q*(z-roots[i])
        return q
    
                    
    #Q1=[q1(k) for k in range(M+1) ]     
    #res=[mpmathify((pz/q1(k)).subs(z,roots[k]).simplify()) for k in range(M) ] 
    if (L,M)==(0,1):
        res=[mpc(-1)/scale];
        roots=[mpc(1)/scale];
    else:
        res=[(pz/q1(k)).subs(z,roots[k]) for k in range(M) ]     
        roots=[r/scale for r in roots]
        res=[mpmathify((r/scale).simplify()) for r in res]
    #res=[(r/scale).simplify() for r in res]
    
        
    
    return (roots,res,err) if (dtype is None) else (np.array(to_number(roots,dtype)),np.array(to_number(res,dtype)),np.double(err))


def exp_poles_res_dt_ZD(LM,D,dt=1,extraprec=55,maxsteps=150,dtype=complex):
    
    
    def to_a(x):        
        return np.array(to_number(x,dtype=dtype),dtype=dtype)
    
    
    poles1,res1,err=exp_poles_res(LM,extraprec=extraprec,maxsteps=maxsteps);
    
    poles,res=[r/dt for r in poles1],[r/dt for r in res1]
    O=jsobj.jsc(poles=poles,res=res,poles1=poles1,res1=res1)
    R=jsobj.jsc(poles=to_a(poles),res=to_a(res))
    
    dt=Number(dt);
    resD=[];
    ZD=[];
    #for 
    #rd=[r for r in res]
    
    
    M=len(poles);
    pd=[p for p in poles];
    ZD=[[Number(1) for m in range(M)]]
    #pd=[Number(1) for m in range(M)];
    for d in range(D):
        rd=[res[m]*pd[m] for m in range(M)];
        resD+=[rd]
        ZD+=[pd]
        pd=[poles[m]*pd[m] for m in range(M)];
        
        
    O.resD=resD;
    O.ZD=ZD;
    R.resD=to_a(resD);
    R.ZD=to_a(ZD);
    return R,O
    
class exp_poles_res_cache_t(object):
    def __init__(this):
        import tempfile
        this._cache={};
        this.fn=tempfile.tempdir+'/~mpmath_dps['+str(mp.dps)+']_exp_poles_res_dps=.~cache'
        
    def io(this,rw=False):
        
        try:
            import dill
            dill.settings['recurse'] = True
            if not rw:
                this._cache=dill.load(open(this.fn,'rb'))
            else:            
                dill.dump(this._cache,open(this.fn,'wb'))
        except:
            pass;
            
    def reset(this):
        this._cache={};
        this.io(True);
        
    def  get_LM(this,LM):
        LM=tuple(LM)
        try:
            return this._cache[LM];
        except:
            try:
                this.io();
                return this._cache[LM];
            except:
                r=this._cache[LM]=exp_poles_res(LM);
                this.io(True);
                return r;
            pass
    def __call__(this,LM):
        return  to_number(this.get_LM(LM));
        
exp_poles_res_cache=exp_poles_res_cache_t()    

def pade_roots_error(LM):
    p,r,e=exp_poles_res_cache.get_LM(LM);
    return np.float64(e);
    
def pade_exp_poles_res(LM,t=1,fhalf=False,dtype=complex):
    
    def to_a(x):        
        return np.array(to_number(x,dtype=dtype),dtype=dtype)
    
    p,r,e=exp_poles_res_cache.get_LM(LM);
    
    if fhalf:        
       pf,rf=p,r
       p,r=[],[]
       
       for m in range(LM[1]):           
           pi=complex(pf[m]).imag
           if pi<-1e-7:
               continue           
           sc= 2 if np.abs(pi)>1e-4  else 1;           
           p+=[pf[m]/t]
           r+=[sc*rf[m]/t]
    else:
        p=[v/t for v in p]
        r=[v/t for v in r]
        
    #print('LM=',LM)
           
    return to_a(p),to_a(r)

#q=exp_polus_res([3,3])    
if __name__=='__main__':
    
    from utils import *
    #tic();r,o=exp_poles_res_dt_ZD([8,8],3,0.001);toc(':')
    #p,r,e=exp_poles_res([4,7])
    p,r=pade_exp_poles_res([12,13],t=1,fhalf=1);
    '''
    for M in range(1,21):
        for L in range(M+1):
            tic()
            exp_poles_res_cache([L,M])
            toc(sprintf('[%d,%d]: ',L,M))
    '''