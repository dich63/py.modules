# -*- coding: utf-8 -*-
"""
"""
from NLS.nls_pade import *

def nls_rescale_to(tt,sname='NLS',g=1,nm=[4,4]):
    tt=np.array(tt,dtype='d').flatten();
    N=tt.size;
    h=(tt[-1]-tt[0])/(N-1);
    h2=h*h;
    gr=h2*g;
    dzu=1.0/h2;
    solver= NLS_ex(sname,N,g=gr,dt=dzu,nm=nm);    
    #solver.dz=dzu;
    #solver.h=h;
    return (solver,dzu,h); 

def js_nls_rescale_to(tt,sname='NLS',g=1,nm=[4,4]):
    from jsobj import jso
    o=jso();
    (o.solver,o.dzu,o.h)=nls_rescale_to(tt=tt,sname=sname,g=g,nm=nm)
    return o; 


def fiber_rescale_to(N,winT=None,fmax=1.0,beta=1.0,gamma=1.0,alpha=0.0):
    
    N=int(N);
    
    if winT is None:
        dT=0.5/fmax;
    else:
        dT=np.float64(winT)/N;
        fmax=0.5/dT;  
        
    rs=beta/(dT**2);    
    g_r=gamma/rs;
    alpha_r=alpha/rs; 
    dzu=rs;        
    return (dzu,rs,g_r,alpha_r,dT,N,fmax);


def DBpKM2Alpha(dB_per_km=0.025):
    return (np.log(10)/20)*dB_per_km;

    
def fiber_rescale_to_jso(N,winT=None,fmax=1.0,beta=1.0,gamma=1.0,alpha=None,alpha_DBpKM=0.0):
    
    if alpha is None:
        alpha=DBpKM2Alpha(alpha_DBpKM); 
        
    from jsobj import jso
    o=jso();
    (o.dzu,o.rs,o.g_r,o.alpha_r,o.dT,o.N,o.fmax)=fiber_rescale_to(N,winT=winT,fmax=fmax,beta=beta,gamma=gamma,alpha=alpha);
    
    return o;


def manakov_solver_create_from_jso(name,o,nm=[4,4]):
    
    def parse_name(f):
        if type(f)==str:            
            return {'nls':False,'pade':False,'ssf':True}[f.lower()];
        else:
            return not not f;    
    
    solver=MS(o.N,g=o.g_r,dt=o.dzu,alpha=o.alpha_r,nm=nm,fSSF=parse_name(name));    
    return solver
    



def pow_norm(x,dT=1.0,energy=1.0):
    x=np.array(x,copy=False);
    Z=np.sqrt(energy/dT)/scipy.linalg.norm(x);
    return Z*x;

def nls_fiber_rescale_to(N,sname='NLS',fmax=1.0,beta=1.0,g=1.0,alpha=0.0,nm=[4,4]):
    
    (dzu,rs,g_r,alpha_r,dT)=fiber_rescale_to(N,fmax=fmax,beta=beta,g=g,alpha=alpha)
    solver= NLS_ex(sname,N,g=g_r,dt=dzu,alpha=alpha_r,nm=nm);
    return (solver,dzu,rs,g_r,alpha_r,dT);



def manakov_fiber_rescale_to(N,fSSF=False,fmax=1.0,beta=1.0,g=1.0,alpha=0.0,nm=[4,4]):
    
    def parsefSSF():
        if type(fSSF)==str:            
            return {'pade':False,'Pade':False,'SSF':True,'ssf':True}[fSSF];
        else:
            return not not fSSF;                  
        
    (dzu,rs,g_r,alpha_r,dT)=fiber_rescale_to(N,fmax=fmax,beta=beta,g=g,alpha=alpha)
    
    solver=MS(N,g=g_r,dt=dzu,alpha=alpha_r,nm=nm,fSSF=parsefSSF());
    
    return (solver,dzu,rs,g_r,alpha_r,dT);
    


    
def js_nls_fiber_rescale_to(N,sname='NLS',fmax=1.0,beta=1.0,g=1.0,alpha=0.0,nm=[4,4]):
    from jsobj import jso
    o=jso();
    (o.solver,o.dzu,o.rs,o.g_r,o.alpha_r,o.dT)=nls_fiber_rescale_to(N=N,sname=sname,fmax=fmax,beta=beta,g=g,alpha=alpha,nm=nm)
    return o; 

def js_manakov_fiber_rescale_to(N,fSSF=False,fmax=1.0,beta=1.0,g=1.0,alpha=0.0,nm=[4,4]):
    from jsobj import jso
    o=jso();
    (o.solver,o.dzu,o.rs,o.g_r,o.alpha_r,o.dT)=manakov_fiber_rescale_to(N=N,fSSF=fSSF,fmax=fmax,beta=beta,g=g,alpha=alpha,nm=nm)
    return o; 





if __name__=='__main__':
    
    psec=10**-12;
    THZ=10**12;
    
    N=2*64*1024;
    
    beta=21*(psec)**2 ;# 21-28 ps**2 /km
    alpha=DBpKM2Alpha(0.2) # 0.025 -0.057
    fmax=100*THZ; # ???
    gamma=1.3 ;#  ~ 1-1.3 1/(W *km)
    alpha_DBpKM=0.025
    
    o=fiber_rescale_to_jso(N,fmax=fmax,beta=beta,alpha_DBpKM=alpha_DBpKM,gamma=gamma)
    solver=manakov_solver_create_from_jso('ssf',o)
    
    #(solver,dzu,rs,g_r,alpha_r,dT)=manakov_fiber_rescale_to(16,'nls');