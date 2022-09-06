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


def fiber_rescale_to(N,fmax=1.0,beta=1.0,g=1.0,alpha=0.0):    
    rs=beta*(fmax**2);
    g_r=g/rs;
    alpha_r=alpha/rs;
    dT=1.0/fmax;
    dzu=1.0*rs;        
    return (dzu,rs,g_r,alpha_r,dT);
    
def js_fiber_rescale_to(N,fmax=1.0,beta=1.0,g=1.0,alpha=0.0):
    from jsobj import jso
    o=jso();
    (o.dzu,o.rs,o.g_r,o.alpha_r,o.dT)=fiber_rescale_to(N,fmax=fmax,beta=beta,g=g,alpha=alpha);
    return o;

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
            return {'NLS':False,'nls':False,'SSF':True,'ssf':True}[fSSF];
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



def DBpKM2Alpha(dB_per_km=0.2):
    return (np.log(10)/20)*dB_per_km;


if __name__=='__main__':
    
    psec=10**-12;
    THZ=10**12;
    
    N=2*64*1024;
    
    beta=21*(psec)**2 ;# 21-28 ps**2 /km
    alpha=DBpKM2Alpha(0.2) # 0.025 -0.057
    fmax=100*THZ; # ???
    gamma=1.3 ;#  ~ 1-1.3 1/(W *km)
    o=js_fiber_rescale_to(N,fmax=fmax,beta=beta,alpha=alpha,g=gamma)
    
    (solver,dzu,rs,g_r,alpha_r,dT)=manakov_fiber_rescale_to(16,'nls');