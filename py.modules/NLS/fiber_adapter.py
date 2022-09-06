# -*- coding: utf-8 -*-
"""
"""
from NLS.nls_pade import *
from jsobj import *    


def fiber_rescale_to(N,dT=None,winT=None,fmax=1.0,beta=1.0,gamma=1.0,alpha=0.0):
    
    N=int(N);
    
    
    if dT is None:
        if winT is None:
            dT=0.5/fmax;
            
        else:
            dT=np.float64(winT)/N;
    else:
        dT=np.float64(dT)
        
    winT=dT*N
    fmax=0.5/dT;
    
    
    half_beta=beta/2.0
    half_alpha=alpha/2.0   
    
    
    rs=half_beta/(dT**2);    
    
    g_r=gamma/rs;
    alpha_r=half_alpha/rs; 
    dzu=rs;        
    return (dzu,rs,g_r,alpha_r,dT,N,fmax);


def DBpKM2Alpha(dB_per_km=0.025):
    return (np.log(10)/10)*dB_per_km;

    
def fiber_rescale_to_jso(N,dT=None,winT=None,fmax=1.0,beta=1.0,gamma=1.0,alpha=None,alpha_DBpKM=0.0):  
    
    
    info=jso(locals());
    
    #beta=beta/2.0
    
    if alpha is None:
        alpha=DBpKM2Alpha(alpha_DBpKM); 
    
        
    
    o=jso();
    (o.dzu,o.rs,o.g_r,o.alpha_r,o.dT,o.N,o.fmax)=fiber_rescale_to(N,dT=dT,winT=winT,fmax=fmax,beta=beta,gamma=gamma,alpha=alpha);
    o.info=info
    return o;


def rescale_dz(dz,dzu):
    sgn=np.sign(dz);
    dz=np.abs(dz);
    
    rep=int(dzu/dz);
    if (dzu-rep*dz) > 1e-11*dz:
        rep+=1;
        
    dz=dzu/rep;   
    
    return (sgn*dz,rep);


def manakov_solver_create_from_jso(o,name='nls',dz_km=1.0,nm=[4,4]):
    
    def parse_name(f):
        if type(f)==str:            
            return {'nls':False,'pade':False,'ssf':True}[f.lower()];
        else:
            return not not f;

    o=to_jso(o)        
    (dz,rep)=rescale_dz(dz_km*o.dzu,o.dzu);
    
    solver=MS(o.N,g=o.g_r,dt=dz,alpha=o.alpha_r,nm=nm,fSSF=parse_name(name));  
    
    return (solver,dz,rep)
    

def manakov_solver_create_from_jso2(o,name='nls',dz_km=1.0,nm=[4,4],falpha=True):    
    
    
    
    def parse_name(f):
        if type(f)==str:            
            return {'nls':False,'pade':False,'ssf':True}[f.lower()];
        else:
            return not not f;

    o=to_jso(o)        
    
    
    
    (dz,rep)=rescale_dz(dz_km*o.dzu,o.dzu);
    
    kwd={
        'pade_factory': lambda f,N,pp : create_nls_context_ex(b"manakov-dec",f,N,pp)
        }
    
    alpha=o.alpha_r if falpha else 0.0; 
    
    solver=MS(o.N,g=o.g_r,dt=dz,alpha=alpha,nm=nm,fSSF=parse_name(name),**kwd);  
    
    return (solver,dz,rep)



def nls_solver_create_from_jso(o,name,dz_km=1.0,nm=[4,4]):
    
    o=to_jso(o)
    (dz,rep)=rescale_dz(dz_km*o.dzu,o.dzu);
    
    solver= NLS_ex(name,N=o.N,g=o.g_r,dt=dz,alpha=o.alpha_r,nm=nm);   
    
    return (solver,dz,rep)


def fiber_rescale_from(N,dT=None,winT=None,fmax=1.0,beta=1.0,gamma=1.0,alpha=None,alpha_DBpKM=0.0):
    
    
    o=fiber_rescale_to_jso(N,dT=dT,winT=winT,fmax=fmax,beta=beta,gamma=gamma,alpha=alpha,alpha_DBpKM=alpha_DBpKM)
    
    o.nls_solver_create= lambda name,dz_km=1.0,nm=[4,4]: nls_solver_create_from_jso(o,name,dz_km=dz_km,nm=nm);
    o.manakov_solver_create= lambda name,dz_km=1.0,nm=[4,4]: manakov_solver_create_from_jso(o,name,dz_km=dz_km,nm=nm);
    o.manakov_solver_create2= lambda name,dz_km=1.0,nm=[4,4]: manakov_solver_create_from_jso2(o,name,dz_km=dz_km,nm=nm);
    
    return o;
    


def pow_norm(x,dT=1.0,energy=1.0):
    x=np.array(x,copy=False);
    Z=np.sqrt(energy/dT)/scipy.linalg.norm(x);
    return Z*x;



if __name__=='__main__':
    
    psec=10**-12;
    THZ=10**12;
    
    N=2*64*1024;
    mr=1e-3
    beta=21*(psec)**2 ;# 21-28 ps**2 /km
    alpha=DBpKM2Alpha(0.2) # 0.025 -0.057
    fmax=100*THZ; # ???
    gamma=1.3 ;#  ~ 1-1.3 1/(W *km)
    alpha_DBpKM=0.025
    dT=1e-12
    o1=fiber_rescale_to_jso(N,fmax=fmax,beta=beta,alpha_DBpKM=alpha_DBpKM,gamma=gamma)
    o=fiber_rescale_from(N,dT=dT,fmax=fmax,beta=beta,alpha_DBpKM=alpha_DBpKM,gamma=gamma)
    
    (solver,dz,rep)=o.nls_solver_create('nls',dz_km=0.001)
    (solver,dz,rep)=o.manakov_solver_create('nls',dz_km=0.001)
    (solver,dz,rep)=o.manakov_solver_create2('nls',dz_km=0.001)
    
    #solver=manakov_solver_create_from_jso('ssf',o)    
    #(solver,dzu,rs,g_r,alpha_r,dT)=manakov_fiber_rescale_to(16,'nls');