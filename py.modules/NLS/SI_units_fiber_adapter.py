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



def DBpMeter2Alpha(dB_per_m):
    return (np.log(10)/10)*dB_per_m;


def DBpKM2AlphaSI(dB_per_Km=0.2171472409516259):
    return DBpMeter2Alpha(dB_per_Km/1000.0);


    
def fiber_rescale_to_jso(N,dT=None,winT=None,fmax=1.0,beta=1.0,gamma=1.0,alpha=None,alpha_DBpKM=0.0):  
    
    
    info=jso(locals());
    
    #beta=beta/2.0
    
    if alpha is None:
        alpha=DBpKM2AlphaSI(alpha_DBpKM); 
    
        
    
    o=jso();
    (o.dzu,o.rs,o.g_r,o.alpha_r,o.dT,o.N,o.fmax)=fiber_rescale_to(N,dT=dT,winT=winT,fmax=fmax,beta=beta,gamma=gamma,alpha=alpha);
    o.info=info
    return o;


def rescale_dz(dz,dzu):
    sgn=np.sign(dz);
    dz=np.abs(dz);
    
    
    
    rep=int(dz/dzu);
    
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
    

def manakov_solver_create_from_jso2(o,name='nls',dz_m=None,dz_km=1.0,nm=[4,4],falpha=True):    
    
    
    if dz_m is None:
        dz_m=dz_km*1000.0;
    
    
    def parse_name(f):
        if type(f)==str:            
            return {'nls':False,'pade':False,'ssf':True}[f.lower()];
        else:
            return not not f;
        
    def rescale_rep(dz,dzu):   
        sgn,dz=np.sign(dz),np.abs(dz);              
        
        rep=int(dz/dzu);
        
        if (dz-rep*dzu) > 1e-11*dzu:
            rep+=1;
            
        dzu=dz/rep; 
             
        return (dzu,rep);
        
        

    o=to_jso(o)        
    
    
    dzu=o.dzu;
    
    dz_m*=dzu;
    
    (dz_r,rep)=rescale_rep(dz_m,dzu);
    
    
    
    kwd={
        'pade_factory': lambda f,N,pp : create_nls_context_ex(b"manakov-dec",f,N,pp)
        }
    
    alpha=o.alpha_r if falpha else 0.0; 
    
    solver=MS(o.N,g=o.g_r,dt=dz_m,alpha=alpha,nm=nm,fSSF=parse_name(name),**kwd);  
    
    return (solver,dz_m,rep)



def nls_solver_create_from_jso(o,name,dz_km=1.0,nm=[4,4]):
    
    o=to_jso(o)
    (dz,rep)=rescale_dz(dz_km*o.dzu,o.dzu);
    
    solver= NLS_ex(name,N=o.N,g=o.g_r,dt=dz,alpha=o.alpha_r,nm=nm);   
    
    return (solver,dz,rep)


def fiber_rescale_from(N,dT=None,winT=None,fmax=1.0,beta=1.0,gamma=1.0,alpha=None,alpha_DBpKM=0.0):
    
    
    o=fiber_rescale_to_jso(N,dT=dT,winT=winT,fmax=fmax,beta=beta,gamma=gamma,alpha=alpha,alpha_DBpKM=alpha_DBpKM)
    
    o.nls_solver_create= lambda name,dz_km=1.0,nm=[4,4]: nls_solver_create_from_jso(o,name,dz_km=dz_km,nm=nm);
    o.manakov_solver_create= lambda name,dz_km=1.0,nm=[4,4]: manakov_solver_create_from_jso(o,name,dz_km=dz_km,nm=nm);
    o.manakov_solver_create2= lambda name,dz_km=1.0,dz_m=None,nm=[4,4]: manakov_solver_create_from_jso2(o,name,dz_km=dz_km,dz_m=dz_m,nm=nm);
    
    return o;
    


def pow_norm(x,dT=1.0,energy=1.0):
    x=np.array(x,copy=False);
    Z=np.sqrt(energy/dT)/scipy.linalg.norm(x);
    return Z*x;



if __name__=='__main__':
    
    psec=10**-12;
    THZ=10**12;
    
    N=128*1024;
    kM=1000;
    
    mr=1e-3
    beta=20*(psec)**2 /kM  ;# 21-28 ps**2 /km
    alpha=0.05 /kM   # 0.05 -0.057    
    gamma=1.3 /kM ;#  ~ 1-1.3 1/(W *km)
    
    dT=1e-12
    
    
    o1=fiber_rescale_to_jso(N,dT=dT,beta=beta,alpha=alpha,gamma=gamma)
    o=fiber_rescale_from(N,dT=dT,beta=beta,alpha=alpha,gamma=gamma)
    
    manakov_solver_create2=o.manakov_solver_create2;
    (solver,dz,rep)=manakov_solver_create2('nls',dz_km=11)
    print('[dz,rep]=',[dz,rep])
    
    
    #(solver,dz,rep)=o.nls_solver_create('nls',dz_km=0.001)
    #(solver,dz,rep)=o.manakov_solver_create('nls',dz_km=0.001)
    
    
    #solver=manakov_solver_create_from_jso('ssf',o)    
    #(solver,dzu,rs,g_r,alpha_r,dT)=manakov_fiber_rescale_to(16,'nls');