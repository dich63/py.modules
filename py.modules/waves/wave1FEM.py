import numpy as np
from  lipa.trisparse import trimesh2sparse,coo_scheme
from lipa.parallel import LIPA_solver_st,LIPA_solver_pp,lu_solver_factory,solver_test_pp,Tic
from scipy.sparse import coo_matrix,csc_matrix
import jsonrpc.jsonclass as jsncls
import numbers

LIPA_solver=LIPA_solver_st;


print('import: ',__file__)

def extend(d1,d2=None):
    d={}
    if d2:
        d.update(d2);
    if d1:
        d.update(d1);
    return d


def to_array(a,N):
    if isinstance(a,numbers.Number):
        return a*np.ones(N,dtype='d');    
    else:
        z=np.array(a,dtype='d').flatten();
        return z[:N];


def cast_str(o):
    if type(o) in jsncls.string_types:
        o=eval(o)
    return o;


Np=50;# vertex count



def _make_sf_1D( 
        Lc,
        alpha,
        M,
        eta,
        eta_p,
        kappa,
        rho,
        rho_f,
        rho_m):
    
    
    aM=alpha*M;
    
    mK=[[Lc,aM],
        [aM,M]]
    
    mM=[[rho ,rho_f],
        [rho_f,rho_m]]
    
    mKd=[[eta_p*Lc,0], [0,0]]
    
    mMd=[[0,0], [0,eta/kappa]]
    
    return [ np.array(mm) for mm in [mK,mM,mKd,mMd]]
    
    
def make_sf_1D(         
        m=0.28e0,
        gamma=0.1e1,        
        eta=0.2e11,
        eta_p=0.0e0,
        kappa=0.1e1,
        K_s=0.39e11,
        K_f=0.57e11,
        K_m=0.85e10,
        rho_s=0.265e4,
        rho_f=0.7e3,
        eps=1e-15
        ):
    
    alpha=1.0-K_m/K_s;
    
    KKfs=(K_f+K_s*m)*(K_s-K_m)    
    
    M=K_s**2*K_f/KKfs;
    
    K_c=K_s*(K_m*m*(K_s-K_f)+K_f*(K_s-K_m) )/(K_s*m*(K_s-K_f)+K_f*(K_s-K_m));
    
    Lc=K_c;
    
    rho=m*rho_f+(1-m)*rho_s;
    
    rho_m=gamma*rho_f/m if rho_f>eps else 0.0 ;
    
    
    
    
    return _make_sf_1D(Lc,alpha,M,
        eta,eta_p,kappa,
        rho,
        rho_f,
        rho_m)    
    
    



def make_wave1FEM(Np=50,L=1.0,E=10,rho=10,alpha=0,beta=0,**kws):
    
    
    opts=extend(kws,{'dt':1,'nd':2,'scheme':0,'fcycle':False})
    '''
    LIPA_solver=cast_str(opts.pop('LIPA_solver',LIPA_solver_pp));
    lu_factory=cast_str(opts.pop('lu_factory',lu_solver_factory()));        
    print(LIPA_solver)    
    '''
    

    fcycle=opts['fcycle'];
    scheme=opts['scheme'];
    
    print(opts)
    
    if not fcycle:
        Nt=Np-1;    
        nidx=np.array(range(Nt));
        trs=np.array([nidx,nidx+1],dtype=np.int32);
        trs1=trs.transpose();
    else:
        Nt=Np;
        nidx=np.array(range(Nt));
        nidx2=np.array(range(Nt))+1;
        nidx2[Nt-1]=0;        
        trs=np.array([nidx,nidx2],dtype=np.int32);
        trs1=trs.transpose();
            
    

    (E,rho,alpha,beta)=[ to_array(o,Nt) for o in (E,rho,alpha,beta)]

    if scheme:
        K0=np.array([[1,-1],[-1,1]],dtype='d');
        M0=np.array([[1,1],[1,1]],dtype='d')/2;
    else:
        K0=np.array([[1,-1],[-1,1]],dtype='d');
        M0=np.array([[3,1],[1,3]],dtype='d')/4;
    
    '''
    mKx=np.zeros([Nt,2,2],dtype='d');
    mMx=np.zeros([Nt,2,2],dtype='d');
    mGx=np.zeros([Nt,2,2],dtype='d');
    '''
    mMx,mKx,mGx=[np.zeros([Nt,2,2],dtype='d') for k in range(3) ]

    for k in range(Nt):
        M=rho[k]*L*M0/2;
        K=E[k]*K0/L;
        mKx[k][:][:]=K;
        mMx[k][:][:]=M;
        mGx[k][:][:]=alpha[k]*K+beta[k]*M;
    
    #dt=opts['dt'];
    
    print('Assembly Start ')
    tic=Tic()    
    
    '''
    smM=trimesh2sparse(trs1,Np,mMx,2);
    smK=trimesh2sparse(trs1,Np,mKx,2);
    smG=trimesh2sparse(trs1,Np,mGx,2);
    '''
    
    smM,smK,smG =[ trimesh2sparse(trs1,Np,mx,2).tocoo()  for mx in [mMx,mKx,mGx] ]
    
    print('Assembly End ',tic.sec())
    
    return smM,smK,smG
    
    #    smG=coo_matrix((Np,Np))
    #G0=coo_scheme((Np,Np))
    
    
'''    
    tic=Tic()    

    
    
    solver=LIPA_solver((smK,smG,smM),SolverFactory=lu_factory,**opts);     
    t=tic.sec();
    solver.tLU=t;
    print('tLu=',t);

    def energy(fsum=1):
        xx=solver.xn;
        x0,x1=(xx[0],xx[1]);
        dx0=(x0[trs[0]]+x0[trs[1]])/2;
        dx1=(x1[trs[0]]-x1[trs[1]])/2;
        xx0,xx1=(dx0*dx0,dx1*dx1);
        e=rho*xx0+E*xx1;
        #e=rho*xx1+E*xx0;
        return  np.sum(e)  if fsum else e;

    def momentum(fsum=1):
        xx=solver.xn;
        x1=xx[1];
        x1=(x1[trs[0]]-x1[trs[1]])/2;        
        p=rho*x1;
        #e=rho*xx1+E*xx0;
        return  np.sum(p)  if fsum else p;

    solver.energy=energy;
    solver.momentum=momentum;
    
    return solver
    
'''    

if __name__=='__main__':
    from utils import *
    from LIPA2.qp_solver import *
    from klu_lipa import *
    
    norm=np.linalg.norm
    normm= lambda x:  norm(x,ord=np.inf)
    
    N=1024*16;
    x0=np.random.rand(N)
    dt=0.1;
    
    [smM,smK,smG]=make_wave1FEM(N)
    
    pade_nm=[8,8]
    nd=2;
    lb=lipa_qp_base_t([smK,smG,smM],pade_nm,'sp_LU_factory',nd=nd);
    tic()
    lb.reset(dt);
    toc('lu:')
    lb.x=x0;
    
    nT=1000
    
    tic()
    lb(nT)
    toc('tic')
    
    x=lb.x
    
    tic()
    solver=LIPA_solver_ctx([smK,smG,smM],
                           pade_nm=pade_nm,dt=dt,
                           fcomplex=1,
                           parallel=1,
                           asyn=0,
                           nd=nd)
    
    toc('lu-ctx:')
    
    solver.x=x0;
    tic()
    solver.step(nT)
    toc('tic')
    
    y=solver.x
    err=normm(y-x)/np.mean(np.abs(x))
    printf("err=%f%%",err)

