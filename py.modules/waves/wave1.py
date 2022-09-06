import numpy as np
from  lipa.trisparse import trimesh2sparse,coo_scheme
from lipa.parallel import LIPA_solver_st,LIPA_solver_pp,lu_solver_factory,solver_test_pp,Tic
from scipy.sparse import coo_matrix,csc_matrix
import jsonrpc.jsonclass as jsncls
import numbers

LIPA_solver=LIPA_solver_st;


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

def make_wave1(Np=50,L=1.0,E=10,rho=10,alpha=0,beta=0,**kws):
    
    
    opts=extend(kws,{'dt':1,'nd':2,'tic_corr':0,'fcycle':False})
    LIPA_solver=cast_str(opts.pop('LIPA_solver',LIPA_solver_pp));
    lu_factory=cast_str(opts.pop('lu_factory',lu_solver_factory()));        
    print(LIPA_solver)    
    
    

    fcycle=opts['fcycle'];
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

    K0=np.array([[1,-1],[-1,1]],dtype='d');
    M0=np.array([[1,0],[0,1]],dtype='d');
    
    
    mKx=np.zeros([Nt,2,2],dtype='d');
    mMx=np.zeros([Nt,2,2],dtype='d');
    mGx=np.zeros([Nt,2,2],dtype='d');

    for k in range(Nt):
        M=rho[k]*L*M0/2;
        K=E[k]*K0/L;
        mKx[k][:][:]=K;
        mMx[k][:][:]=M;
        mGx[k][:][:]=alpha[k]*K+beta[k]*M;
    
    #dt=opts['dt'];
    
    print('Assembly Start ')
    tic=Tic()    
    smM=trimesh2sparse(trs1,Np,mMx,2);
    smK=trimesh2sparse(trs1,Np,mKx,2);
    smG=trimesh2sparse(trs1,Np,mGx,2);
    
    print('Assembly End ',tic.sec())
    
    #    smG=coo_matrix((Np,Np))
    #G0=coo_scheme((Np,Np))
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
    
    
