import numpy as np
from  lipa.trisparse import trimesh2sparse,coo_scheme
from lipa.parallel import LIPA_solver_st,LIPA_solver_pp,lu_solver_factory,solver_test_pp,Tic

from lipa.kernel import calc_Az,convert_tocsc,bicgstab_solver_factory
from scipy.sparse import coo_matrix,csc_matrix
import jsonrpc.jsonclass as jsncls
import numbers
import types
import copy
import sys

from klu_lipa import *

LIPA_solver=LIPA_solver_st;
norm=np.linalg.norm;

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
    
    
    opts=extend(kws,{'dt':1,'nd':3,'tic_corr':0,'fcycle':False,'scheme':0,'fCp':0})
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

    
    scheme=opts['scheme'];
    print('scheme=',scheme)
    if scheme==1:    
        K0=np.array([[1,-1],[-1,1]],dtype='d');
        M0=np.array([[3./4,1./4],[1./4,3./4]],dtype='d');
        K00=np.array([[1,0],[-1,1]],dtype='d');
        M00=np.array([[0,0],[1./4,3./4]],dtype='d');
    elif scheme==0:
        K0=np.array([[1,-1],[-1,1]],dtype='d');
        M0=np.array([[1,0],[0,1]],dtype='d');
        K00=np.array([[1,0],[-1,1]],dtype='d');
        M00=np.array([[0,0],[0,1]],dtype='d');
    else: 
        raise Exception('unknown scheme')
        

    G00=np.array([[1e-5,0],[0,0]],dtype='d');

    mKx=np.zeros([Nt,2,2],dtype='d');
    mMx=np.zeros([Nt,2,2],dtype='d');
    mGx=np.zeros([Nt,2,2],dtype='d');
    

    mDx=np.zeros([Nt,2,2],dtype='d');

    mMx3=np.zeros([Nt,2,2],dtype='d');

    D=K0/L;
    ML2=(L/2.)*M0;
    G0=G00/L;
    D0=K00/L;
    ML20=(L/2.)*M00;

    for k in range(Nt):
        M=rho[k]*ML2;        
        K=E[k]*D;
        mDx[k][:][:]=D;
        mKx[k][:][:]=K;
        mMx[k][:][:]=M;
        mGx[k][:][:]=alpha[k]*K+beta[k]*M;
    
    #dt=opts['dt'];
     
    fbcr=opts.get('fbcr',0);
    if fbcr:
        M0=rho[0]*ML20;        
        K0=E[0]*D0;
        mKx[0][:][:]=K0;
        mMx[0][:][:]=M0;
        mGx[0][:][:]=G0;



    print('Assembly Start ')
    tic=Tic()    
    smM=trimesh2sparse(trs1,Np,mMx,2);
    smK=trimesh2sparse(trs1,Np,mKx,2);
    smG=trimesh2sparse(trs1,Np,mGx,2);

    smD=trimesh2sparse(trs1,Np,mDx,2);

    smM3=trimesh2sparse(trs1,Np,mMx3,2);
    
    print('Assembly End ',tic.sec())
    
    #    smG=coo_matrix((Np,Np))
    #G0=coo_scheme((Np,Np))
    tic=Tic()    

    
    
    solver=LIPA_solver((smK,smG,smM),SolverFactory=lu_factory,**opts);     
    
            
    t=tic.sec();
    solver.tLU=t;
    print('tLu=',t);

    sM=smM.tocsc();
    sK=smK.tocsc();
    sD=smD.tocsc();

    def energy(fsum=1):
        xx=solver.xn;
        x0,x1=(xx[0],xx[1]);
        #dx0=(x0[trs[0]]+x0[trs[1]])/2;
        #dx1=(x1[trs[0]]-x1[trs[1]])/2;
        #xx0,xx1=(dx0*dx0,dx1*dx1);
        #e=rho*xx0+E*xx1;
        k=sM*x1;
        v=sK*x0;
        dx=sD*x0;
        ek=x1*k
        ev=dx*v;        
        e=ek+ev;
        #e=rho*xx1+E*xx0;
        return  np.sum(e)  if fsum else e;

    def momentum(fsum=1):
        xx=solver.xn;
        x1=xx[1];
        p=sM*x1;
        #x1=(x1[trs[0]]-x1[trs[1]])/2;        
        #p=rho*x1;
        #e=rho*xx1+E*xx0;
        return  np.sum(p)  if fsum else p;

    solver.energy=energy;
    solver.momentum=momentum;
    return solver


def step_delta(self):
    pass
    xn0=xx_save=np.copy(self.xn);

    self.vxn[:]=self.dxn;
    nd=xn0.shape[0];
    J=self.zeros();
    for k in range(nd):
        J-=self.dAC[k]*xn0[k];
    self.reset_J(J);
    self.step();
    self.dxn[:]=self.vxn;
    self.vxn[:]=xx_save;

    self.reset_J(0);

    return self.dxn;
    
def step_and_delta(self,repeat=1):
    for r in range(repeat):
        self.step();
        self.step_delta();
    return (self.vxn,self.dxn);

    
def make_wave1delta(Np=50,L=1.0,E=10,rho=10,dE=0,drho=0,alpha=0,beta=0,**kws):
    
    
    opts=extend(kws,{'dt':1,'nd':3,'tic_corr':0,'fcycle':False,'scheme':0,'fCp':0,'fbcr':0,'fbcl':0})
    LIPA_solver=LIPA_solver_st #cast_str(opts.pop('LIPA_solver',LIPA_solver_pp));
    LIPA_solver=cast_str(opts.pop('LIPA_solver',LIPA_solver_pp));
    lu_factory=cast_str(opts.pop('SolverFactory',lu_solver_factory));        
    lu_factory=lu_factory()
    opts['SolverFactory']=lu_factory;
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
            
    

    (E,rho,alpha,beta,dE,drho)=[ to_array(o,Nt) for o in (E,rho,alpha,beta,dE,drho)]

    scheme=opts['scheme'];
    print('scheme=',scheme)
    if scheme==1:    
        K0=np.array([[1,-1],[-1,1]],dtype='d');
        M0=np.array([[3./4,1./4],[1./4,3./4]],dtype='d');
    elif scheme==0:
        K0=np.array([[1,-1],[-1,1]],dtype='d');
        M0=np.array([[1,0],[0,1]],dtype='d');
    else: 
        raise Exception('unknown scheme')
        
    M0R=np.array([[0,0],[0,1]],dtype='d');
    K0R=np.array([[1,-1],[0,1]],dtype='d');
    G0R=np.array([[0,0],[0,1]],dtype='d');


    M0L=np.array([[1,0],[0,0]],dtype='d');
    K0L=np.array([[1,0],[1,-1]],dtype='d');
    G0L=np.array([[1,0],[0,0]],dtype='d');

        
    mKx=np.zeros([Nt,2,2],dtype='d');
    mdKx=np.zeros([Nt,2,2],dtype='d');
    mMx=np.zeros([Nt,2,2],dtype='d');
    mdMx=np.zeros([Nt,2,2],dtype='d');
    mGx=np.zeros([Nt,2,2],dtype='d');
    mdGx=np.zeros([Nt,2,2],dtype='d');
    mDx=np.zeros([Nt,2,2],dtype='d');

    D=K0/L;
    ML2=(L/2.)*M0;

    for k in range(Nt):
        M=rho[k]*ML2;        
        K=E[k]*D;
        dM=drho[k]*ML2;        
        dK=dE[k]*D;
        mDx[k][:][:]=D;
        mKx[k][:][:]=K;
        mMx[k][:][:]=M;
        mdKx[k][:][:]=dK;
        mdMx[k][:][:]=dM;

        mGx[k][:][:]=alpha[k]*K+beta[k]*M;
    
    #dt=opts['dt'];
    
    if opts['fbcr']:
        k=Nt-1;
        M=M0R;        
        K=1*E[k]*K0R/L;
        G=G0R;
        mKx[k][:][:]=K;
        mMx[k][:][:]=M;
        mGx[k][:][:]=G0R;

    if opts['fbcl']:
        k=0;
        M=M0L;        
        K=1*E[k]*K0L/L;
        G=G0L;
        mKx[k][:][:]=K;
        mMx[k][:][:]=M;
        mGx[k][:][:]=G0L;
        
        
    print('Assembly Start ')
    tic=Tic()    
    smM=trimesh2sparse(trs1,Np,mMx,2);
    smK=trimesh2sparse(trs1,Np,mKx,2);
    smG=trimesh2sparse(trs1,Np,mGx,2);
    smD=trimesh2sparse(trs1,Np,mDx,2);


    smdM=trimesh2sparse(trs1,Np,mdMx,2);
    smdK=trimesh2sparse(trs1,Np,mdKx,2);
    smdG=trimesh2sparse(trs1,Np,mdGx,2);
    
    print('Assembly End ',tic.sec())
    
    #    smG=coo_matrix((Np,Np))
    #G0=coo_scheme((Np,Np))
    tic=Tic()    

    
    (smK,smG,smM)=convert_tocsc((smK,smG,smM))
    #solver=LIPA_solver((smK,smG,smM),SolverFactory=lu_factory,**opts);     
    #
    solver=LIPA_solver((smK,smG,smM),**opts);     
    #    solver=LIPA_solver_ctx((smK,smG,smM),**opts);     
    t=tic.sec();
    solver.reset_corr(0*1);
    solver.dxn=np.zeros_like(solver.xn);
    
    solver.dAC=convert_tocsc((smdK,smdG,smdM));
    if sys.version_info.major==3:
        solver.step_delta=types.MethodType(step_delta,solver);
        solver.step_and_delta=types.MethodType(step_and_delta,solver);
    else:
        solver.step_delta=types.MethodType(step_delta,solver,LIPA_solver);
        solver.step_and_delta=types.MethodType(step_and_delta,solver,LIPA_solver);

    
    #js.fu=types.MethodType(fu, js, jsncls.jsobject) 
    
    solver.tLU=t;
    print('tLu=',t);

    sM=smM.tocsc();
    sK=smK.tocsc();
    sD=smD.tocsc();

    def energy(fsum=1):
        xx=solver.xn;
        x0,x1=(xx[0],xx[1]);
        #dx0=(x0[trs[0]]+x0[trs[1]])/2;
        #dx1=(x1[trs[0]]-x1[trs[1]])/2;
        #xx0,xx1=(dx0*dx0,dx1*dx1);
        #e=rho*xx0+E*xx1;
        k=sM*x1;
        v=sK*x0;
        dx=sD*x0;
        ek=x1*k
        ev=dx*v;        
        e=ek+ev;
        #e=rho*xx1+E*xx0;
        return  np.sum(e)  if fsum else e;

    def momentum(fsum=1):
        xx=solver.xn;
        x1=xx[1];
        p=sM*x1;
        #x1=(x1[trs[0]]-x1[trs[1]])/2;        
        #p=rho*x1;
        #e=rho*xx1+E*xx0;
        return  np.sum(p)  if fsum else p;

    solver.energy=energy;
    solver.momentum=momentum;
    return solver
    
