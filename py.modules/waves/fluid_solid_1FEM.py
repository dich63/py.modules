# -*- coding: utf-8 -*-
"""
Created on Tue Jan 12 07:41:53 2016

@author: dich6_000
"""


import numpy as np
from  lipa.trisparse import trimesh2sparse,coo_scheme
#from lipa.parallel import LIPA_solver_st,LIPA_solver_pp,lu_solver_factory,solver_test_pp,Tic
from scipy.sparse import coo_matrix,csc_matrix
import jsonrpc.jsonclass as jsncls
import numbers




print('import: ',__file__)

def extend(d1,d2=None):
    
    d={}
    if d2:
        d.update(d2);
    if d1:
        d.update(d1);
    return d

def ext_def(opts,**dflt):    
    t={}
    t.update(dflt)   
    t.update(opts)
     
    return t;

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

def _get_locals(l):
    d={};
    d.update(l)
    return d;

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


    
def default_params_sf_1D():
    return dict(         
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
        );

def make_params_sf_1D(**opts):
    return ext_def(default_params_sf_1D(),**opts)

def make_params_list_sf_1D(*lst):    
    o=default_params_sf_1D();
    for l in lst:
        o=ext_def(o,l)
    return o;
    
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
    
    if np.abs(m)>eps:
        K_c=K_s*(K_m*m*(K_s-K_f)+K_f*(K_s-K_m) )/(K_s*m*(K_s-K_f)+K_f*(K_s-K_m));
    else:
        K_c=K_s;
        
    
    Lc=K_c;
    
    
    if m<eps:
        rho_f=0.0;
        
    
    rho=m*rho_f+(1-m)*rho_s;
    
    
    
    rho_m=gamma*rho_f/m if rho_f>eps  else 0.0 ;
    
    
    
    
    return _make_sf_1D(Lc,alpha,M,
        eta,eta_p,kappa,
        rho,
        rho_f,
        rho_m)    
    

def make_sf_FEM_data_old(Nt,FSM,dx=1.0,**opts):
    
    
    def _to_arr(a,N):
        a=np.array(a,copy=False);
        A=np.empty((N,2,2),dtype=a.dtype)
        A[:]=a;
        return A;
        
        
    mK,mM,mKd,mMd =[ _to_arr(mm,Nt) for mm in FSM]
    
    opts=ext_def(opts,
                 K1=np.array([[1,-1],[-1,1]],dtype=np.float64)/dx,
                 M1=np.array([[3,1],[1,3]],dtype=np.float64)*dx/4
                 )
    
    K1,M1=opts['K1'],opts['M1']
    
    #oneT=np.ones([Np,1],dtype=np.float64)
    
    mMx,mKx,mGx=[np.empty((Nt,4,4),dtype=np.float64) for k in range(3)]
    
    
    
    for n in range(Nt):
        mKx[n]=np.kron(K1,mK[n]);
        mMx[n]=np.kron(M1,mM[n]);
        mGx[n]=np.kron(K1,mKd[n])+np.kron(M1,mMd[n]);       
        
    # set breakpoint
    #import pdb; pdb.set_trace()
    return (mKx,mMx,mGx)

def make_sf_FEM_data(Nt,FSM,dx=1.0,**opts):
    
    
    def _to_arr(a,N):
        a=np.array(a,copy=False);
        A=np.empty((N,2,2),dtype=a.dtype)
        A[:]=a;
        return A;
        
        
    mK,mM,mKd,mMd =[ _to_arr(mm,Nt) for mm in FSM]
    
    opts=ext_def(opts,
                 K1=np.array([[1,-1],[-1,1]],dtype=np.float64)/dx,
                 M1=np.array([[3,1],[1,3]],dtype=np.float64)*dx/4
                 )
    
    K1,M1=opts['K1'],opts['M1']
    
    #oneT=np.ones([Np,1],dtype=np.float64)
    
    mKx=np.kron(K1,mK);
    mMx=np.kron(M1,mM);
    mGx=np.kron(K1,mKd)+np.kron(M1,mMd);       

    
    
    '''
    mMx,mKx,mGx=[np.empty((Nt,4,4),dtype=np.float64) for k in range(3)]
    
    
    
    for n in range(Nt):
        mKx[n]=np.kron(K1,mK[n]);
        mMx[n]=np.kron(M1,mM[n]);
        mGx[n]=np.kron(K1,mKd[n])+np.kron(M1,mMd[n]);       
    '''    
    # set breakpoint
    #import pdb; pdb.set_trace()
    
    return (mKx,mMx,mGx)


def make_trs_1D(Nt,C=1,fcycle=False,dtype=np.int32):
    
    if fcycle:
        Np=Nt;
        Nt-=1;
        trs=np.empty((Np,2*C),dtype=dtype);        
        #trs[Nt]=[2*Nt,2*Nt+1,0,1];
        ac=np.arange(C)        
        trs[Nt][C:]=ac
        trs[Nt][:C]= C*Nt+ ac
        
    else:
        Np=Nt+1;
        trs=np.empty((Nt,2*C),dtype=dtype);
        
    #t0=np.array([0,1,2,3],dtype=dtype)    
    t0=np.arange(2*C,dtype=dtype)    
    
    for n in range(Nt):
        trs[n]=t0;
        t0+=C; 
        
    return trs;

from FEM.trs  import *
from FEM.FEM2sparse  import *
#make_trs_1D=FEM.trs.make_trs_1D    
    
def make_sf_FEM(FSM,fcycle=False,fmt='coo',fbatch=False):  
    
    fl= type(FSM) in (tuple,list);
    
    if not fl:
        FSM=(FSM,)
        
    Nt=FSM[0].shape[0];
    '''
    if fcycle:
        Np=Nt;
        Nt-=1;
        trs=np.empty((Np,4),dtype=np.int32);
        trs[Nt]=[2*Nt,2*Nt+1,0,1];        
    else:
        Np=Nt+1;
        trs=np.empty((Nt,4),dtype=np.int32);
        
    t0=np.array([0,1,2,3],dtype=np.int32)    
    
    for n in range(Nt):
        trs[n]=t0;
        t0+=2;        
        
    FSM =[ trimesh2sparse(trs,2*Np,mx,4)  for mx in FSM]
    
    if fmt=='coo':
        FSM=[sm.tocoo() for sm in FSM ]
    elif fmt=='csr':
        FSM=[sm.tocsr() for sm in FSM ]
    elif fmt=='csc':
        FSM=[sm.tocsc() for sm in FSM ]
    '''
    
    Np = Nt if fcycle else Nt+1
    
    trs=make_trs_1D(Np,C=2,fcycle=fcycle,dtype=np.int32)    
    
    FSM =FEM2data_coo(trs,FSM)
    #FSM =[ trimesh2sparse(trs,2*Np,mx,4)  for mx in FSM]       
    if fbatch:
        return FSM;
    
    if fmt=='coo':
        FSM=[sm for sm in FSM ]
    elif fmt=='csr':
        FSM=[sm for sm in FSM.tocsr() ]
    elif fmt=='csc':
        FSM=[sm.tocsc() for sm in FSM ]
    
    
    return FSM if fl else FSM[0];


def make_sf_1D_rgn(
        Np,
        dx=1,        
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
        eps=1e-15,
        **kwd
        ):
    
    
    (mKx,mMx,mKdx,mMdx)=make_sf_1D(
        m,gamma,     
        eta,eta_p,
        kappa,
        K_s,K_f,K_m,
        rho_s,rho_f,eps);
    
    
    params=_get_locals(locals());
    
    
    
    mK,mM,mG=make_sf_FEM_data(Np,(mKx,mMx,mKdx,mMdx),dx)
    
    return mK,mM,mG,params
    

def make_sf_1D_full(
        Np,
        dx=1,
        fcycle=0,
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
        eps=1e-15,
        **kwd
        ):
    
    N=Np if fcycle else Np-1
    
    
    
    mK,mM,mG,params=make_sf_1D_rgn(
            N,dx,
            m,gamma,     
            eta,eta_p,
            kappa,
            K_s,K_f,K_m,
            rho_s,rho_f,eps);
    
    '''
    (mKx,mMx,mKdx,mMdx)=make_sf_1D(
        m,gamma,     
        eta,eta_p,
        kappa,
        K_s,K_f,K_m,
        rho_s,rho_f,eps);
    
    
    params=_get_locals(locals());
    
    N=Np if fcycle else Np-1
    
    mK,mM,mG=make_sf_FEM_data(N,(mKx,mMx,mKdx,mMdx),dx)
    #'''
    
    smK,smG,smM=make_sf_FEM([mK,mG,mM],fcycle=fcycle)
    
    
    return smK,smG,smM,params;




if __name__=='__main__':
    
    
    from lipa.parallel import LIPA_solver_st,LIPA_solver_pp,lu_solver_factory,solver_test_pp,Tic
    from utils import *
    from LIPA2.qp_solver import *
    from klu_lipa import *
    
    LIPA_solver=LIPA_solver_st;
    
    norm=np.linalg.norm
    normm= lambda x:  norm(np.array(x,copy=False).reshape(-1),ord=np.inf)
    #normm= lambda x:  norm(x,ord=np.inf)
    
    t=make_trs_1D(4,C=1,fcycle=1)
    
    N=5;
    ff=np.empty((N,4,4));
    ff[:]=np.kron([[1,-1],[-1,1]],[[1,0.1],[0.1,1]])
    #ff[:]=np.kron([[1,-1],[-1,1]],[[1,1],[1,1]])
    
    
    mK,mM,mG,params=make_sf_1D_rgn(2)
    
    sff=make_sf_FEM(ff,fcycle=1)
    
    print(sff.todense(),sff.shape)
    #raise SystemExit(0);
    
    N=1024*8;
    N=3;
    N=4*1024
    N=128*1024
    N=64*1024
    #N=2**18
    #N=3
    #N=4*1024*1024
    #x0=np.random.rand(N)
    dt=0.1;
    
    make_sf_1D_full(2)
    
    FSM=make_sf_1D(rho_f=0);
    print('start...')
    tic()    
    mK,mM,mG=make_sf_FEM_data_old(N,FSM)
    toc('assemble old:')
    tic()
    mKo,mMo,mGo=make_sf_FEM_data(N,FSM)
    toc('assemble new:')
        
    tic()
    [smK,smG,smM]=make_sf_FEM([mK,mG,mM])
    
    
    
    if 1:
        import jsonrpc.json_io as jio
        fn=sprintf("o:/__ss/matrix/sfFEM%dk.json",int(N/1024));         
        jio.encode({'M':smM,'K':smK,'G':smG},fn)
    
    
    
    toc('sparse:')
    
    lb=lipa_qp_base_t([smK,smG,smM])
    
    
    raise SystemExit(0);
    
    