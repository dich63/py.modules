import numpy as np
from klu_lipa import *
from scipy.sparse import spdiags,dia_matrix

def init_nls_eq(r2,g=1.0):#,lipa_solver=LIPA_solver_ctx):
    
    N=np.size(r2);
    o=np.ones(N,dtype = np.complex);
    on=g*r2-1.0;
    #on[:]=range(1,N+1)
    oz=np.zeros(N+2,dtype = np.complex);
    oz[1:-1]=on;
    #ou,od=on[0],on[-1];
    #H=spdiags([on,2*o,oz[:-1]],[-1,0,1],N,N);
    H=spdiags((oz[2:],2*o,oz[:-2]),[-1,0,1],N,N);
    H=H.tocsc();

    #H[0,-1]=od;
    #H[-1,0]=od;
    M=spdiags([-1j*o],[0],N,N);
    return (H,M.tocsc());


def init_nls_eqU(r2,g=1.0,dec=0):#,lipa_solver=LIPA_solver_ctx):
    
    N=np.size(r2);
    N1=N-1;
    o=np.ones(N,dtype = np.complex);
    on=-g*r2-2.0-dec;
    

    H=dia_matrix((N,N),dtype = np.complex);
    H.setdiag(on,0);
    H.setdiag(o,1);
    H.setdiag(o,-1);
    H.setdiag(o[0],N1);
    H.setdiag(o[0],-N1);


    #on[:]=range(1,N+1)
    
    #ou,od=on[0],on[-1];
    #H=spdiags([o,on,o],[-1,0,1],N,N);


    
    

    #H[0,-1]=1;
    #H[-1,0]=1;
    M=spdiags([1j*o],[0],N,N);
    H=H.tocsc();
    M=M.tocsc();
    return (H,M);


def nls(dt,ff,pade_nm=(4,4),g=1.0,rep=1,dec=0.0,LIPA_solver=LIPA_solver_ctx):
    
    ff=ff.astype(np.complex,copy=False);
    if not ff.shape[0]==1:
      ff=ff.reshape((1,ff.size));

    #print(dt,pade_nm,g)
    #print(ff)
    dec=1j*dec;
    rep=int(rep)
    while rep>0:
        r2=ff[0]*ff[0].conj();
        (H,M)=init_nls_eqU(r2,g,dec);
        solver=LIPA_solver((H,M),dt=dt,pade_nm=pade_nm,xx=ff,fcomplex=True);
        solver.xn=ff;
        solver.step();
        ff=solver.xn;
        rep-=1;
    return ff;
    pass

def nls_e(dt,ff,pade_nm=(4,4),g=1.0,rep=1,LIPA_solver=LIPA_solver_ctx):
    
    ff=ff.astype(np.complex,copy=False);
    if not ff.shape[0]==1:
      ff=ff.reshape((1,ff.size));

    #print(dt,pade_nm,g)
    #print(ff)
    rep=int(rep)
    while rep>0:
        r2=ff[0]*ff[0].conj();
        r2l=np.roll(r2,1);
        r2r=np.roll(r2,-1);
        r2=(r2+r2l+r2r)/3.0;
        (H,M)=init_nls_eqU(r2,g);
        solver=LIPA_solver((H,M),dt=dt,pade_nm=pade_nm,xx=ff,fcomplex=True);
        solver.xn=ff;
        solver.step();
        ff=solver.xn;
        rep-=1;
    return ff;
    pass


if __name__=='__main__':
    from utils import *
    norm=np.linalg.norm;

    #(H,M)=init_nls_eq(np.array([0,0,777.0,0,0]))
    N=5;
    f=np.zeros((1,N),dtype = np.complex );
    f[0,(N//2):((3*N)//4)]=11;
    #f[0,(N//2):(N//2+80)]=11;
    dt=4;
    f=(10+0j)*np.random.randn(1,N)
    tic()
    g=0
    print('norm(0):',norm(f)**2)
    ff=nls(dt,f,g=g);
    print('norm(t):',norm(ff-f)**2)
    ff=nls(dt,ff,g=g);
    print('norm(2t):',norm(ff)**2)
    ff=nls(dt,ff,g=g);
    print('norm(3t):',norm(ff)**2)
    ff=nls(dt,ff,g=g);
    print('norm(4t):',norm(ff)**2)
    t=toc();
    print('t sec:',t)
    pass


