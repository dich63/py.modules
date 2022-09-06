#
from FEM.FEM2sparse import *
from FEM.trs import *

from jet.jet_csr_iklu import *


class FEM_jet_t(object):
    def __init__(this,trs,N=None,maskDiff=None,Dmax=None,lbound=0,dtype=complex):
        
        trs=np.array(trs,copy=False);        
        
        N=np.max(trs)+1-lbound if N is None else int(N)
        
        this.N=N;          
        this.maskDiff=maskDiff;          
        this.Dmax=Dmax;          
        this.dtype=dtype;          
        
        this.jetsolver=None;
        smbcoo=this.smbcoo=FEM2data_coo(trs,lbound=lbound,shape=(N,N));        
        smbcoo.set_csr_scheme();        
        
        pass
    
    def set_FEM_data(this,datas,jetsolver=jet_csr_klu_t,**opts):
        
        spmb=this.smbcoo.tocsr(datas=datas,dtype=this.dtype)
        
        jetsolver=this.jetsolver
        
        if jetsolver is None:                   
            jetsolver=this.jetsolver=jet_csr_klu_t(spmb,
                                         maskDiff=this.maskDiff,
                                         Dmax=this.Dmax,**opts)
        else:
            jetsolver.spmb=spmb;
        
        pass
    
    def preset(this,dt=1,LM=(4,4)):
        this.jetsolver.preset(dt=dt,LM=LM)
        return this
    def factorize(this): 
        this.jetsolver.factorize();
        return this;
    
    def reset(this,dt=1,LM=(4,4)):
        this.jetsolver.reset(dt=dt,LM=LM)
        return this
    
    def reset_sources(this,source=(None,None),singular_source=(None,None)):        
        def _p(s):
            return (None,None) if s is None else s;
                              
        this.jetsolver.reset_source(*_p(source)).reset_singular_source(*_p(singular_source));        
        return this;
    @property 
    def poles(this):
        return this.jetsolver.poles_res[0];
    
    def add_deform_Hz(this,data_z):                          
        csr=this.smbcoo.tocsr(data_z);        
        this.jetsolver.Hz+=csr.data;        
        return this;
    
    def __call__(this,rep=1):
        return this.jetsolver(rep=rep);
    
    def dump(this,N,rep=1):
        return this.jetsolver(N=N,rep=rep);
    
    @property    
    def xx(this):
        return this.jetsolver.xx;
    
    @xx.setter
    def xx(this,v):
        this.jetsolver.xx=v;    
         
    @property    
    def x(this):
        return this.jetsolver.x;
    
    @x.setter
    def x(this,v):
        this.jetsolver.x=v;
        
        
if __name__=='__main__':
    
    
    
    
    from utils import *
    
    from jet.tools import *
    import LIPA2.tools as l2ts
    from  LIPA2.qp_solver import *
    from FEM.trs import *
    from FEM.FEM_wave_fft import *
    
    from jsonrpc.json_io import *
    import os
    import scipy
    import scipy.sparse as sp
    
    norm=np.linalg.norm
        
    N=100000
    N=128*1024
    
    
    N=8*1024
    N=64*1024
    N=1024*1024
    N=4*64*1024
    N=1024*1024
    #N=7
    trs=make_trs_1D(N,fcycle=1)
    Np=trs.shape[0];
    k,m=FEM_KM_def()
    K0,M0=make_datas_FEM(Np)
    L=100;
    c=.5;
    c=1;
    dt=0.8;
    #dt=0.2;
    dx=L/Np
    
    M=(1+0)*M0*dx;
    K=(1+0)*K0/dx*c**2;
    
    
    LM=(8,8)
    LM=(6,8)
    LM=(12,14)
    LM=(10,12)
    LM=(8,8)
    #LM=(12,12)
    #FEM2data_coo(trs,datas)
    
    
    jet=FEM_jet_t(trs,Dmax=1)#,maskDiff=(0,2));
    jet2=FEM_jet_t(trs,Dmax=1,maskDiff=(0,2));
    
    
    
    
    jet.set_FEM_data([K,0.0*K,M])
    #jet.set_FEM_data([K,M])
    jet2.set_FEM_data([K,M],freal=1,dtype=np.double)
    
    mm=jet.jetsolver.spmb.mm;
    
    #mm2=[mm[0],None,mm[1]]
    lq=lipa_qp_base_t(mm,LM=LM,LU_factory=sp_LU_factory,nd=-1);
    '''
    tic()
    lq.reset(dt=dt);        
    mtoc('lq.reset:')
    tic()
    jet.reset(dt=dt,LM=LM);    
    mtoc('jet.reset:')
    '''
    tic()
    jet2.preset(dt=dt,LM=LM);    
    mtoc('jet2.preset:')
    tic()
    jet2.factorize();    
    mtoc('jet2.factorize:')
    
    r=np.linspace(-L/2,L/2-dx,N);
    gauss=lambda x,d: np.exp(-0.5*(x/d)**2)/(np.sqrt(2*np.pi)*d)
    
    gaussm=lambda x,d,m: np.exp(-0.5*(x/d)**m)    
    gaussm1 = lambda x,d,m: np.exp(-0.5*(x/d)**m)*(0.5*m)*x**(m-1)*(1./d)**m
    gr=gaussm(r,0.1,4);
    
    m=2;d=0.3
    d=0.5
    gr=gaussm(r,d,m);
    gr1=gaussm1(r,d,m);
    #gr=gaussm(r,0.15,2);
    #gr=gaussm(r,0.4,2);
    
    jet.x=gr;      
    lq.x=gr;
    
    jet2.x=gr;
    
    
    jet2.x=gr;
    jet2.xx[1]=-c*gr1;
    
    xx0=jet2.xx.copy()
    
    rep=3025
    '''
    tic() 
    lq(rep=rep)
    mtoc('lq rep:%d'%(rep))
    tic() 
    jet(rep=rep)
    mtoc('jet rep:%d'%(rep))
    '''
    tic()
    jet2(rep=rep)
    mtoc('jet2 rep:%d'%(rep))
    y=jet.x
    y2=jet2.x
    yy2=jet2.xx
    yy=jet.xx
    yl=lq.x
    yyl=lq.xx
    tic()    
    fw=FEM_wave_U_t(N,c=c,dx=dx)
    yyf=fw(xx0,dt*rep)
    yf=yyf[0];
    mtoc('FEM_wave_U_t')
    
    
    #print('errl=',norm(y-yl)/norm(yl))
    
    #raise SystemExit(0)
    import matplotlib.pyplot as plt

    #%matplotlib auto
    fig=plt.figure(figsize=(18,18))
    #fig.suptitle('sensors data + pure response ', fontsize=16)
    plt.grid(True,which='major')
    plt.minorticks_on()
    plt.grid(True,which='minor',alpha=0.2)
    
    #plt.plot([dt,dt],[-ma,ma],color='#777777')
    
    
    fn=lambda x: np.real(x)
    

    #plt.plot(r,fn(gr),label='$x_{0}$')
    #plt.plot(r,fn(yl),label='$x_{q}$')
    #plt.plot(r,fn(y),label='$x_{j}$')
    plt.plot(r,fn(y2),label='$x_{j2}$')
    plt.plot(r,fn(yf),label='$x_{fft}$')
    #legend = plt.legend(loc='upper right', shadow=True, fontsize='x-large')
    legend = plt.legend(loc='upper right', shadow=True, fontsize=20)
    plt.show()
    
    '''
    printf('errl=%0.2e%%\n',norm(y-yl)/norm(yl)*100)
    printf('errl2=%0.2e%%\n',norm(y2-yl)/norm(yl)*100)
    printf('err02=%0.2e%%\n',norm(y2-y)/norm(y)*100)
    '''
    printf('err_FFT=%0.2e%%\n',norm(y2-yf)/norm(yf)*100)



    
    
    
       