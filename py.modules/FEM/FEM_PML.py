# -*- coding: utf-8 -*-
"""
Created on Sat Aug 27 19:04:18 2022

@author: DICH
"""
import numpy as np

def create_profile(N,h,L=1,xc=0,win=np.hanning):
    Nw=int((N*h)/L)
    nc=int((N*xc)/L)
    fw=win(Nw)
    Nwl=Nw//2
    Nwr=Nw-Nwl
    f=np.zeros(N,dtype=np.float64);
    i=np.arange(N+nc-Nwl,N+nc+Nwr)% N;    
    f[i]=fw
    return f,i


def create_PML_profile(zz,profile):    
    
    pz= [profile/(profile+z) for z in zz]
    #pz= [profile*z for z in zz]
    pz=np.array(pz,dtype=complex);
    PML_p= -pz*(2-pz);    
    
    PML_p=np.array(pz,dtype=complex);
    
    return PML_p;
    

def make_PML(PMLz,K):
    # PMLz --> (m,n)     
    # D.shape --> (n,k,j)
    # rout -->(m,n,k,j)   
    rout=np.einsum('mn,nkj ->mnkj',PMLz,K)
    return rout

def create_PML_data(zz,K,h,xc=0,L=1,L0=0,df=1.0e11,win=np.hanning):
    N=K.shape[0];
    pf,ii=create_profile(N,h,L=L-L0,xc=xc-L0,win=win);
    pf*=df;
    pml=create_PML_profile(zz,pf)
    Kpml=make_PML(pml,K)    
    return Kpml,ii;
    
    

if __name__=='__main__':
    
    from FEM.FEM2sparse import *
    from FEM.trs import *
    from utils import *
    from FEM.FEM_jet import *
    from FEM.FEM_wave_fft import *
    import scipy.signal.windows as win
    norm=np.linalg.norm
    '''
    N=64
    L=16
    h=4
    pf=create_profile(N,h,L,win=np.hanning)#lambda x: 1)
    zz=[1,1j,10]
    pml=create_PML_profile(zz,pf)
    
    K0,M0=make_datas_FEM(N)
    
    Kp=make_PML(pml,K0)
    '''
    N=100000
    N=128*1024
    
    
    
    N=64*1024
    N=1024*1024
    N=4*64*1024
    N=8*1024
    N=4*64*1024
    #N=7
    trs=make_trs_1D(N,fcycle=1)
    Np=trs.shape[0];
    k,m=FEM_KM_def()
    K0,M0=make_datas_FEM(Np)
    L=100;
    
    c=.5;
    c=1;
    mf=5
    dt=0.5/mf;
    
    
    dx=L/Np
    
    M=(1+0)*M0*dx;
    K=(1+0)*K0/dx*c**2;
    
    
    LM=(8,8)
    LM=(6,8)
    LM=(12,14)
    LM=(10,12)
    LM=(4,4)
    #LM=(8,8)
    #LM=(4,6)
    #FEM2data_coo(trs,datas)
    
    
    
    jet=FEM_jet_t(trs,Dmax=1,maskDiff=(0,2));
    
    jet.set_FEM_data([K,M],freal=1,dtype=np.double)
    
    mm=jet.jetsolver.spmb.mm;

    tic()
    jet.preset(dt=dt,LM=LM);    
    #def create_PML_data(zz,K,h,xc=0,L=1,L0=0,df=1.0,win=np.hanning):
    zz=jet.poles;
    h=24;
    #wif=lambda M: win.exponential(M, 0, -(M-1) / np.log(0.01), False)
    wif=lambda M: win.exponential(M, None, -(M-1) / np.log(1e-12),sym= True);wif(111)
    pmlHz,ii=create_PML_data(zz,K,h,xc=-30,L=L/2,L0=-L/2,df=1e1,win=wif)#win.triang);    
    jet.add_deform_Hz(pmlHz)
    mtoc('jet.preset:')
    tic()
    jet.factorize();
    mtoc('jet.factorize:')
    
    r=np.linspace(-L/2,L/2-dx,N);
    gauss=lambda x,d: np.exp(-0.5*(x/d)**2)/(np.sqrt(2*np.pi)*d)
    
    gaussm=lambda x,d,m: np.exp(-0.5*(x/d)**m)    
    gaussm1 = lambda x,d,m: np.exp(-0.5*(x/d)**m)*(0.5*m)*x**(m-1)*(1./d)**m
    
    
    m=2;d=1.3
    gr=gaussm(r,d,m);
    gr1=gaussm1(r,d,m);
    
    jet.x=gr;
    jet.xx[1]=-c*gr1;
    xx0=jet.xx.copy()
    rep=45*mf
    #lq();jet();jet2()
    
    tic() 
    jet(rep=rep)
    mtoc('jet rep:%d'%(rep))
    
    
    y=jet.x
    #y[ii]=np.nan
    tic()    
    fw=FEM_wave_U_t(N,c=c,dx=dx)
    yyf=fw(xx0,dt*rep)
    yf=yyf[0];
    mtoc('FEM_wave_U_t')
    
    #y[np.abs(y)>max(np.abs(yf))]=np.nan;
    import matplotlib.pyplot as plt

    #%matplotlib auto
    %matplotlib inline
    fig=plt.figure(figsize=(18,18))
    #fig.suptitle('sensors data + pure response ', fontsize=16)
    plt.grid(True,which='major')
    plt.minorticks_on()
    plt.grid(True,which='minor',alpha=0.9)
    
    #plt.plot([dt,dt],[-ma,ma],color='#777777')
    
    
    fn=lambda x: np.real(x)
    

    plt.plot(r,fn(y),label='$x_{j2}$')
    plt.plot(r,fn(yf),label='$x_{fft}$')
    #legend = plt.legend(loc='upper right', shadow=True, fontsize='x-large')
    legend = plt.legend(loc='upper right', shadow=True, fontsize=20)
    plt.show()
    