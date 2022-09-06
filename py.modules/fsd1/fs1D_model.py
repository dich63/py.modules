# -*- coding: utf-8 -*-
"""
Created on Sat Jun  4 14:28:04 2022

@author: wwww
"""
'''
from utils import *
from jsobj import *
import numpy as np

'''
from waves.fluid_solid_regions import *
from LIPA2.qp_solver import *
norm=np.linalg.norm

class fs1D_t(object):
    def __init__(self,model):
        self.model=reparse_mesh(model);
        
        
    def solve(self,dt,dz,Tm,ndmax=4,pade_nm=[4,6],fcycle=False):
        self.dt,self.dz,self.Tm,self.ndmax=dt,dz,Tm,ndmax;
        self._xxtt,self.solver=generate_signals(self.model,dt,dz,Tm,fcycle=fcycle,ndmax=ndmax,pade_nm=pade_nm)
        
    @property
    def xx(self):
        return self._xxtt;
    

    def frame(self,t,fs=0,nd=3):
        xxtt=self._xxtt
        dt=self.dt
        iz=int(t/dt);
        
    @property    
    def base_velocity(self):
        prms=to_dict(self.fullmodel.main_region.params);
        (Lc,alpha,M,eta,eta_p,kappa,rho,rho_f,rho_m)=make_sf_1D_params(**prms);
        return np.sqrt(Lc/rho);
        
       
    def signal_at(self,z,fs=0,nd=3):
        xxtt=self._xxtt
        Nz=xxtt.shape[0]
        iz=int(z/self.dz);
        if iz>=0 and iz<Nz:
            return xxtt[iz,nd,:,fs];
        else:
            raise Exception('region out of range')
    
    def spectrum_at(self,z,fs=0,nd=3):
        y=self.signal_at(z=z,fs=fs,nd=nd);
        return np.fft.rfft(y);
        
        
        

def generate_signals(model,dt,dx,Tm,zs=None,fcycle=False,ndmax=3,pade_nm=[4,6]):
    nd=ndmax+1;
    tic();
    mK,mM,mG,dx,Ne,model=model2FEM1D(model,dx=dx);
    toc('mK,mM,mG,dx,N,model=model2FEM1D:')
    tic();
    [smK,smG,smM]=make_sf_FEM([mK,mG,mM],fcycle=fcycle)
    toc('[smK,smG,smM]=make_sf_FEM([mK,mG,mM]):')
    solver=lipa_qp_base_t([smK,smG,smM],pade_nm,'sp_LU_factory',nd)
    tic();
    solver.reset(dt)
    toc('solver.reset(dt):')
    
    N=int(smK.shape[0]/2);
    
    f=np.zeros([N,2],dtype=np.complex);
    f[0,0]=1
    f[0,1]=0
    solver.x=0
    solver.f=f.reshape(-1)
    NT=int(Tm/dt+0.9999);
    tic();
    yy=solver.dump(NT);
    toc('solver.dump(NT):')
    xxtt=yy.reshape(N,nd,-1,2)
    if zs is None:
        return xxtt,solver
    else:
        izs=np.array(zs)/dx;
        izs=np.array(izs,dtype=int)
        izs[izs>N]=N-1;
        return xxtt[izs],solver;


if __name__=='__main__': 
    
    from skimage.transform import resize
    import k3d
    
    model=':file:'+__file__+'/../mod0.json'
    #model=':file:'+__file__+'/../model1D.json'
    '''
    xxtt=generate_signals(model,dt=1/100,dx=0.1,Tm=2,zs=[100,1000.1,1200])
    print('max(xxtt.reshape(-1))=',max(xxtt.reshape(-1)))
    '''
    fs=fs1D_t(model);
    fs.solve(dt=1/100,dz=0.1,Tm=1)
    
    
    import matplotlib.pyplot as plt
    fig=plt.figure(figsize=(12,8))
    #fig.suptitle('sensors data + pure response ', fontsize=16)
    plt.grid(True,which='major')
    plt.minorticks_on()
    plt.grid(True,which='minor',alpha=0.2)
    
    fs.solver.xx=0;
    x=fs.solver(0).reshape(5,-1,2);y=x[1,:,0]
    tt=np.arange(y.size)*fs.dz;
    
    nd=2;x=fs.solver(10).reshape(5,-1,2).real;plt.plot(tt,x[nd,:,0],tt,x[nd,:,1])

    pass



    