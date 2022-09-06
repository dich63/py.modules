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


#def generate_signals(model,dt,dx,Tm,zs=None,fcycle=False,ndmax=3,pade_nm=[4,6]):
ndmax=3    
dx=0.001
model=model_expand(':file:'+__file__+'/../model1D.json')    
nd=ndmax+1;
tic();
mK,mM,mG,dx,Ne,model=model2FEM1D(model,dx=dx);

toc('mK,mM,mG,dx,N,model=model2FEM1D:')
tic();
[smK,smG,smM]=make_sf_FEM([mK,mG,mM],fcycle=0,fmt='csr')


toc('[smK,smG,smM]=make_sf_FEM([mK,mG,mM]):')

s=smK.tocsr()
nnz=s.nnz

'''
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


if __name__=='__main__': 
    
    model=':file:'+__file__+'/../mod0.json'
    xxtt=generate_signals(model,dt=1/10,dx=0.1,Tm=2,zs=[100,1000.1,1200])
    print('max(xxtt.reshape(-1))=',max(xxtt.reshape(-1)))
    pass

'''    