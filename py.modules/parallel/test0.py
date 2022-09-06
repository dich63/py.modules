#
from sparse import *
import numpy as np

sm=spm_context([[0,1],[1,0]])
si=sps_invoker(sm,y=np.array([1j,2]));
sr=hinvoke_batch((si.factorize,si.solve));

A=np.array([[1,0],[0,1]])
C=np.array([[0,-1],[1,0]])
y=[33,1e-7]
z=0;
# calculate matrix pensils : z*C+A
pp=pp_group_ex(mode=1);
N=5
yy=np.empty([N,2],dtype=np.complex128);
yy[:]=y

for k in range(N):       
    pp(sps_invoker(A+k*C,y=yy[k]).isolve);

pp.join();


print(yy)

    
    
    
    
    


