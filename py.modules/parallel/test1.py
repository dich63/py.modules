#
from parallel.sps import *
import numpy as np

#breakpoint()

A=np.array([[1,0],[0,1]])
C=np.array([[0,-1j],[1,0]])

cc=klu_common_t()
struct_def(cc,{'scale':11,'tol':1.134e-2})


s=sps_invoker(A+7*C,common={'scale':0,'tol':1e-2})

# calculate matrix pensils : z*C+A

N=5
yy=np.empty([N,2],dtype=np.complex128);

#yy[:,:]=np.random.randn(N,2)

# sequential
yy[:]=[33,1e-7]
for k in range(N):
    yy[k]=sps_invoker(A+k*C)(yy[k]);       
    
print('seq',yy)

# parallel if mode=1 ; sequential if mode=0
yy[:]=[33,1e-7]
pp=pp_group_ex(mode=1);
for k in range(N):       
    pp(sps_invoker(A+k*C,y=yy[k]).handle);

pp.join();

print('pp',yy)


    

    


