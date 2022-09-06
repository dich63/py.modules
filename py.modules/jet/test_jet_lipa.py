#

from utils import *
from jet.jet_csr_iklu import *
from LIPA2.qp_solver import *
from FEM.FEM2sparse import *
from utils.derr2m import derr2m
norm=np.linalg.norm
normm= lambda x:  norm(np.array(x).reshape(-1),ord=np.inf)



x=np.array(range(12));
xx=x.reshape(3,-1)
a=0.5j
C1=-np.array([[1,0],[0,-1]],dtype='complex')
C=np.array([[1,0],[0,1]],dtype='complex')
D=np.array([[1*a,0],[0,1*a]],dtype='complex')

xx=np.array([[1,-1],[-a,a],[a**2,-a**2],[-a**3,a**3],[a**4,-a**4]],dtype='complex')
xx0=copy.copy(xx)
f=np.array([0,0],dtype='complex')

dt=.4
dt=1
nd=5
LM=[8,10]
#lb=lipa_qp_base_t(np.array((D,C)),LM=LM,nd=nd,LU_factory=sp_LU_factory).reset(dt);
DC=(D,C)

lb=lipa_qp_base_t(np.array(DC),LM=LM,nd=nd).reset(dt);
lb.xx=xx[0:nd];
rep=100
zz=lb(rep)

#print(zz.T.real)
print('pp',list(zz.T[0]))
dtr=dt*rep
xex=np.array([np.exp(-a*dtr),-a*np.exp(-a*dtr),a**2*np.exp(-a*dtr),-a**3*np.exp(-a*dtr),a**4*np.exp(-a*dtr)])
print('ex',xex)

CD=sparse2coo_matrix_batch(*DC)
CDcsr=CD.tocsr(dtype=complex)

j=jet_csr_klu_t(CDcsr,Dmax=nd-1)
#%timeit
j.reset(LM=LM,dt=dt)
j.xx=xx
j(rep)

print('nb',list(j.xx.T[0]))

err=normm(j.xx-lb.xx)/(normm(j.xx)+normm(lb.xx))
printf('err=%3.3g%% ',err)
derr=np.max(np.abs(derr2m(j.xx,lb.xx)))
printf('derr2m=%3.3g%%\n',derr*100)

err=normm(j.xx-lb.xx)/(normm(j.xx)+normm(lb.xx))
printf('err=%3.3g%% ',err)
derr=np.max(np.abs(derr2m(j.xx,lb.xx)))
printf('derr2m=%3.3g%%\n',derr*100)

derrj=np.max(np.abs(derr2m(j.xx.T[0],xex)))
printf('derrj=%3.3g%%\n',derrj*100)
derrq=np.max(np.abs(derr2m(lb.xx.T[0],xex)))
printf('derrq=%3.3g%%\n',derrq*100)