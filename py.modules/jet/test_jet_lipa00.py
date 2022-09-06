#

from utils import *
from jet.jet_csr_iklu import *
from LIPA2.qp_solver import *
from FEM.FEM2sparse import *
from utils.derr2m import derr2m

from jsonrpc.json_io import *
import lipa.pade_exp_poles_res as pepr

def unicoo(AA):
    
    def _zero_datas():
        return [ np.zeros_like(m.data) for m in AA ]
    
    nnzs=[m.nnz for m in AA]
    nnz=np.sum(nnzs);  
    
    
    #zdatas=[ np.zeros_like(m.data) for m in AA ]
    
    rows=np.concatenate([ np.array(m.row) for m in AA ])    
    cols=np.concatenate([ np.array(m.col) for m in AA ])
    rAA=[]
    for m in range(len(AA)):
        #d=zdatas.copy();
        d=_zero_datas();
        d[m][:]=AA[m].data;
        d=np.concatenate(d)
        rAA+=[sp.coo_matrix((d,(rows,cols)),shape=AA[m].shape)]
    
    return rAA;



pepr.get((6,6));

crand=lambda *ls: np.random.randn(*ls)+1j*np.random.randn(*ls)

norm=np.linalg.norm
normm= lambda x:  norm(np.array(x).reshape(-1),ord=np.inf)

N=2
E=np.eye(N,dtype=complex);

g=2j-0.00;
M,K=coo_matrix(E),coo_matrix(-g*E)
sK,sM=unicoo([K,M])

dt=.4
dt=.05
#dt=0.3

LM=(3,4)
Dmax=3
nd=Dmax+1

CD=sparse2coo_matrix_batch(sK,sM)
CDcsr=CD.tocsr(dtype=complex)

lb=lipa_qp_base_t([sK,sM],LM=LM,nd=nd,LU_factory=sp_LU_factory)
x=(1+0j)*np.ones(N)
xx0=np.tensordot([g**m for m in range(Dmax+1)],np.exp(0)*x,axes=0)

lb.reset(dt)    

lb.xx=xx0

j=jet_csr_klu_t(CDcsr,Dmax=Dmax)
j.reset(dt=dt,LM=LM)
j.xx=xx0
rep=100
lb(rep)
j(rep);
exa=np.tensordot([g**m for m in range(Dmax+1)],np.exp(rep*g*dt)*x,axes=0)
print('exp=\n',exa)
print('lb.xx=\n',lb.xx)
print(' j.xx=\n',j.xx)
#raise SystemExit();
printf('errl %2.2e%%\n',normm(exa-lb.xx)/normm(exa)*100)
printf('erri %2.2e%%\n',normm(exa-j.xx)/normm(exa)*100)
printf('derr %2.2e%%\n',normm(lb.xx-j.xx)/normm(exa)*100)