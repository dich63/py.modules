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



pepr.get((8,8));

crand=lambda *ls: np.random.randn(*ls)+1j*np.random.randn(*ls)

norm=np.linalg.norm
normm= lambda x:  norm(np.array(x).reshape(-1),ord=np.inf)


fn=r'O:\__ss\matrix\sfFEM4096k.json'
fn=r'O:\__ss\matrix\sfFEM1024k.json'
#fn=r'O:\__ss\matrix\FEM128k.json'
#fn=r'O:\__ss\matrix\sfFEM128k.json'
#fn=r'O:\__ss\matrix\sfFEM64k.json'
#fn=r'O:\__ss\matrix\sfFEM24k.json'
#fn=r'O:\__ss\matrix\sfFEM4k.json'
#fn=r'O:\__ss\matrix\KGM.json'
d=decode(fn,1)
smK,smG,smM=[m.tocoo() for m in [d.K,d.G,d.M]]
#smK,smG,smM=[m.tocsr() for m in [d.K,d.G,d.M]]
#smK[1000,222]=1e-4
#smK=0*smK
#sK,sG,sM=coo_disjunct_sum([smK,smG,smM])
sK,sG,sM=unicoo([smK,smG,smM])


CD=sparse2coo_matrix_batch(smK,smG,smM)
#smG=1e0*smK;
#CD=sparse2coo_matrix_batch(*coo_disjunct_sum([smK,smG,smM]))


#CD=sparse2coo_matrix_batch(smK,smM)
CDcsr=CD.tocsr(dtype=complex)
#CDscsr=CDs.tocsr(dtype=complex)

dt=.4
dt=1
dt=0.1

LM=(4,4)
Dmax=2
nd=Dmax+1

N=CDcsr.shape[0]
xx0=crand(nd,N);

tic()
lb=lipa_qp_base_t(CDcsr,LM=LM,nd=nd,LU_factory=sp_LU_factory)
mtoc('lipa_qp_base_t')
tic()
j=jet_csr_klu_t(CDcsr,Dmax=Dmax)
mtoc('jet_csr_klu_t')

tic()
lb.reset(dt)    
mtoc('lipa_qp_base_t::reset')

tic()
j.reset(dt=dt,LM=LM)
mtoc('jet_csr_klu_t::reset')
j.xx=xx0
lb.xx=xx0

lb();
j.step(1);
rep=1
tic()
#[lb() for k in range(rep)]
lb(rep)
toc('lipa_qp_base_t::step')
tic()
#[j.step(1) for k in range(rep)]
toc('jet_csr_klu_t::step')


err=normm(j.xx-lb.xx)/(normm(j.xx)+normm(lb.xx))
printf('err=%3.3g%% ',err*100)
derr=np.max(np.abs(derr2m(j.xx,lb.xx)))
printf('derr2m=%3.3g%%\n',derr*100)
derr2m(j.xx,lb.xx)
'''
M=LM[1]
S=5
xx=crand(M,S);
F=csr_matrix((N,S),dtype=complex)
F[:M,:]=xx
xt=xx.T
csr_gaxpy_NM(N,M,F.indptr,F.indices,F.data,xx,j.yz)
'''

