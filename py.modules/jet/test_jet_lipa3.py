#

from utils import *
from jet.jet_csr_iklu import *
from LIPA2.qp_solver import *
from FEM.FEM2sparse import *
from utils.derr2m import derr2m

from jsonrpc.json_io import *
import lipa.pade_exp_poles_res as pepr
#from parallel.sparse import *
from lipa_ctx.LIPA_ctx import *


pepr.get((8,8));

crand=lambda *ls: np.random.randn(*ls)+1j*np.random.randn(*ls)

norm=np.linalg.norm
normm= lambda x:  norm(np.array(x).reshape(-1),ord=np.inf)


fn=r'O:\__ss\matrix\sfFEM4096k.json'
fn=r'O:\__ss\matrix\sfFEM1024k.json'
#fn=r'O:\__ss\matrix\sfFEM128k.json'
fn=r'O:\__ss\matrix\sfFEM64k.json'
#fn=r'O:\__ss\matrix\sfFEM24k.json'
#fn=r'O:\__ss\matrix\sfFEM4k.json'
#fn=r'O:\__ss\matrix\KGM.json'
#fn=r'O:\__ss\matrix\AC_3E.json'
d=decode(fn,1)
'''
smK,smG,smM=[m.tocoo() for m in [d.K,d.G,d.M]]
CD=sparse2coo_matrix_batch(smK,smG,smM)
'''
smK,smG,smM=[m.tocoo() for m in [d.K,d.K,d.M]]
CD=sparse2coo_matrix_batch(smK,smG,smM)
CDcsr=CD.tocsr(dtype=complex)
CDcsc=[m.tocsc() for m in CDcsr]




dt=.4
dt=0.001
#dt=0.01

LM=(8,8)
Dmax=1
nd=Dmax+1

N=CDcsr.shape[0]
xx0=crand(nd,N);
xx0=crand(nd,N);
xx0=np.linspace(-1,1j,N)
tic()
lb=lipa_qp_base_t(CDcsr,LM=LM,nd=nd,LU_factory=sp_LU_factory)
mtoc('lipa_qp_base_t')
tic()
j=jet_csr_klu_t(CDcsr,Dmax=Dmax)
mtoc('jet_csr_klu_t')




tic()
lctx=LIPA_solver_ctx(CDcsc,dt=dt,pade_nm=LM,nd=nd,fcomplex=1)
mtoc('LIPA_solver_ctx[reset]')


tic()
lb.reset(dt=dt)
mtoc('lipa_qp_base_t::reset')


tic()
j.reset(dt=dt,LM=LM)
mtoc('jet_csr_klu_t::reset')



j.xx=xx0
lb.xx=xx0
lctx.xn=xx0

tic()
lb()
mtoc('lipa_qp_base_t::step')

tic()
j.step()
mtoc('jet_csr_klu_t::step')

tic()
lctx.step(1)
mtoc('lctx.step')



err=normm(j.xx-lb.xx)/(normm(j.xx)+normm(lb.xx))
err2=normm(j.xx-lctx.xx)/(normm(j.xx)+normm(lctx.xx))
printf('err=%3.3g%% ',err*100)
printf('err2=%3.3g%% ',err2*100)
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

