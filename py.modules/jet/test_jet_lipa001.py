#

from utils import *
from jet.jet_csr_iklu import *
from LIPA2.qp_solver import *
from FEM.FEM2sparse import *
from utils.derr2m import derr2m

from jsonrpc.json_io import *
import lipa.pade_exp_poles_res as pepr
from lipa_ctx.LIPA_ctx import *

def exp_pp(LM,g,t):
    
    def get_sp(LM):
        L,M=LM
        if L==M:
            if M&1==0:
                return 1;
            else:
                return -1
        else:
            return 0;

    po,re=pepr.get_poles_res_array(LM)
    pot,ret=po/t,re/t;
    
    #e=np.sum(ret/(g-pot))
    '''
    e=0;
    for m in range(len(po)):
        e+=ret[m]/(g-pot[m])
    '''
    e=np.sum(re/(g*t-po))
    
    return get_sp(LM)+e;

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





crand=lambda *ls: np.random.randn(*ls)+1j*np.random.randn(*ls)

norm=np.linalg.norm
normm= lambda x:  norm(np.array(x).reshape(-1),ord=np.inf)

N=2
E=np.eye(N,dtype=complex)*1000;


dt=.4

dt=.0000000022;g=(-.000500+4j)*8;g2=-.500+7j
dt=.022;g=(-0.10+4j)*77*1;g2=-.500+7j
dt=.0522e-0;g=(-0.051+5j)*77*0.4;g2=-.500+7j

M,mK=E,-g*E
#mK[1,1]=-g

M,K=coo_matrix(E),coo_matrix(mK)

sK,sM=unicoo([K,M])

#dt=0.3

LM=(2,2)
#LM=(2,4)


LM=(2,4)
LM=(6,12)
LM=(8,14)
LM=(8,10)
LM=(1,6)
LM=(6,14)
LM=(8,12)
LM=(3,7)
LM=(1,3)
LM=(10,14)        
    
    



Dmax=4
nd=Dmax+1

CD=sparse2coo_matrix_batch(sK,sM)
CDcsr=CD.tocsr(dtype=complex)

#
lb=lipa_qp_base_t([sK,sM],LM=LM,nd=nd,LU_factory=sp_LU_factory)
#lb=lipa_qp_base_t([sK.todense(),sM.todense()],LM=LM,nd=nd)
x=(1+0j)*np.ones(N)
xx0=np.tensordot([g**m for m in range(Dmax+1)],np.exp(0)*x,axes=0)

lb.reset(dt)    

lb.xx=xx0

j=jet_csr_klu_t(CDcsr,Dmax=Dmax)
j.reset(dt=dt,LM=LM)

po,re=j.poles_res
D=j.D
ZD=np.array(j.ZD,copy=True)
pp=(1+0j)*po;
#this.ZD=np.array([poles**d for d in range(1,D+1)]) 
for d in range(D):
    ZD[d]=pp;
    pp*=po;  

#j.ZD[:]=ZD

j.xx=xx0
[A,C]=[m.tocsc() for m in CD]

lctx=LIPA_solver_ctx(CD,dt=dt,pade_nm=LM,fcomplex=1,nd=nd)
#lctx=LIPA_solver_ctx([A,C],dt=dt,pade_mn=LM,fcomplex=1,nd=nd)
#lctx=LIPA_solver_ctx([sK,sM],dt=dt,pade_mn=LM,fcomplex=1,nd=nd)
lctx.xn=xx0



rep=100
lctx.step(rep)
lb(rep)
j(rep);
ep=exp_pp(LM,g,dt);epr=ep**rep
ex1=np.exp(1*g*dt);
ex1p=ex1**rep;

exa=np.tensordot([g**m for m in range(Dmax+1)],np.exp(rep*g*dt)*x,axes=0)
exa=np.tensordot([g**m for m in range(Dmax+1)],np.exp(rep*g*dt)*x,axes=0)
print('exp=\n',exa)
print('lb.xx=\n',lb.xx)
print(' j.xx=\n',j.xx)
print(' lctx.xn=\n',lctx.xx)
#raise SystemExit();
printf('errl %2.2e%%\n',normm(exa-lb.xx)/normm(exa)*100)
printf('erri %2.2e%%\n',normm(exa-j.xx)/normm(exa)*100)
printf('errx %2.2e%%\n',normm(exa-lctx.xx)/normm(exa)*100)

printf('\nderr %2.2e%%  |ex|=%e\n',normm(lb.xx-j.xx)/normm(exa)*100,normm(exa))

err2m=lambda x,y:np.abs(x-y)/np.abs(x)

poles,res=get_poles_res_array(LM)
po,re=j.poles_res
nn=Dmax
pr=np.sum(re*po**nn)

'''
pon=1
for n in range(len(po)):
    pon*=1/dt;
'''    
PR=np.sum(res*poles**nn)

perr=pr-PR*dt**-(nn+1)

printf(' perr = %e\n',np.abs(perr))
printf(' perrn = %e\n',np.abs(perr/pr))
printf(' dt^%d = %e\n',-(nn+1),dt**-(nn+1))

ee=err2m(exa,j.xx)
print('err2m(exa,j.xx):\n',ee)
print('np.abs(g*dt): ',np.abs(g*dt))
