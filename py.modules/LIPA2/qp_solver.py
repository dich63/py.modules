#
import numpy as np
import LIPA2.tools as ts
import LIPA2.qp as qp
import lipa.pade_exp_poles_res as pepr  
import copy
import numpy.linalg
try:
    import scipy.sparse as sp
    import scipy.sparse.linalg
except:
    pass
    
 

class lipa_qp_poles_t(object):
    
    def __init__(self,iDC,ixx0,ixxz,ifz,iCzxx,iFF=None,iqp=None,ig=None,unwrap=lambda x: x):
        
        
        self.DC,self.xx0,self.xxz,self.fz,self.Czxx,self.FF,self.qp,self.g=\
            [unwrap(p) for p in [iDC,ixx0,ixxz,ifz,iCzxx,iFF,iqp,ig]]
            
        
        #self.init(z,Bz);
        
    
    def reset(self,poles_res,LU_factory):
        
        self.solver=None;
        
        z,Bz=poles_res;
        self.z,self.Bz,self.LU_factory=z,Bz,LU_factory;
        
        self.qp_laplace=qp.qp_laplace_t(z,self.qp,self.g);
        
        
        self.Hz=Hz=ts.DCz(self.DC,z);
        
        self.solver=LU_factory(Hz);       
        
        self.zz=ts.get_zz(self.xx0.shape[0],z);
        
        #self.Czxx=np.zeros_like(self.xx0[0]);
        self.xxz[:]=np.zeros_like(self.xx0);
        
    
        return True
        
    def __call__(self):
        
        xxz,Czxx,Bz,zz,qp_laplace=self.xxz,self.Czxx,self.Bz,self.zz,self.qp_laplace;
        
        
        #xxz[:]=np.nan;
        
        
        ts.AzC1(self.xx0,self.DC,self.z,xxz,Czxx);
        
        
        fz=qp_laplace(self.FF) if qp_laplace.fqp else self.fz;            
        # set breakpoint
        #import pdb; pdb.set_trace()

        Czxx+=fz;
        
        xxz[0][:]=xz=self.solver(Czxx);
        
        if xxz.shape[0]>1:        
            xxz[1:][:]+=np.kron(zz,xz);
            
        
        xxz*=Bz;
        
                    
        return True


def LU_factory_def(A):                   
    
    invA=numpy.linalg.pinv(A);    
    return lambda x: invA.dot(x) ;

def sp_LU_factory(A):    
    
    options=dict(Equil=True,PivotGrowth=True
                                     ,PrintStat=False
                                     ,DiagPivotThresh=0.1
                                     ,ColPerm='MMD AT PLUS A'#'MMD_ATA'#
                                     ,ConditionNumber=False,IterRefine='NOREFINE')
    #options={}
    #lu=sp.linalg.splu(A,options=options)
    #
    lu=sp.linalg.splu(A)
    #lu=sp.linalg.splu(A,options=dict(Equil=False, IterRefine='SINGLE'))
    return lambda x: lu.solve(x) ;

def tocscifsp(A):
    if hasattr(A,'tocsc'):
        return A.tocsc();
    else:
        return A;
    
def get_factory(f):
    if type(f)==str:
        return eval(f)
    else:
        return f;
                       
    
class lipa_qp_base_t(object): 
    
    def __init__(self,DC,LM=[2,2],LU_factory=LU_factory_def,nd=-1,FF=None,QP=None,g=None):
        
        LM =[int(v) for v in LM]
        
        DC=[tocscifsp(m) for m in DC ]
        
        ndL=len(DC)-1
        if nd<ndL:
            nd=ndL;        
        
        
        self.poles_res=prs=pepr.get(LM);
        Np=len(prs);
        
        self.assemble=ts.get_assemble(LM);
        
        self.LU_factory=get_factory(LU_factory);
        self.DC=DC;
        
        shape=ts.get_jet_shape(DC);
        
        self._xx=xx=np.zeros((nd,shape[0]),dtype=np.complex128);
        self._f=np.zeros((shape[0],),dtype=np.complex128);
        self.psolvers=[];
        shape=(Np,)+xx.shape;        
        self.xxz=xxz=np.empty(shape,dtype=np.complex128);
        xxz[:,:]=np.nan;
        shape=(Np,)+xx[0].shape;
        self.fz=fz=np.zeros(shape,dtype=xx.dtype);
        
        #self.Czxxs=Czxxs=[np.zeros(xx[0].shape,dtype=np.complex128) for k in range(Np)];
        
        self.Czxxs=Czxxs=np.zeros((Np,xx[0].shape[0]),dtype=np.complex128);
        
        
        self.qp_matrix=qp.qp_matrix_t(QP,g);
        
        QP=self.qp_matrix.qp
        
        self.FF=FF;
        
        self.psolvers=[ lipa_qp_poles_t(DC,xx,xxz[k],fz[k],Czxxs[k],FF,QP,g)  for k in range(Np) ];
        
        self.fpulse=False
        
    
        
        
    def reset(self,dt):
        
        prs=self.poles_res;
        self.dt=dt;
        f=self.f;
        LU_factory=self.LU_factory
        
        self.qp_matrix.reset(dt);
        
        
        pdt=[]
        rdt=[]
        
        for sol,pr,fz in zip(self.psolvers,prs,self.fz):
            
            p,r=pr[0]/dt,pr[1]/dt            
            pdt+=[p]            
            rdt+=[r]
            sol.reset((p,r),LU_factory); 
            
        self.pdt=pdt;    
        self.rdt=rdt;
        
        self.f=self.f;
        
        return self;
    
    def dump(self,N,rep=1):
        xx=self._xx
        xxn=np.empty((N,)+xx.shape,dtype=xx.dtype);
        for n in range(N):
           xxn[n]=self(rep);
        return xxn;

    
        
    def __call__(self,rep=1):        
        
        for r in range(rep):
            ts.AzC0(self.xx,self.DC,self.Czxxs);
            
            for sol in self.psolvers:
                sol();
            xx,xxz=self._xx,self.xxz;      
            self.assemble(xx,xxz);
            self.qp_matrix();
            #self._pulse_reset();
            
        return self._xx
        '''
        xq=1*xx;
        r=self.assemble(xx,xxz);
        #r[1]-=xq[0];
        return r;
        '''
    
    @property
    def x(self):
        return self._xx[0];
    @x.setter
    def x(self,v):
        self._xx[1:]=0;
        self._xx[0]=v;
        
    @property    
    def xx(self):
        return self._xx;
    
    @xx.setter
    def xx(self,v):
        self._xx[:]=v;
    
    @property
    def f(self):
        return self._f;
    
    @f.setter
    def f(self,v):
        
        f=self._f
        
        if id(v) !=id(f):
            f[:]=v;           
                
        for p,fz in zip(self.pdt,self.fz):                        
            fz[:]=f/p;
        
    def pulse(self,f,m=0):
        if m:
            for p,fz in zip(self.pdt,self.fz):                     
                fz[:]=f*(p**m-1);
        else:
            for fz in self.fz:
                fz[:]=f;
            
        self.fpulse=True;
            
    def _pulse_reset(self):
        if self.fpulse:
            self.fpulse=False;
            for fz in self.fz:                        
                fz[:]=0;


def lipa_qp_dense(DC,LM=[2,2],nd=-1):
    def toa(x):
        return np.array(x,dtype=np.complex128)
    
    return lipa_qp_base_t(toa(DC),LM,nd=nd);

def lipa_qp_number(DC,LM=[2,2],nd=-1,FF=None,qp=None,g=None):    
    def toa(x):
        return np.array([x],dtype=np.complex128).reshape([1,1]) if not x is None else None;    
    
    return lipa_qp_base_t([toa(c) for c in DC ],LM,nd=nd,FF=FF,QP=qp,g=g);

if __name__=='__main__':
    
    from utils import *
    '''
    d=7;x=10
    lb=lipa_qp_number([d,1,1],[7,8],nd=2).reset(1)
    lb.x=[x]
    yy=lb()
    print(yy)
    print([x*np.exp(-d*1),-d*np.exp(-d*1)*x])
    '''
    #raise SystemExit(0)
    
    x=np.array(range(12));
    xx=x.reshape(3,-1)
    a=0.5
    C1=np.array([[1,0],[0,1]],dtype='complex')
    C=np.array([[1,0],[0,1]],dtype='complex')
    D=np.array([[1*a,0],[0,1*a]],dtype='complex')
    
    xx=np.array([[1,-1],[-a,a],[a**2,-a**2],[-a**3,a**3],[a**4,-a**4]],dtype='complex')
    xx0=copy.copy(xx)
    f=np.array([0,0],dtype='complex')
    
    dt=.4
    dt=10
    nd=5
    lb=lipa_qp_base_t(np.array((D,C)),LM=[8,8],nd=nd).reset(dt);
    lb.xx=xx[0:nd];
    
    zz=lb()
    
    #print(zz.T.real)
    print('pp',list(zz.T[0].real))
    print('ex',[np.exp(-a*dt),-a*np.exp(-a*dt),a**2*np.exp(-a*dt),-a**3*np.exp(-a*dt),a**4*np.exp(-a*dt)])
    
    raise SystemExit(0)
    tic()
    for k in range(1000):
        zz=lb()
    toc('')
    
    '''
    print(lb())
    print(np.exp(-2*dt))
    
    
    tic();
    for k in range(1000):
        lb();
    t=toc('')
    '''