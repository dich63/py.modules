#
import numpy as np
import LIPA2.tools as ts
import lipa.pade_exp_poles_res as pepr  
import copy
import numpy.linalg
 

class lipa_poles_t(object):
    def __init__(self,iDC,ixx0,ixxz,ifz,iCzxx,unwrap=lambda x: x):
        
        self.DC,self.xx0,self.xxz,self.fz,self.Czxx=\
            [unwrap(p) for p in [iDC,ixx0,ixxz,ifz,iCzxx]]
        
        #self.init(z,Bz);
        
    
    def reset(self,poles_res,LU_factory):
        
        self.solver=None;
        
        z,Bz=poles_res;
        self.z,self.Bz,self.LU_factory=z,Bz,LU_factory;
        
        
        
        Hz=ts.DCz(self.DC,z);
        
        self.solver=LU_factory(Hz);       
        
        self.zz=ts.get_zz(self.xx0.shape[0],z);
        
        #self.Czxx=np.zeros_like(self.xx0[0]);
        self.xxz[:]=np.zeros_like(self.xx0);
        
    
        return True
        
    def __call__(self):
        #AzC(xx0,CC,z,xxz_out=None,Czxx_out=None):
        xxz,Czxx,Bz,zz=self.xxz,self.Czxx,self.Bz,self.zz;
        
        #ts.AzC(self.xx0[:-1],self.DC[1:],self.z,xxz[1:],Czxx)
        xxz[:]=np.nan;
        
        #ts.AzC0(self.xx0,self.DC,Czxx);        
        ts.AzC1(self.xx0,self.DC,self.z,xxz,Czxx);
        
        '''
        ts.AzCn(self.xx0,self.DC,self.z,xxz,Czxx)
        
        
        xxy=np.nan*xxz
        Cxxx=np.nan*Czxx
        
        ts.AzC0(self.xx0,self.DC,Cxxx);        
        ts.AzC1(self.xx0,self.DC,self.z,xxy,Cxxx);
        '''
        
        Czxx+=self.fz;
        
        xxz[0][:]=xz=self.solver(Czxx);
        
        if xxz.shape[0]>1:
            #xxz[1:]=np.kron(zz,xz)-xxz[1:];
            xxz[1:][:]+=np.kron(zz,xz);
            
        
        xxz*=Bz;
        
                    
        return True


def LU_factory_def(A):                   
    
    invA=numpy.linalg.pinv(A);    
    return lambda x: invA.dot(x) ;

    
    
class lipa_base_t(object): 
    
    def __init__(self,DC,LM=[2,2],LU_factory=LU_factory_def,nd=-1):
        
        LM =[int(v) for v in LM]
        
        ndL=len(DC)-1
        if nd<ndL:
            nd=ndL;
            
        
        
        self.poles_res=prs=pepr.get(LM);
        Np=len(prs);
        
        self.assemble=ts.get_assemble(LM);
        
        self.LU_factory=LU_factory;
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
        
        self.psolvers=[ lipa_poles_t(DC,xx,xxz[k],fz[k],Czxxs[k])  for k in range(Np) ];
        
        
    
        
        
    def reset(self,dt):
        
        prs=self.poles_res;
        self.dt=dt;
        f=self.f;
        LU_factory=self.LU_factory
        
        pdt=[]
        
        for sol,pr,fz in zip(self.psolvers,prs,self.fz):
            
            p,r=pr[0]/dt,pr[1]/dt            
            pdt+=[p]            
            sol.reset((p,r),LU_factory); 
            
        self.pdt=pdt;    
        
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
        
        return xx;
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
        
    
        


def lipa_dense(DC,LM=[2,2],nd=-1):
    def toa(x):
        return np.array(x,dtype=np.complex128)
    
    return lipa_base_t(toa(DC),LM,nd=nd);

def lipa_number(DC,LM=[2,2],nd=-1):    
    def toa(x):
        return np.array(x,dtype=np.complex128).reshape([1,1]);    
    
    return lipa_base_t([toa([c]) for c in DC ],LM,nd=nd);

if __name__=='__main__':
    
    from utils import *
    '''
    d=7;x=10
    lb=lipa_number([d,1,1],[7,8],nd=2).reset(1)
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
    nd=5
    lb=lipa_base_t(np.array((D,C)),LM=[5,5],nd=nd).reset(dt);
    lb.xx=xx[0:nd];
    
    zz=lb()
    
    print(zz)
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