#
import numpy as np
from scipy.linalg import expm
import LIPA2.tools as ts
from collections.abc import Iterable

def set_class_params(self,lcs,delname='self'):
    self.__dict__.update(lcs)    
    del self.self;
    return self

def list2array(L,dtype=np.complex128):
    if type(L) in (tuple,list):
        ls= [len(v) for v in L];        
        nL=np.zeros([len(L),np.max(ls)],dtype=dtype)
        for n,d,l in zip(nL,L,ls):
            n[:l]=d;
        return nL
    else:
        return L
    

class schz_matrix_t(object):
    def __init__(self,qp=None,t0=None,fperiod=False,dtype=np.complex128):    
        
        fqp= not qp is None
        
        if fqp:
            qp=list2array(qp,dtype=dtype)           
                                    
            ft0=not t0 is None            
            if ft0:
                tmp=t0
                t0=np.zeros(qp.shape[0],dtype=np.float64);
                t0[:len(tmp)]=tmp;
                
            del tmp;    
                
        
        set_class_params(self,locals());
                
    
        
    def reset(self,*lp,**kp):                                
        return self
    
    def __bool__(self):
        return self.fqp
    
    def __call__(self):
        
        return self.qp;
        
    def laplace(self,z,dtype=np.complex128):
        return schz_laplace_t(z,qp=self.qp,gs=self.g,dtype=dtype);
        
        
# #
#  list of   hev(t-t0[k])exp(g[k]*(t-t0[k]))*sum c_i* (t-t0[k])**i  
# #
class qp_matrix_t(object):
    def __init__(self,qp=None,g=None,t0=None,fperiod=False,dtype=np.complex128):  
        
        fqp= not qp is None
        
        if fqp:
            qp=list2array(qp,dtype=dtype)
            tmp=[] if g is None else g  
            g=np.zeros(qp.shape[0],dtype=dtype)
            g[:len(tmp)]=tmp;
                                    
            ft0=not t0 is None            
            if ft0:
                tmp=t0
                t0=np.zeros(qp.shape[0],dtype=np.float64);
                t0[:len(tmp)]=tmp;
                
            del tmp;    
                
        
        set_class_params(self,locals());
        
    
        
    def reset(self,t):        
        
            
        
        if self.fqp:
            L,K=self.qp.shape
            g=self.g;
            Nt=expm(np.diag(t*np.ones([K-1]),-1));                
            #G=np.eye(L) if g is None else np.diag(np.exp(t*np.array(g)))
            G=np.diag(np.exp(t*g))
            self.Nt,self.G=Nt,G
            
            if self.ft0:
                t0=self.t0;
                Nt0=np.zeros((L,K,K),dtype=np.float64)                
                #G0=np.ones(L) if g is None else np.exp((t-t0)*np.array(g))
                G0=np.exp((t-t0)*g)
                for l in range(L):                   
                    Nt0[l]=expm(np.diag((t-t0[l])*np.ones([K-1]),-1));
                    
                self.Nt0,self.G0=Nt0,G0
                    
            '''
            qp,g=self.qp,self.g;    
            Nt=expm(np.diag(t*np.ones([qp.shape[1]-1]),-1));                
            G=np.eye(qp.shape[0]) if g is None else np.diag(np.exp(t*np.array(g)))
            self.G,self.Nt=G,Nt;
            '''
        return self
    
    def __bool__(self):
        return self.fqp
    
    def __call__(self):       
        
        if self.fqp and (not self.fperiod ):
            qp=self.qp
            if not self.ft0:
                G,Nt=self.G,self.Nt
                qp=G@qp@Nt
                self.qp[:]=qp;        
            else:
                self.ft0=False;
                G0,Nt0=self.G0,self.Nt0
                L=G0.shape[0]
                for l in range(L):
                    qp[l]=G0[l]* qp[l]@Nt0[l];
                
            
        return self.qp;
        
    def laplace(self,z,dtype=np.complex128):
        return qp_laplace_t(z,qp=self.qp,gs=self.g,t0=self.t0,fperiod=self.fperiod,dtype=dtype);
        
  
class qp_laplace_t(object):
    
    def __init__(self,z,qp=None,gs=None,t0=None,fperiod=False,dtype=np.complex128):     
        
        
        zi=isinstance(z, Iterable)
        if zi:
            z=np.array(z);
            
            
            
        
        self.bz=bz=z**-1
        self.bzn=np.array([bz]);
        self.fqp=fqp= not ( qp is None);
        self.ft0=ft0= not ( t0 is None);
        self.fperiod=fperiod
        
        if not fqp:
            return;
        
        if gs is None:
            gs=np.zeros(qp.shape[1],dtype=dtype);
            
        #nz=len(gs)+1;   
        nz=qp.shape[1];      
        def getzz(x):
            pp=np.empty(nz,dtype=dtype);
            p=np.complex128(1.0)/x;
            for k in range(nz):
                pp[k]=p;
                p*=p;
            return pp
        
        self.qp=qp;
        
        sqp=qp.shape
        
        if zi:            
            self.pzs=pzs=np.full(z.shape+sqp, fill_value=np.nan,dtype=complex) 
            for m in range(len(z)):
                for pz,g in zip(pzs[m],gs):
                    pz[:]=getzz(z[m]-g);
        else:    
            self.pzs=pzs=np.full(sqp, fill_value=np.nan,dtype=complex)         
            for pz,g in zip(pzs,gs):
                pz[:]=getzz(z-g);
        if ft0:
            if zi:
                self.pzs0=pzs0=np.full(z.shape+sqp, fill_value=np.nan,dtype=complex)
                for m in range(len(z)):
                    for pz,pz0,dt in zip(pzs,pzs0,t0):
                        pz0[:]=np.exp(-dt*z[m])*pz;
            else:
                self.pzs0=pzs0=np.full(sqp, fill_value=np.nan,dtype=complex)                
                for pz,pz0,dt in zip(pzs,pzs0,t0):
                    pz0[:]=np.exp(-dt*z)*pz;
                    
                
                
    @property        
    def imz(self):        
        
        if self.ft0:
            self.ft0=self.fperiod
            return np.sum(self.pzs0*self.qp,-1);
        else:
            return np.sum(self.pzs*self.qp,-1);
        
            
        
    
    def __bool__(self):
        return self.fqp
        
    def __call__(self,F):
        #return np.sum(self.pzs*self.qp,1)
    
        if self.fqp and np.isfinite(self.qp[0,0]):            
            #return F@self.imz;
            return F@np.sum(self.pzs*self.qp,-1);
        else:
            return self.bz*F;
    
        
      
if __name__=='__main__':
    from utils import *
    import LIPA2.qp as qp0
    
    norm=np.linalg.norm
    normm= lambda x:  norm(x.reshape(-1),ord=np.inf)
    crand=lambda *ls: np.random.randn(*ls)+1j*np.random.randn(*ls)
    
    t=1;
    F=crand(1024,3)
    
    q=qp_matrix_t([[0,1],[0,0,3j+3,0,0,0,1],[0,0,0,1]],[0,2j,3])    
    q.reset(t);
    q()
    o=qp_laplace_t(1+10j,q.qp,q.g)
    q()
    #f=o([[11,11,1],[1,1,1],[1,1,1]])
    f=o(F)
    
    q0=qp0.qp_matrix_t([[0,1],[0,0,3j+3,0,0,0,1],[0,0,0,1]],[0,2j,3])    
    q0.reset(t);
    q0()    
    o0=qp0.qp_laplace_t(1+10j,q0.qp,q0.g)
    q0()
    #f0=o0([[11,11,1],[1,1,1],[1,1,1]])
    f0=o0(F)
    print('norm(f-f0)=',norm(f-f0))
    
    
    
    
    raise SystemExit()
    #q=qp_matrix_t([[1],[1,0,0,1],[1]],[1,2,3])
    
    q=qp_matrix_t([[0,1],[0,0,0,0,1],[0,0,0,1]],[0,2,3],[0,0.5,0.3]) 
    q.reset(t);
    qq0=q()
    qq1=q()
    
    q=qp_matrix_t([[0,1],[0,0,0,0,0,0,1],[0,0,0,1]],[0,2,3])    
    q.reset(t);
    q()
    ql=q.laplace(np.array([1+1j,2j]))
    ql=q.laplace(1+1j)
    q=qp_matrix_t([[0,1],[0,0,0,0,0,0,1],[0,0,0,1]])
    
    q.reset(t);
    o=qp_laplace_t(10j,q.qp,q.g)
    
    f=o([[11,11,1],[1,1,1],[1,1,1]])   
    q=qp_matrix_t([[0,1],[0,0,0,0,0,1],[0,0,0,1],[0,10,0,1]])
    q.reset(t);
    o=qp_laplace_t([10j,11,3],q.qp,q.g)
    ww=o.imz;
    f=o([[11,11,1],[1,1,1],[1,1,1]])   
    print(1/q.qp.real)
    a=q();
    print(1/q.qp.real)    
    
    '''
    o=q.laplace(10j)
    f=o([[11,11,1],[1,1,1],[1,1,1]])   
    print(1/q.qp.real)
    a=q();
    print(1/q.qp.real)    
    '''        
        
        
    