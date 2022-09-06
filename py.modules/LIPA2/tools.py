#
import numpy as np
import copy
#from numba import njit,jit

'''
def AzC1(xx0,CC,z,xxz_out,Czxx_out):
    
    nC=len(CC);
    nd=len(xx0);
    
    x=xxz_out[0];
    
    for k in range(1,nC):
        x=xx0[k]+z*x;   
        xxz_out[k]=x
        Ck=CC[k];
        if not Ck is None:
            Czxx_out[:]+=Ck.dot(x);
        
    for k in range(nC,nd):    
        xxz_out[k]=x=xx0[k]+z*x;
        
    return (xxz_out,Czxx_out)        
'''

def AzC0(xx0,DC,C1xx_out):  
     
    d=DC[1].dot(xx0[0]) if len(DC)>1 and not DC[1] is None   else 0+0j
    
    '''
    if type(C1xx_out) is np.ndarray:
        C1xx_out[:,:]=d;             
    '''      
    if type(C1xx_out) in (tuple,list):
        for c in C1xx_out:
            c[:]=d;
    else:
        C1xx_out[:]= d
        
pass    

def xxzjet(xx0,z,xxz_out):
    nD=len(xx0);   
    xxz_out[0]=xx0[0];
    for k in range(1,nD):
        xxz_out[k]=z*xx0[k]-xxz_out[k-1];
    return xxz_out;
    
#@njit
#@jit(nogil=True,fastmath=False)
def AzC1(xx0,DC,z,xxz_out,Czxx_out,buf=None):
    
    nC=len(DC);
    nD=len(xx0);   
    if buf is None:
        buf=np.empty(xx0.shape[-1],dtype=np.complex128)
    
    x=buf
    if nD>1:        
        
                
        x[:]=xx0[0]
        xxz_out[1]=-x;
        x[:]=xx0[1]+z*x;
        
        for k in range(2,nC):        
            Ck=DC[k];
            if not Ck is None:                
                Czxx_out[:]+=Ck.dot(x);
            
            if k<nD:
                xxz_out[k]=-x
                x[:]=xx0[k]+z*x;
            
        
        for k in range(nC,nD):
            xxz_out[k]=-x
            if k<(nD-1):
                x[:]=xx0[k]+z*x;                
                



def AzCn(xx0,DC,z,xxz_out,Czxx_out):
    
    nC=len(DC);
    nD=len(xx0);
    
    x=xx0[0];
    Czxx_out[:]=0
    
    for k in range(1,nC):
        
        Ck=DC[k];
        if not Ck is None:
            Czxx_out[:]+=Ck.dot(x);            
            
        if k<nD:
            xxz_out[k]=-x
            x=xx0[k]+z*x;        

    
    for k in range(nC,nD):
        xxz_out[k]=-x
        if k<(nD-1):
            x=xx0[k]+z*x;        
        
        
    return (xxz_out,Czxx_out)        


def AzC(xx0,CC,z,xxz_out=None,Czxx_out=None):
    
    if Czxx_out is None:
        Czxx_out=np.empty_like(xx0[0]);
    if xxz_out is None:
            xxz_out=np.empty_like(xx0);
    
    AzC0(xx0,CC,Czxx_out)
    AzC1(xx0,CC,z,xxz_out,Czxx_out)
    #return xxz_out,Czxx_out

def DCz0(DC,z):    
    
    Hz=copy.copy(DC[0]);
    f=False;
    for C in DC[1:]:
        if not C is None:            
                Hz+=z*C
        z*=z;
        
    return Hz;

def DCz(DC,z):            
    
    f=True;
    
    zz=1+0j;
    
    for C in DC:
        if not C is None:
            if f:
                f=False;
                Hz=zz*C                
            else:                
                Hz+=zz*C
                
        zz*=z;
        
    return Hz;


def get_zz(nd,z):
    zz=np.zeros((nd-1,1),dtype=np.complex128);
    zn=z;
    for k in range(nd-1):
        zz[k]=zn;
        zn*=z;        
    return zz;

    

    
def get_assemble(LM): 
    
    
    def assemble0(xx,xxz):
        xx[:]=-np.sum(xxz,0);
        return xx
    
    def assemble1(xx,xxz):
        xx-=np.sum(xxz,0);
        return xx

    def assemble2(xx,xxz):
        xx[:]=-xx-np.sum(xxz,0);
        return xx

    
    
    L,M=LM;        
    if L<M:
        assemble=assemble0
    elif M&1==0:
        assemble=assemble1
    else:
        assemble=assemble2
    
    return assemble    

def get_jet_shape(LM):
    for d in LM:
        try:
            shape=d.shape;
            return shape;
        except:
            pass
          
    raise Exception('Jet Error shape not found')
    
class keeper_t(object):
    def __init(self):
        self._buf=[];
    def __call__(self,v):
        self._buf.append(v);
        return v;