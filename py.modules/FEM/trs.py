#
import numpy as np
#import numba

def ext_def(opts,**dflt):    
    t={}
    t.update(dflt)   
    t.update(opts)
     
    return t;



def _make_trs(Nt,C,trs,t0):
    for n in range(Nt):
        trs[n]=t0;
        t0+=C; 
    pass

try:
    import numba
    _make_trs=numba.njit(_make_trs);
except:    
    pass


    

def make_trs_1D(Np,C=1,fcycle=False,dtype=np.int32):
    
    Nt=Np-1;
    
    
    if fcycle:
        
        trs=np.empty((Np,2*C),dtype=dtype);        
        #trs[Nt]=[2*Nt,2*Nt+1,0,1];
        ac=np.arange(C)        
        trs[Nt][C:]=ac
        trs[Nt][:C]= C*Nt+ ac
        
    else:
        #Np=Nt+1;
        trs=np.empty((Nt,2*C),dtype=dtype);
        
    #t0=np.array([0,1,2,3],dtype=dtype)    
    t0=np.arange(2*C,dtype=dtype)    
    '''
    for n in range(Nt):
        trs[n]=t0;
        t0+=C; 
    '''
    _make_trs(Nt,C,trs,t0);
    
    return trs;


def FEM_KM_def():
    return [
        np.array([[1,-1],[-1,1]],dtype=np.float64),
        np.array([[3,1],[1,3]],dtype=np.float64)/8
           ]


def make_datas_FEM(N,dtype=np.float64,**ls):
    
    if len(ls)==0:
        ls=FEM_KM_def();
    else:
        ls=[np.array(l,dtype=dtype) for l in ls ];
        
    sh=ls[0].shape;    
    M=len(ls);    
    datas=np.empty((M,N)+sh,dtype=dtype);
    for m in range(M):
        datas[m][:]=ls[m];
        
    return datas
        
def make_trs2datas_FEM(trs,dtype=np.float64,**ls):
    return make_datas_FEM(len(trs),dtype=dtype,**ls)
    
    
    
    
