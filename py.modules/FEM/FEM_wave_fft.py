#from FEM.FEM_jet import *
import numpy as np

pi2=2*np.pi;
pi2j=1j*pi2;

def FEM_G(N,c):
    G=np.empty(N,dtype=complex);
    
    for n in range(N):
        k=(pi2*n)/N;
        cosk=np.cos(k)
        Dk=8*(1-cosk)/(3+cosk)
        G[n]=c*np.sqrt(Dk)
    return G;

def FDM_G(N,c):
    G=np.empty(N,dtype=complex);
    
    for n in range(N):
        k=(pi2*n)/N;
        cosk=np.cos(k)
        Dk=2*(1-cosk)
        G[n]=c*np.sqrt(Dk)
    return G;

        
def native_G(N,c):
    G=np.empty(N,dtype=complex);
    N2=N//2;
    for n in range(N):
        n2= n if n<N2 else N-n;
        k=(pi2*n2)/N;        
        G[n]=c*k
    return G;



class FEM_wave_U_t(object):
    def __init__(this,N,c=1.0,method=FEM_G,dx=1.0):
        this.N=N
        c=c/dx;
        this.g=g=method(N,c);
        this.G=G=np.empty((N,2,2),dtype=complex);
        this.iG=iG=np.empty((N,2,2),dtype=complex);
        
        for n in range(N):
            gn=g[n]
            '''
            G[n]=[[1,1],[1j*gn,-1j*gn]]
            ign=[[0.5,-0.5j/gn],[0.5,0.5j/gn]] if np.abs(gn)>0 else [[0.5,0],[0.5,0]]
            '''
            if np.abs(gn)>0:
                G[n]=[[1,1],[1j*gn,-1j*gn]]
                ign=[[0.5,-0.5j/gn],[0.5,0.5j/gn]]
                iG[n]=ign
            else:                
                iG[n]=G[n]=[[1.0, 0],[0,1.0]]
                

            
            
    def exp(this,t):        
        N,g,G,iG=this.N,this.g,this.G,this.iG
        eU=np.zeros((N,2,2),dtype=complex);
        eU[:,0,0]=np.exp(1j*t*g)
        eU[:,1,1]=eU[:,0,0].conj();
        
        for n in range(N):
            if np.abs(g[n])>0:
                eU[n]=G[n]@eU[n]@iG[n]        
            else:
                eU[n]=[[1,t],[0,1.0]]
        
        this.expU=eU; 
        return this
        
        
    def make(this,xx1):
        ff1=np.fft.fft(xx1,axis=1).T
        Ut=this.expU;
        for n in range(this.N):
            ff1[n]=Ut[n]@ff1[n]
            
        return np.fft.ifft(ff1.T,axis=1)
        
    
    def __call__(this,xx1,t=1):
        return this.exp(t).make(xx1);
        
if __name__=='__main__':
    
    N=1024*256
    L=100
    dx=L/N
    r=np.linspace(-L/2,L/2-dx,N);
    gauss=lambda x,d: np.exp(-0.5*(x/d)**2)/(np.sqrt(2*np.pi)*d)    
    gaussm=lambda x,d,m: np.exp(-0.5*(x/d)**m)    
    gr=gaussm(r,0.1,16);
    
    wFEM=FEM_wave_U_t(N,c=2,dx=dx,method=FEM_G)
    wFDM=FEM_wave_U_t(N,c=2,dx=dx,method=FDM_G)
    wnat=FEM_wave_U_t(N,c=2,dx=dx,method=native_G)
    xx=np.zeros((2,N),dtype=complex)
    xx[0]=gr
    t=20
    yy=wFEM(xx,t)
    yyd=wFDM(xx,t)
    yyn=wnat(xx,t)
    
    
    import matplotlib.pyplot as plt

    #%matplotlib auto
    fig=plt.figure(figsize=(18,18))
    plt.grid(True,which='major')
    plt.minorticks_on()
    plt.grid(True,which='minor',alpha=0.2)
    
    fn=lambda x: np.real(x)
    
    plt.plot(r,xx[0],label='$y_{0}$',color='#777777')
    plt.plot(r,yy[0],label='$y_{FE}$')
    plt.plot(r,yyd[0],label='$y_{FD}$')
    plt.plot(r,yyn[0],label='$y_{nat}$')
    legend = plt.legend(loc='upper right', shadow=True, fontsize=30)
    plt.show()

        