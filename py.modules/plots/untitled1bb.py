# -*- coding: utf-8 -*-
#from utils import *
from NLS.brizers import *
from NLS.nls_adapter import *
#import cmath
#import numpy as np
#from skimage.transform import resize

#from matplotlib import pyplot as plt
#from matplotlib.animation import FuncAnimation
from plots.anima import Animation, asio_figure,Amt
from plots.asyn_fig import asyn_figure
#
from plots.k3dsurf import surf,cmps,resize
from matplotlib.animation import FuncAnimation
#from plots.animclosure import *


#import k3d



A=1-1/8
A=1
N=16*1024
NZ=1000;
NZ=20;

mm=256+256;
mm=2;
mm=32
mm=8

omagnus=1
#mm=256;
#omagnus=0
flags=-1;

flags=0;
alpha=0*1e-8;

sname='nls-exp-dec';
sname='nls';
#sname='ssf';omagnus=0
[pT,pZ]=brizerKM_periods(A);
LT=5*abs(pT);
LZ=8/11*2*2.5*abs(pZ);
NvT=200;
NvZ=300;
tLv=np.linspace(-LT,LT,NvT);
zLv=np.linspace(-LZ*0,LZ,NvZ);
[TTv,ZZv]=np.meshgrid(tLv,zLv);
#figure
#tic




FFv=brizerKM(TTv,ZZv,A);

tt=np.linspace(-LT,LT,N+1);
tt=tt[0:N]
g=1;
[nls,dzu,h]=nls_rescale_to(tt,sname,g);
hz=LZ/NZ; 
dz=hz* dzu/mm

nls.reset(dt=dz,nm=[4,4],omagnus=omagnus,alpha=alpha,w=1/3,flags=flags);




x0=brizerKM(tt,0,A);

nr=np.linalg.norm(x0)*(np.max(tt)-np.min(tt))/N;

z=np.array([0.0,nr,-np.inf]);

xxr=np.zeros((NZ,N), dtype=complex);




def callback(af,i,A,tt,nls,z,hz,mm):           
    
    
    
    xe=brizerKM(tt,z[0],A);
    if i==0:
        nls.x=xr=xe;
    else:
        xr=nls(rep=int(mm),pp=1);
        
    tcpu,niter=nls.elapsed;
    
    xxr[i,:]=xr;
    
    z[0]+=hz;
    
    
    
    
    
    
    axr=np.abs(xr);
    axe=np.abs(xe);
    
    af[1](1,(tt,axr),color='#33f',legend='$|\\Psi_{num}|$',linewidth=1.5)
    af[1](2,(tt,axe),color='#007',legend='$|\\Psi_{an}|$',linewidth=0.5)
    
    af[1](3,(tt,np.real(xr)),color='#3f3',legend='$Re\\Psi_{num}$',linewidth=1.5)
    af[1](4,(tt,np.real(xe)),color='#070',legend='$Re\\Psi_{an}$',linewidth=0.5)
    
    af[1].set_title("[%d] elapsed: $\\tau_{cpu}$=%3.3f sec iter=%d, ",i+1,tcpu,niter)
    
    nr=z[1];
    
    err_r=(np.abs(axr-axe)/nr)**2;
    err_f=(np.abs(xr-xe)/nr)**2;
    
    af[2](1,(tt,err_r),color='#33f',legend='$\\delta\\rho$')
    af[2](2,(tt,err_f),color='#3f3',legend='$\\delta\\Psi$')
    
    em=10*np.log10(np.max(err_f));
    if z[2]<em:
        z[2]=em;
    
    af[2].set_title("$\\epsilon=%3.2f [\\epsilon_{max}:%3.2f][dB]$  ",em,z[2]);
    
    af[2].set_ylim((1e-12,100))
    



af=asyn_figure(figsize=(8,5))
af[1].reset(211);
af[2].reset(212, yscale='log');
print('Start....')
af(NZ,lambda af,i: callback(af,i,A,tt,nls,z,hz,mm),lambda : print("+++++End+++!")).display()

    
'''        
def fa(i):
    print('fa:',i)

        

        




fun = lambda tt,z: [brizerKM(tt,z,A),nls(rep=int(mm),pp=1),nls.elapsed] 
#fun = lambda tt,z: [brizerKM(tt,z,A),brizerKM(tt,z,A),nls.elapsed] 



xr=x0;
xe=x0;

a=Amt()
a(187).display();




errmax=-1e100;
a=Animation(xxr,tt,hz,fun)
a(187).display();

b=Animation();
animate1=b.animclosure(xxr,tt,hz,fun)

#anim = FuncAnimation(a.fig, animate1, frames=int(NZ), interval=200, blit=True, repeat=False)
#animate1(100);
b(200).show();

'''
#plt.show()
