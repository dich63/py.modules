# -*- coding: utf-8 -*-
#from utils import *
from NLS.brizers import *
from NLS.nls_adapter import *
#import cmath
#import numpy as np
#from skimage.transform import resize

#from matplotlib import pyplot as plt
#from matplotlib.animation import FuncAnimation
from plots.anima import Animation
#
from plots.k3dsurf import surf,cmps,resize
#from plots.animclosure import *


#import k3d



A=1-1/8
A=1
N=16*1024
NZ=1000;
NZ=200;

mm=256+256;
mm=2;
mm=32
mm=8

omagnus=0
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

z=0;
x0=brizerKM(tt,z,A);


fun = lambda tt,z: [brizerKM(tt,z,A),nls(rep=int(mm),pp=1),nls.elapsed] 
#fun = lambda tt,z: [brizerKM(tt,z,A),brizerKM(tt,z,A),nls.elapsed] 


nls.x=x0;
xr=x0;
xe=x0;

xxr=np.zeros((NZ,N), dtype=complex);
xxr[0,:]=xr;



errmax=-1e100;
a=Animation(xxr,tt,hz,fun)
a(187).display();
'''
b=Animation();
animate1=b.animclosure(xxr,tt,hz,fun)
'
#anim = FuncAnimation(a.fig, animate1, frames=int(NZ), interval=200, blit=True, repeat=False)
#animate1(100);
b(200).show();
'''
#plt.show()
