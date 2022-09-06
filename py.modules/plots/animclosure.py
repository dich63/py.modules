# -*- coding: utf-8 -*-

from utils import *
import cmath
import numpy as np
from skimage.transform import resize
from matplotlib import pyplot as plt
from matplotlib.animation import FuncAnimation

import k3d

def animclosure(fig,xxr,tt,hz,fun):
    
    x0=xxr[0,:];
    NZ,N=xxr.shape
    
    nr=np.linalg.norm(x0)*(np.max(tt)-np.min(tt))/N;
    z=0
    errmax=-1e100;
    

    #ax = plt.axes(xlim=(np.min(tt), np.max(tt)), ylim=(-2,2))
    ax = fig.add_subplot(211,xlim=(np.min(tt), np.max(tt)), ylim=(-4,4))
    
    ax.legend(loc=4) 
    
    
    
    ax2=fig.add_subplot(212,xlim=(np.min(tt), np.max(tt)), ylim=(1e-13,100), yscale='log');
    
    
    line, = ax.plot([], [], linewidth=0.5)
    line2, = ax.plot([], [], linewidth=1)
    line11,=ax.plot([], [], linewidth=1)
    line21,=ax.plot([], [], linewidth=1)
    
    
    
    line3, = ax2.semilogy([], [], linewidth=0.5)
    line4, = ax2.semilogy([], [], linewidth=1)
    
    #xr=np.arange(np.size(t)*int(REP), dtype=np.clongdouble).reshape(int(REP),np.size(t));
    #xr=np.arange(np.size(tt)*int(NZ), dtype=np.clongdouble).reshape(int(NZ),np.size(tt));
    xe=np.arange(np.size(tt)*int(NZ), dtype=np.clongdouble).reshape(int(NZ),np.size(tt));
    axe=np.arange(np.size(tt)*int(NZ), dtype=np.clongdouble).reshape(int(NZ),np.size(tt));
    #for n in range(2,NZ):



    def animate1(i):
        
        #te=nls.elapsed[0];
        #ne=nls.elapsed[1];
        #xxr[n]=xr[i+1];
        global hz;
        nonlocal z,nr,errmax,ax,ax2,xe,axe,fig,xxr;
        nonlocal line,line2,line11,line21,line3,line4;
        
    
        z=z+hz;
        #xe[i+1]=np.conj(brizerKM(tt,z,A));  
        #xxr[i+1]=nls(rep=int(mm),pp=1);
        #xe[i+1]=(brizerKM(tt,z,A));
        xxr[i+1],xe[i+1],elapsed=fun(tt,z);
        tcpu,niter=elapsed;
        ax.set_title(sprintf("elapsed: $\\tau_{cpu}$=%3.3f sec iter=%d",tcpu,niter))
        
        axe[i+1]=np.abs(xe[i+1]);
        line.set_data(tt, axe[i+1]);
        line.set_label('Psi')
        
        
        line2.set_data(tt, np.abs(xxr[i+1]));
        line11.set_data(tt, np.real(xe[i+1]));
        line21.set_data(tt, np.real(xxr[i+1]));

        line3.set_data(tt,(np.abs(np.abs(xe[i+1])-np.abs(xxr[i+1]))/nr)**2);
        line4.set_data(tt,np.abs(((xe[i+1]-xxr[i+1])/nr))**2);

        #set(0, 'CurrentFigure', fh);
        #subplot(2,1,1)

        #plot(tt,abs(xe),'b',tt,abs(xr),'m',tt,gfun(xe),'k',tt,gfun(xr),'g');grid minor;legend('|\Psi_{an}|','|\Psi_{num}|','Re\Psi_{an}','Re\Psi_{num}')
        #title(sprintf('[%d] z=%f elapsed=[sec:%3.3f iter:%d ]',n,z,te,ne));
        # subplot(2,1,2)
        #semilogy(tt,(abs(abs(xe)-abs(xr))/nr).^2,tt,(abs(xe-xr)/nr).^2);grid minor;legend('|\delta\rho|','|\delta\Psi|')
        #axis([-LT,LT,1e-14,5]);
        em=20*np.log10(np.max(np.abs(xe[i+1]-xxr[i+1])/nr));
        if errmax<em:
            errmax=em;
        ax2.set_title(sprintf("$\\epsilon=%3.2f [\\epsilon_{max}:%3.2f][dB]$  ",em,errmax))
        ax2.grid(True)
        ax.grid(True)
        legend = plt.legend()
        return [line,line2,line11,line21,line3,line4]+[legend]
    
    
    #
    
    #animate1(0);
    
    
    return animate1

