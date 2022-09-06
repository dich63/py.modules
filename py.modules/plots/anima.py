import numpy as np
import cmath
from utils import *
import time
from matplotlib import pyplot as plt
from matplotlib.animation import FuncAnimation


def fake_proc(*ls,**kw):
    pass


class minmax(object):
    def __init__(self,**kw):
        self.reset(**kw);
        
    def reset(self,mi=np.inf,ma=-np.inf):
        self.mm=[mi,ma];
    def __call__(self,x):
        x=np.array(x,copy=False,dtype=np.float64).flatten();
        mi,ma=self.mm;
        xmi,xma=np.min(x),np.max(x);
        if xmi<mi:
            self.mm[0]=xmi;
        if xma>ma:
            self.mm[1]=xma;
        return self;
    
        
class multiplot(object):
    def __init__(self,fig,*ls,**kw):
        self.fig=fig 
        self._f=0;
        
    def reset(self,*ls,**kw):        
        self.ax=self.fig.add_subplot(*ls,**kw);
        self.pm={};
        self.mmx=minmax();
        self.mmy=minmax();
        
    def set_plot(self,key,data,linewidth=0.5,color='k',legend=None):
        pm=self.pm
        if not key in pm:
            pm[key]=pls=[self.ax.plot([], [])[0],key];
        else:
            pls=pm[key];    
                        
        if not legend is None:
            pls[1]=legend;
            
        pl,t=pls;    
            
        pl.set_data(*data);
        pl.set_color(color);
        pl.set_linewidth(linewidth);
        self.mmx(data[0]);
        self.mmy(data[1]);
        self._f=1;
        return self;
    
    def __call__(self,key,data,linewidth=0.5,color='k',legend=None):
        return self.set_plot(key,data,linewidth=linewidth,color=color,legend=legend);
    
    def set_title(self,*ls):
        self.ax.set_title(sprintf(*ls));
        return self;
        
    
    def update(self,title=None,grid=True,loc="upper right"):
        
        if not self._f:
            return self;
            
        ax,mmx,mmy=self.ax,self.mmx,self.mmy;
        ax.grid(grid)
        
        ax.set_xlim(mmx.mm);
        mmx.reset();
        ax.set_ylim(mmy.mm);
        mmy.reset();
        
        p,l=[],[]
        for k,v in self.pm.items():
            p.append(v[0])
            l.append(str(v[1]))
            
        ax.legend(p,l,loc=loc);
        if not title is None: 
            ax.set_title(str(title))
        
        return self;
        


class asio_figure(object):
        def __init__(self,*ls,**kw):
            self.fig=plt.figure(*ls,**kw);
            self.sp={};
            self.init_callback=fake_proc;
            self.callback=fake_proc;
            
        def __getitem__(self,tag):
            sp=self.sp;
            if not tag in sp:
                sp[tag]=mp=multiplot(self.fig);
            else:
                mp=sp[tag];
            return mp;
        def i_init_callback(self,i):
            print('init:',i)
            return self.init_callback(self,i);
        
        def i_callback(self,i):
            print('initcb:',i)
            r=self.callback(self,i);
            
            for k,v in self.sp.items():
                v.update();
                
            return r;
        
        def empt(self):
            print('empt')
            pass
        
        def __call__(self,rep,callback=None,interval=200):
            self.callback=callback;
            
            def cbf(i):
                print('cbf i:',i)
                return self.i_callback(i);
            self.anima=cbf
            self.FAH=FuncAnimation(self.fig,cbf,init_func=self.empt, frames=int(rep), interval=int(interval), blit=False, repeat=False)
            #FuncAnimation(self.fig, self.anima,init_func=self.empt, frames=int(rep), interval=200, blit=False, repeat=False)
            return self
            
        def display(self):
            plt.show()
            return self
    
            
class Amt(object):
    def __init__(self):
        self.fig = plt.figure(figsize=(8,5));
    
    def empt(self):
        print('empt')
        pass
    
    def empt2(self,i):
        print('empt2',i)
        pass
    
    def __call__(self,NZ):
        self.NZ=NZ;
        self.nz=0;
        #        self.anim = FuncAnimation(self.fig, self.anima,init_func=self.empt, frames=int(NZ), interval=200, blit=False, repeat=False)
        #
        self.anim = FuncAnimation(self.fig,self.empt2,init_func=self.empt, frames=int(NZ), interval=200, blit=False, repeat=False)
        return self
    
    def display(self):
        plt.show()
        return self
        
    def w_call__(self,NZ):
        return self.StartAnim(NZ)
    
        

class Animation(object):
    numlines=0;
    Lines=list();
    def __init__(self,xxr,tt,hz,fun):
        self.fig = plt.figure(figsize=(8,5));
        
        self.animclosure(xxr,tt,hz,fun);
        '''
        if (subplots>1):
            self.ax=list(range(subplots));
            for i in range(0,subplots):
                self.ax[i] = self.fig.add_subplot(subplots,1,i+1,xlim=(-25, 25), ylim=(1e-13,100), yscale='log');
        '''
               
    def initiate(self,subplots=1 ):
        del self.fig;
        #del self.ax
        self.fig = plt.figure(1,figsize=(8,5));
        if (subplots>1):
            self.ax=list(range(subplots));
            for i in range(0,subplots):
                if i>0:
                    self.ax[i] = self.fig.add_subplot(subplots,1,i+1,xlim=(-25, 25), ylim=(1e-13,100), yscale='log');   
                else:
                    self.ax[i] = self.fig.add_subplot(subplots,1,i+1,xlim=(-25, 25),ylim=(-4,4));#,autoscaley_on=True)#, );    
               
    def set_data(self,subplot=1, yscale='linear',lw=1):
        if (yscale=='linear'):
            line=self.ax[subplot-1].plot([], [], linewidth=lw)[0]
            self.Lines.append(line);
        elif(yscale == 'log'):
            #line=self.ax[subplot-1].semilogy([], [], linewidth=lw)[0]
            line=self.ax[subplot-1].plot([], [], linewidth=lw)[0]
            self.Lines.append(line);
        self.numlines=self.numlines+1;
        
    def animclosure(self,xxr,tt,hz,fun):
    
        self.xxr=xxr;
        self.tt=tt;
        self.hz=hz;
        self.fun=fun;
        NZ,N=xxr.shape
        self.nr=np.linalg.norm(xxr[0,:])*(np.max(tt)-np.min(tt))/ N;
        self.z=0
        errmax=-1e100;
        xr=xxr[0,:];
        
        self.Lines=[];
        #fig = plt.figure(1,figsize=(8,5))
        self.initiate(2);
        for i in range(4):
            self.set_data(1,'linear', (lambda lw: 0.5 if(lw==0) else 1)(i));
        for i in range(2):
            self.set_data(2,'log', (lambda lw: 0.5 if(lw==0) else 1)(i));
        #ax = plt.axes(xlim=(np.min(tt), np.max(tt)), ylim=(-2,2))
        #ax = fig.add_subplot(211,xlim=(np.min(tt), np.max(tt)), ylim=(-4,4))
       # line, = ax.plot([], [], linewidth=0.5)
       # line2, = ax.plot([], [], linewidth=1)
       # line11,=ax.plot([], [], linewidth=1)
       # line21,=ax.plot([], [], linewidth=1)

       # ax2=fig.add_subplot(212,xlim=(np.min(tt), np.max(tt)), ylim=(1e-13,100), yscale='log');
       # line3, = ax2.semilogy([], [], linewidth=0.5)
       # line4, = ax2.semilogy([], [], linewidth=1)
        #xr=np.arange(np.size(t)*int(REP), dtype=np.clongdouble).reshape(int(REP),np.size(t));
        #xr=np.arange(np.size(tt)*int(NZ), dtype=np.clongdouble).reshape(int(NZ),np.size(tt));
       # self.xe=np.arange(np.size(tt)*int(NZ), dtype=np.clongdouble).reshape(int(NZ),np.size(tt));
        #self.axe=np.arange(np.size(tt)*int(NZ), dtype=np.clongdouble).reshape(int(NZ),np.size(tt));
        #for n in range(2,NZ):
        print(self.hz)
        def animate1(i):
            
            print('i,nz1',i,self.nz)
            
            self.z=self.z+self.hz;
            [self.xxr,self.xe,elapsed]=fun(self.tt,self.z);
            tcpu,niter=elapsed;
            
            if i!=self.nz:
                print('i,nz',i,self.nz)
            
            n,N=self.nz,self.NZ
            self.ax[0].set_title(sprintf("[i=%d %d:%d] elapsed: $\\tau_{cpu}$=%3.3f sec iter=%d, %d",i,n,N,tcpu,niter,i))
            #
            #print('i,nz1',i,self.nz)
            self.nz=self.nz+1;
            #print('i,nz2',i,self.nz)
            self.axe=np.abs(self.xe);
            self.Lines[0].set_data(self.tt, self.axe);
            self.Lines[0].set_label('Psi')


            self.Lines[1].set_data(self.tt, np.abs(self.xxr));
            self.Lines[2].set_data(self.tt, np.real(self.xe));
            self.Lines[3].set_data(self.tt, np.real(self.xxr));

            self.Lines[4].set_data(self.tt,(np.abs(np.abs(self.xe)-np.abs(self.xxr))/self.nr)**2);
            self.Lines[5].set_data(self.tt,np.abs(((self.xe-self.xxr)/self.nr))**2);
         
            #self.em=20*np.log10(np.max(np.abs(self.xe[i+1]-self.xxr[i+1])/self.nr));
            #if errmax<self.em:
            #    errmax=self.em;
            #self.ax[1].set_title(sprintf("$\\epsilon=%3.2f [\\epsilon_{max}:%3.2f][dB]$  ",em,errmax))
            self.ax[1].grid(True)
            self.ax[0].grid(True)
            #self.legend = plt.legend()
            self.ax[0].legend(self.Lines[0:4],['$|\\Psi_{an}|$','$|\\Psi_{num}|$','Re$\\Psi_{an}$','Re$\\Psi_{num}$'],loc="upper right")
            self.ax[1].legend(self.Lines[4:6],['$\\delta|\\Psi|$','$\\delta\\Psi$'],loc="upper right")

            #return self.Lines[0],self.Lines[1],self.Lines[2],self.Lines[3],self.Lines[4],self.Lines[5]
        
        self.anima=animate1;
        return animate1
    
    def empt(self):
        print('empt')
        pass
    
    def empt2(self,i):
        print('empt2',i)
        pass
    
    def StartAnim(self,NZ):
        self.NZ=NZ;
        self.nz=0;
        #
        self.anim = FuncAnimation(self.fig, self.anima,init_func=self.empt, frames=int(NZ), interval=200, blit=False, repeat=False)
        #self.anim = FuncAnimation(self.fig,self.empt2,init_func=self.empt, frames=int(NZ), interval=200, blit=False, repeat=False)
        return self
    
    def display(self):
        plt.show()
        return self
        
    def __call__(self,NZ):
        return self.StartAnim(NZ)
        
if __name__=='__main__':
    
    fig = plt.figure(figsize=(8,5))
    #axx=fig.add_subplot(211)
    mp=multiplot(fig)
    mp.reset(211);
    tt1=np.linspace(0,5,20);
    tt2=np.linspace(-5,3,20);
    mp.set_plot('a122a',[tt1,2+np.sin(tt1)],color='r')
    mp.set_plot('a22a',[tt2,-np.cos(tt2)],color='g')
    mp.update()
    
    #axx2=fig.add_subplot(212, yscale='log')
    mp2=multiplot(fig)
    mp2.reset(212, yscale='log')
    tt1=np.linspace(0,5,20);
    tt2=np.linspace(-5,3,20);
    mp2.set_plot(1,[tt1,2+np.sin(tt1)],color='r',legend='$|\\Psi_{an}|$')
    mp2.set_plot(2,[tt2,-np.cos(tt2)],color='g')
    mp2.update('AAAAAAAAAAAAA')
    
    
    