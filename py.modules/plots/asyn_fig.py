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
        self.flxy=[1,1];
        
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
        
    def set_xlim(self,vv):        
        self.ax.set_xlim(vv)
        self.flxy[0]=0;
        return self;
    
    def set_ylim(self,vv):        
        self.ax.set_ylim(vv)
        self.flxy[1]=0;
        return self;
    
    def update(self,title=None,grid=True,loc="upper right"):
        
        if not self._f:
            return self;
            
        ax,mmx,mmy=self.ax,self.mmx,self.mmy;
        ax.grid(grid)
        if self.flxy[0]:
            ax.set_xlim(mmx.mm);
            
        mmx.reset();
        
        if self.flxy[1]:
            ax.set_ylim(mmy.mm);
            
        mmy.reset();
        
        self.flxy=[1,1];
        
        p,l=[],[]
        for k,v in self.pm.items():
            p.append(v[0])
            l.append(str(v[1]))
            
        ax.legend(p,l,loc=loc);
        if not title is None: 
            ax.set_title(str(title))
        
        return self;
        


class asyn_figure(object):
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
            #print('init:',i)
            return self.init_callback(self,i);
        
        def i_callback(self,i):
            #print('initcb:',i)
            r=self.callback(self,i);
            
            for k,v in self.sp.items():
                v.update();
            #print('III....=',i);   
            if (not (self.oncomplete is None) ) and (i==self.Ne):
                
                self.oncomplete();
                
            return r;
        
        def empt(self):
            print('empt')
            pass
        
        def __call__(self,rep,callback=None,oncomplete=None,interval=200):
            self.callback=callback;
            self.oncomplete=oncomplete;
            self.Ne=int(rep-1);
            self.FAH=FuncAnimation(self.fig,self.i_callback,init_func=self.empt, frames=int(rep), interval=int(interval), blit=False, repeat=False)
            #FuncAnimation(self.fig, self.anima,init_func=self.empt, frames=int(rep), interval=200, blit=False, repeat=False)
            return self
            
        def display(self):
            plt.show()
            return self
    
            

        
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
    
    
    