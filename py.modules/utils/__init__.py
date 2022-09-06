#
import time
import sys,os
import copy

from .printf import printf,sprintf


__cb_dict__={}

def set_cb_log(key,fun=None):
    __cb_dict__[key]=fun;
    if fun is None:
        del __cb_dict__[key];
    
def call_cb_log(key,*lp):
    try:
        return __cb_dict__[key](*lp);
    except:
        return None;
           

class nul_t(object):
    def __init__(this):
        pass
    @property 
    def v(this):
        return None
    
    @v.setter
    def v(this,v):
        pass
    
    def __getattr__(self,name):
        return nul_t();
    
    def __setattr__(self,name,value):
        pass
    def __call__(self,*kw,**dw):
        return nul_t();
    
    def __getitem__(self,*kw):        
        return nul_t();
    
    def __setitem__(self,*kw):         
        pass
    def __bool__(self):
        return False
    
    
nul=nul_t()


def check_error(f,msg='Error'):
        if not f :
            raise Exception(msg);

def safe_call(f,kw=[]):
    try:
        f(*kw)
    except Exception:
        pass
    

def isempty(x):
    if hasattr(x,'__len__'):
        return len(x)==0;
    return  not x;
    
def val_def(v,d):
    return d if isempty(v) else v;

def opt_def(o,d={},fsametype=False):
    

    
    def cast_dict(x):
        if hasattr(x,'__dict__'):
            x=x.__dict__;
        return copy.copy(x);
    
    
    t_o=type(o);       
    
    o,d=[ cast_dict(val_def(x,{})) for x in (o,d) ];
        
    #o=copy.deepcopy(o);
    #d=copy.deepcopy(d);
    d.update(o);    
        
    return t_o(d) if fsametype else d;

def setifNone(o,newo):
    if o is None:
        o=newo;
    return o;

class _to(object):
    def __init(self):
        self.reset();
        
    def reset(self):
        self.t=time.perf_counter();
        
    def sec(self):
        return time.perf_counter()-self.t;

__t__=_to();        
        
def __getTimer(t=None):
    global __t__
    return __t__ if t is None else t;

def __resetTimer(t=None):
    global __t__
    if type(t)==_to:
        return t;
    else:  
        return __t__ if t is None else _to();    

def tic(s=None,t=None):       
    tm=__resetTimer(t);    
    if not s is None:
        printf(s)
    tm.reset();   
    return t if t is None else tm;

def toc(s=None,t=None):
    
    t=__getTimer(t).sec();
    if not s is None:
        printf('%s %g sec\n',s,t); 
    return t;
    pass

def mtoc(s=None,t=None):
    
    t=__getTimer(t).sec()*1000;
    if not s is None:
        printf('%s %3.3g Ms\n',s,t); 
    return t;
    pass


def readln():
    s='';
    while True:
        c=sys.stdin.read(1);
        #sys.stdout.write(c);
        if (c=='') or (ord(c)==10)  :
            break
        s+=c;
    return s; 
    