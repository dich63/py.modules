# @@DICH 2020

import copy
from jsobj import *

def _rem_none(d):
    r={}
    for k,v in d.items():
        if not v is None:
            r[k]=v;
    return r;

class js_getter_setter(object):
    def __init__(self, getter=lambda k: 'None' ,setter=lambda k,v: 'sNone',caller=lambda *d : '??'):
        self.__dict__['__getter__']=getter;
        self.__dict__['__setter__']=setter;
        self.__dict__['__caller__']=caller;
        
    def __call__(self,*d):
        return self.__caller__(*d);
    
    def __getattr__(self,name):
        if name.startswith('__'):
            #print('Name:',name)
            return self.__dict__[name];
        else:     
            return self.__getter__(name);
    
    def __getitem__(self,name):     
        return self.__getattr__(name);    
    
    
    def __setattr__(self,name,value):
        #print('::',name,value);#,self.__dict__)                
        return self.__setter__(name,value);         
   
    def __setitem__(self,name,value):         
        return self.__setattr__(name,value);
    


class _jsbase(object):
    def _def(self,name,val=None):
        return self.__dict__.get(str(name),val);
    
    def __setattr__(self,name,value):
        #print('::',name,value,self.__dict__)                
        if (value is None):
            if (name in self.__dict__) :
                del self.__dict__[name];
        else:
            self.__dict__[name]=value;
    def __getstate__(self):
            return self.__dict__;
    def __setstate__(self,state):
            self.__dict__.update(state);
            

class jsobject(_jsbase):
    def __init__(self, entries={}):
        self.__dict__.update(_rem_none(get_dict(entries)))
        
    def __repr__(self):
        d=self.__dict__
        if len(d):
            s='jso:{\n'
            for k,v in d.items():
                s+=str(k)+':'+repr(v)+'\n'
            
            s+='}'            
        else:
            s='jso:{}';            
        return s;
    def __str__(self):
        return 'jso:'+str(self.__dict__)

            
    def __getattribute__(self,attr):
        if attr.startswith('__') or attr=='_def':                            
            return super(jsobject,self).__getattribute__(attr)
        else:
            return self.__dict__.get(attr,None)
    def __getitem__(self,name):        
        return self.__getattribute__(str(name))
    
    def __setitem__(self,name,value):         
        super(jsobject,self).__setattr__(str(name),value)
        #print('post::',name,value,self.__dict__)       
            
    def __iter__(self):
        #print('iter...')
        return iter(self.__dict__)
    
    def __def(self,name,value):        
        return self.__dict__.get(str(name),value);
       
        
   
    
    def __iadd__(self,o):        
        self.__dict__.update(_rem_none(o.__dict__ if type(o)==jsobject else o))
        return self;
    def __add__(self,o):        
        d=jsobject()
        d+=self
        d+=o
        return d;



def get_dict(o):
    return o.__dict__ if type(o)==jsobject else o;
    

def ext_def(o,dfl={},fdeep=False):
    dfl=copy.deepcopy(dfl);
    if fdeep:
        o=copy.deepcopy(o);
    dfl.update(o)
    return dfl;
    
    
def as_dict(o):
    if type(o) in (dict,):
        return o;
    return o.__dict__        
    
def to_dict(o):
    d={};
    d.update(as_dict(o))
    return d
    
def to_obj(d):
    return jsobject(d);

def to_jso(o):
    return o if type(o)==jsobject else jsobject(o);

def to_obj_r(o):
    t=type(o)
    if t is dict:  
        o=jsobject(o)         
        for k in o:
            o[k]=to_obj_r(o[k]);
            
    elif t in (list,tuple):  
        l=len(o);
        for k in range(l):
            o[k]=to_obj_r(o[k]);
        
    return o;

dict2jso=to_obj_r
    
def jso2dict(o):
    t=type(o)
    if t is jsobject:  
        o=o.__dict__;         
        for k in o:
            o[k]=jso2dict(o[k]);
            
    elif t in (list,tuple):  
        l=len(o);
        for k in range(l):
            o[k]=jso2dict(o[k]);
        
    return o;

def arg2jso(**kws):
    return jsobject(kws);

    

jsc=arg2jso;
jso=jsobject;

