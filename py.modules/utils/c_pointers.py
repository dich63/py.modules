from .c_structs import *


import numpy as np
from ctypes import *

#c_void_p.from_buffer(n).value



def c_is_null(pv):
    return c_void_p.from_buffer(n).value==None
    #return (cast(pv, c_void_p).value == None)

class link_ptr_t(object):
    def __init__(self,p=c_void_p(),*link):
        self._ptr,self._link=p,link
    @property
    def ptr(self):
        return self._ptr;
    
    @property
    def iptr(self):
        return addressof(self._ptr.contents);
    @property
    def value(self):
        return self._ptr.contents;
    @property
    def ibyref(self):
        return addressof(self._ptr);
    @property
    def byref(self):
        return byref(self._ptr);
    
smart_ptr_t=link_ptr_t

class cache_t(object):
    def __init__(self):
        self._list=[];
    def __call__(self,o):
        self._list+=[o];
        return o;

        


def ptr_sizeof():
    return sizeof(c_void_p)

def iptr_type():
    return np.int64 if sizeof(c_void_p)==8 else np.int32

def iptr(L):
    return np.array(L,dtype=iptr_type())

def ptr_array(p,*link):    
    return link_ptr_t(p.ctypes.data,p,*link)    

def ptr_ptr_array(pp,*link):
    ip=iptr([p.ctypes.data for p in pp])         
    return ptr_array(ip,pp,*link);    

    '''
    L=[p.ctypes.data for p in pp]         
    ip=np.array(L,dtype=_iptr_type());        
    return link_ptr_t(ip.ctypes.data,ip,pp);         
    '''


def ptr_ptr_ptr_array(ppp,*link):
    L=[ptr_ptr_array(pp) for pp in ppp];
    ip=iptr([l.ptr for l in L]);
    return ptr_array(ip,L,ppp,*link)


def _foreachNone(fun,xx):
    if not fun is None:
            for x in xx:
                fun(x);
    return xx;

class handle_array_t(object):
    def __init__(self,xx=None,nbatch=0,release=None,addref=None ):
        
        self._release,self._addref=release,addref;
        if xx is None:
            self._handles=np.zeros(nbatch,dtype=iptr_type());
        else:
            self._handles=iptr([x.ctypes.data for x in xx]);
        
                
    def __del__(self):
        
        _foreachNone(self._release,self._handles)
        
            
    def __len__(self):
        return len(self._handles)       
        
    @property
    def values(self):
        return self._handles        
    @property
    def addresses(self):
        return self._handles        
    @property
    def ptr(self):
        return self._handles.ctypes.data
    
    def __getitem__(self,i):        
        return int(self._handles[i]);
    
    def subset(self,mask):
        xx=self._handles[mask];
        addref=self._addref;
        
        hm=handle_array_t(release=self._release,addref=addref)        
        
        hm._handles=_foreachNone(addref,xx)
        
        return hm;
    
    def __add__(self,o):       
        release,addref=self._release,self._addref;
        if o._release==release and o._addref==addref:
            xx=np.concatenate((self._handles,o._handles));         
            hm=handle_array_t(release=release,addref=addref);              
            hm._handles=_foreachNone(addref,xx);
            return hm;
        
        
        
            
        



    '''
    pL=np.array([l.ptr for l in L],dtype=_iptr_type())
    return link_ptr_t(pL.ctypes.data,pL,L,ppp,*link); 
    '''

class rhs_t(Structure):
    _pack_=8
    _fields_ = [    
    ("xx",c_void_p),  
    ("nrhs",c_long),        
    ]


p_rhs_t  = POINTER(rhs_t)


class xx_buffer_t(object):   
    
    def __init__(self,n=-1,nbatch=-1,xx=None,nrhs=1,dtype=np.complex128):
        
        #xxl=[np.ascontiguousarray(np.zeros((nrsh,n),dtype=dtype)) for k in range(nbatch)];
        if xx is None:
            xx=np.empty(nbatch,dtype=object);     
            for k in range(nbatch):
                xx[k]=np.ascontiguousarray(np.zeros((nrhs,n),dtype=dtype))       
        else:
            nm=np.prod(xx[0].shape);
            if n>0:                
                if nm%n:
                    raise Exception('bad xx data');
                nrhs=nm//n
            else:
                if nm%nrhs:
                    raise Exception('bad xx data');
                n= nm//nrhs
                             
        self._xx=xx;        
        self._ppxx=ppxx=ptr_ptr_array(xx);        
        self._rhs=rhs_t(xx=ppxx.ptr,nrhs=nrhs)
        
    @property
    def ptr(self):
        return pointer(self._rhs);
    
    @property
    def xx(self):
        return self._xx;
    
    def __getitem__(self,i):        
        return self._xx[i];
    
    def __setitem__(self,i,value):
        self._xx[i][:]=value;
        
    def __iter__(self):    
        return iter(self._xx)
    
    def __len__(self):
        return len(self._xx)

