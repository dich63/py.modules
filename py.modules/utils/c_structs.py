

def c_update(cs,**d):
    fs=cs._fields_;    
    for f in fs:
        n=f[0];
        if n in d:            
            cs.__setattr__(n,d[n]);
    return cs



def c_2_dict(cs):
    d={};
    fs=cs._fields_;    
    for f in fs:
        d[f[0]]=cs.__getattribute__(f[0])
    return d

def c_copy(cs,**kw):
    d=c_2_dict(cs);
    d.update(kw)
    return cs.__class__(**d);

def ptr2dict(p):    
    return c_2_dict(p.contents);

'''
try:
    import numpy as np
    import ctypes
    
    class link_ptr_t(object):
        def __init__(self,p,*link):
            self._ptr,self._link=p,link
        @property
        def ptr(self):
            return self._ptr;
    
    
    
    
    def _iptr_type():
        return np.uint64 if ctypes.sizeof(ctypes.c_void_p)==8 else np.uint32
    
    def iptr(L):
        return np.array(L,dtype=_iptr_type())
    
    def ptr_array(p,*link):    
        return link_ptr_t(p.ctypes.data,p,*link)    
    
    def ptr_ptr_array(pp,*link):
        ip=iptr([p.ctypes.data for p in pp])         
        return ptr_array(ip,pp,*link);    
    
    
    
    def ptr_ptr_ptr_array(ppp,*link):       
        L=[ptr_ptr_array(pp) for pp in ppp];
        ip=iptr([l.ptr for l in L]);
        return ptr_array(ip,L,ppp,*link)
    
   
    
    
    
except:
    pass
'''