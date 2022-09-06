import ctypes
import numpy as np
import multiprocessing  as mp

def create_raw_buffer(size):
    size=np.int64(size);
    if size<0x080000000:
        return mp.RawArray('b',int(size));
    else:
        sp=size>>3;
        r=size & 0x07
        if r:
            sp+=1;
        return mp.RawArray('d',int(sp));

class schema:
    def __init__(self,shape,dtype=np.float64,capacity=-1):
        self.shape=shape;
        self.dtype=dtype;
        self.nbytes=cb=np.dtype(dtype).itemsize*np.prod(shape);
        self.capacity=cb if capacity<cb else capacity
        #self.nbytes=np.empty((1,),dtype=dtype).nbytes*np.prod(shape);
        #self.nbytes=np.empty(shape,dtype=dtype).nbytes;

def ipc_unpack(a):
    if (type(a)==tuple) or (type(a)==list):
        a=[i.value for i in a]
    else:
        a=a.value;
    return a;


class ipc_array:
    
    @staticmethod
    def cast_to_ndarray(a):
        if type(a)==ipc_array:
            return a.value
        elif type(a)==np.ndarray:
            return a;
        else:
            raise TypeError('Can''t convert to ndarray from '+str(type(a)))
        
    def __init__(self,a,fakemode=False):
        f=type(a)==np.ndarray;
        if fakemode==False:
            #self.flat=a.flat if type(a)==ipc_array else mp.RawArray('b',a.nbytes);
            dtype=np.dtype(a.dtype)
            self.count=np.int64(a.nbytes/dtype.itemsize);
            #self._ra=mp.RawArray('b',a.nbytes);
            cb=a.nbytes if f else a.capacity
            self._ra=create_raw_buffer(cb);
            self.shape=a.shape;
            self.dtype=dtype;

            if(f):
                self.flat[:]=a.flat;

            #self._value=None;
        else:
            if f:
                self._value=a;
            else:
                self._value=np.zeros(a.shape,a.dtype);

    @property
    def flat(self,offset=0):
        return np.frombuffer(self._ra,count=self.count,dtype=self.dtype,offset=offset);
    @property
    def value(self):
        if not hasattr(self,'_value'):
            v=self.flat;
            v.resize(self.shape);
            return v
        else:
            #self._value.resize(self.shape)
            return self._value
    @value.setter
    def value(self,v):
        self.value[:]=v


class ipc_cache(object):
    def __init__(self):
        self._links=[];
        
    def link_array(self,a):
        ia=ipc_array(a);
        self._links.append(ia);
        return ia;
    
    def __call__(self,o):
        t=type(o)
        if t is ipc_array:
            self._links.append(o);
        elif t is dict or  hasattr(o,'__dict__'):
            d= o   if  t is dict else  o.__dict__;
            for k in d:
                self(d[k]);
        elif t in (list,tuple):
            for v in o:
                self(v);
            
            
        return o;
    
if __name__=='__main__':
    
    
    from asyn.ipc_matrix import *
    from scipy.sparse import csc_matrix    
    A = csc_matrix([[1., 0., 0.], [5., 0., 2.], [0., -1., 0.]], dtype=float)
    ia=ipc_matrix(A)
    c=ipc_cache()
    c([ia,ia,[ia]])
    
