import ctypes
import numpy as np
import multiprocessing  as mp
#from multiprocessing.reduction import ForkingPickler

from scipy.sparse import csc_matrix,coo_matrix
#from .shmp import sh_array
"""
def sh_array(count):    
          sh=count>>32;
          sl=count and 0xFFFFFFFF;
          np.dtype(np.complex).itemsize
          ct=c_byte*sh
          w=mp.heap.BufferWrapper(2**31)
          buf = w.create_memoryview()
            o=ct.from_buffer(buf)            
            w=mp.heap.BufferWrapper(2**31)
            buf = w.create_memoryview()
            ct=c_byte*2*g
            o=ct.from_buffer(buf)
"""    

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
    def __init__(self,shape,dtype=np.float64):
        self.shape=shape;
        self.dtype=dtype;
        self.nbytes=np.empty(shape,dtype=dtype).nbytes;

class ipc_array:
    def __init__(self,a,fakemode=False):
        f=type(a)==np.ndarray;
        if fakemode==False:
            #self.flat=a.flat if type(a)==ipc_array else mp.RawArray('b',a.nbytes);
            dtype=np.dtype(a.dtype)
            self.count=np.int64(a.nbytes/dtype.itemsize);
            #self._ra=mp.RawArray('b',a.nbytes);
            self._ra=create_raw_buffer(a.nbytes);
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
            self._value=self.flat
            self._value.resize(self.shape)
        return self._value


class ipc_csc_matrix:
    def __init__(self,a,fakemode=False):
        if fakemode==False:
            fmt=a.format;
            if fmt=='coo':
                self.format='coo';
                self.col=ipc_array(a.col);
                self.row=ipc_array(a.row);
                self.data=ipc_array(a.data);
            else:
                a=a.tocsc();
                self.format=a.format;
                self.indices=ipc_array(a.indices);
                self.indptr=ipc_array(a.indptr);
                self.data=ipc_array(a.data);
                self.shape=a.shape;
                self._value=None;
        else:
            self._value=a;

    @property
    def value(self):
        if self._value==None:
            if self.format=='coo':
                v=coo_matrix((self.data,(self.col,self.row)),shape=self.shape,copy=False);
                self._value=v.tocsc();
            else:
                self._value=csc_matrix((self.data.value,self.indices.value,self.indptr.value),shape=self.shape,copy=False);
        return self._value
    def tocsc(self):
        return self.value.tocsc();




