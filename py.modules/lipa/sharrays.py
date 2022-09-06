import numpy as np
from scipy.sparse import csc_matrix,coo_matrix
import multiprocessing  as mp

class schema:
    def __init__(self,shape,dtype=np.float64):
        self.shape=shape;
        self.dtype=dtype;
        self.nbytes=np.empty(shape,dtype=dtype).nbytes;

class ipc_array:
    def __init__(self,a,fakemode=False):
        f=type(a)==np.ndarray;
        if fakemode==False:
            self.ra=mp.RawArray('b',a.nbytes);
            self.shape=a.shape;
            self.dtype=a.dtype
            if(f):
                self.flat[:]=a.flat;
            self._value=None;
        else:
            if f:
                self._value=a;
            else:
                self._value=np.zeros(a.shape,a.dtype);

    @property
    def flat(self):
        return np.frombuffer(self.ra,dtype=self.dtype)
    @property
    def value(self):
        if self._value==None:
            self._value=self.flat
            self._value.resize(self.shape)
        return self._value


class ipc_csc_matrix:
    def __init__(self,a,fakemode=False):
        if fakemode==False:
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
            self._value=csc_matrix((self.data.value,self.indices.value,self.indptr.value),shape=self.shape,copy=False);
        return self._value




