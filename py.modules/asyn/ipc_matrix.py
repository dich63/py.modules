
from .ipc_array import ipc_array
from .ipc_array import ipc_unpack
from .ipc_array import schema
from scipy.sparse import csc_matrix,coo_matrix

import numbers

class ipc_csc_matrix:
    def __init__(self,a,fakemode=False):
        
        if fakemode or isinstance(a,numbers.Number):
            self._value=a;
        else:
            fmt=a.format;
            if fmt=='coo':
                self.format='coo';
                self.col=ipc_array(a.col);
                self.row=ipc_array(a.row);
                self.data=ipc_array(a.data);
                if hasattr(a,'shape'):
                    self.shape=a.shape;
            else:
                a=a.tocsc();
                self.format=a.format;
                self.indices=ipc_array(a.indices);
                self.indptr=ipc_array(a.indptr);
                self.data=ipc_array(a.data);
                if hasattr(a,'shape'):
                    self.shape=a.shape;
            
    @property
    def value(self):
        if not hasattr(self,'_value'):
            if self.format=='coo':
                v=coo_matrix((self.data.value,(self.col.value,self.row.value)),shape=self.shape,copy=False);
                v=v.tocsc();
            else:
                v=csc_matrix((self.data.value,self.indices.value,self.indptr.value),shape=self.shape,copy=False);
            return v
        else:
            return self._value
            

ipc_matrix=ipc_csc_matrix


