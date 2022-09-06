# -*- coding: utf-8 -*-
"""
Created on Sun Mar 20 20:09:56 2022

@author: wwww
"""
import numpy as np
from scipy.sparse import csc_matrix
from scipy.sparse.linalg import splu
A = csc_matrix([[1., 0., 0.], [5., 0., 2.], [0., -1., 0.]], dtype=float)
B = splu(A)
x = np.array([1., 2., 3.], dtype=float)
B.solve(x)



from asyn.SharedWorker import *
from asyn.ipc_matrix import *
from jsobj import *
import time

class Test(object):
    def __init__(self,d):
        
        ix=d.d.ix;
        
        self.ix=ix;
        self.x=ix.value;
        #time.sleep(1)
    def pp(self):
        
        self.x[:]=3*self.x;
        #time.sleep(1)
        return 'AAA'
    



if __name__=='__main__':
    from multiprocessing import freeze_support
    freeze_support()
    x=np.array([1j,2,3,4]);
    ix=ipc_array(x);
    iy=ipc_array(np.array([1,2]));
    y=ix.value;
    print(y)
    d=arg2jso(ix=ix);
    d=arg2jso(d=d)
    ws=SharedWorker(Test)(d);
    #f=ws(d).response;
    #print(f.result,f.error)
    r=ws.call('pp').result;
    y=ix.value;
    print(y)
