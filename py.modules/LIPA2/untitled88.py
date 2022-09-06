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


class Tricky:
    def __init__(self,x):
        self.data=x

    def __setstate__(self,d):
        self.data=-10

    def __getstate__(self):
        return {}    

class Test(object):
    def __init__(self,d):
        
        ix=d.d.ix;
        
        self.ix=ix;
        self.x=ix.value;
        #time.sleep(1)
    def pp(self,n=4,d=None):
        self.dd=d;
        #
        #self.x[:]=self.x+2;
        #
        self.x[0,0]+=n;
        #self.x[:,:]=0.9*self.x+0.1;
        #time.sleep(1)
        return 'AAA'
    
    def getM(self,n=9):
        tt=Tricky(1111);
        return 11;
        




if __name__=='__main__':
    from multiprocessing import freeze_support
    freeze_support()
    
    from utils import *
    #e=np.array([1j,2,3,4]);
    x=1j*np.random.rand(28,28)
    
    ix=ipc_array(x);
    iy=ipc_array(np.array([1,2]));
    y=ix.value;
    print(y)
    d=arg2jso(ix=ix);
    d=arg2jso(d=d)
    SW=SharedWorker(Test)   
    
    ws=SW(d);
    
    r=ws.call('pp',0).result;
    print('pp:',r)
    jj=arg2jso(p=11,g=33)
    #f=ws(d).response;
    #print(f.result,f.error)
    r=ws.call('getM').result;
    print('getM:',r)
    r=ws.call('getM').result;
    print('getM:',r)
    print('')
    tic()
    for nn in range(10000):
        r=ws.recall(nn,jj).result;
    toc('mp::')
    
    tt=Test(d)
    tic()
    for nn in range(10000):
        r=tt.pp(nn,jj);
    toc('in::')
    
    y=ix.value;
    #print(y)
    #toc('2:')