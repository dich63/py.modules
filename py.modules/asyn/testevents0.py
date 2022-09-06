# -*- coding: utf-8 -*-
"""
Created on Thu Apr  7 02:16:08 2022

@author: wwww
"""

import os
import multiprocessing  as mp
import time
import inspect
import numpy as np
from asyn.ipc_matrix import *
from jsobj import *



def worker_proc(b_event,e_event,ix,iy):
    x=ix.value;
    y=iy.value;
    
    while(True):
        b_event.wait();
        b_event.clear();
        y[:]+=1;
        e_event.set()

if __name__=='__main__':
    from multiprocessing import freeze_support
    freeze_support()
    
    
    
    
    
    from utils import *
    
    
    
    
    
    
    #e=np.array([1j,2,3,4]);
    x=1j*np.random.rand(1028,128)
    y=np.zeros((14,14),dtype='double')
    
    ix=ipc_array(x);
    iy=ipc_array(y);
    y=iy.value;
    b_event=mp.Event();
    e_event=mp.Event();
    
    process=mp.Process(target=worker_proc,args=(b_event,e_event,ix,iy));
    #process.daemon=daemon;
    b_event.clear();
    
    
    process.start();
    for k in range(10):
        e_event.clear();
        b_event.set();
        
        e_event.wait();
    tic()
    N=20000;
    for k in range(N):
        e_event.clear();
        b_event.set();
        
        e_event.wait();
    t=toc('')
    print(y[1:4,1:4])
    print('tmp[Ms]=',1000*t/N)
    tic()
    for k in range(N):
        y[:]+=1;
    t=toc('')
    print(y[1:4,1:4])
    print('t[Ms]=',1000*t/N)
