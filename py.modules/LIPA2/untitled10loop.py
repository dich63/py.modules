# -*- coding: utf-8 -*-
"""
Created on Tue Apr 12 19:01:58 2022

@author: wwww
"""

import pickle
import multiprocessing

class Tricky:
    def __init__(self,x):
        self.data=x

    def __setstate__(self, d):
        print('setstate happening')
        self.data = 10

    def __getstate__(self):
        return self.data
        print('getstate happening')

def report(ar,q):
    q.put(ar.data)

global_ctx=globals()
def geval(iexpr,*_lp,**_kp): 
    #print('lp=',_lp)
    global_ctx['_r']=None
    global_ctx['_lp']=_lp
    global_ctx['_kp']=_kp
    exec(iexpr,global_ctx);
    r=global_ctx['_r'];
    del global_ctx['_r'],global_ctx['_lp'],global_ctx['_kp']
    return r


def gelo(*_lp,**_kp): 
    return [_lp,_kp]

    
def eval_loop(qread,qwrite):
    while(1):
        pass
    
    

if __name__ == '__main__':
    
    from utils import *
    tic()
    
    for k in range(1000):
        geval('_r=gelo(*_lp,**_kp)',1111,112,ttt=11,r=2)
        
    print('t=',toc()*1000)
    
    ii=compile('_r=gelo(*_lp,**_kp)','','exec',optimize=-1)
    tic()
    
    for k in range(1000*1000):
        geval(ii,1111,112,ttt=11,r=2)
        
    print('tc=',toc())
    
    tic()
    for k in range(1000*1000):
        _r=gelo(1111,112,ttt=11,r=2)        
    print('te=',toc())
    raise SystemExit(0)
        
    
    
    geval('ddq=111');
    ar = Tricky(5)
    
    q = multiprocessing.JoinableQueue()
    p = multiprocessing.Process(target=report, args=(ar, q))
    print('now starting process')
    p.start()
    print('now joining process')
    p.join()
    print('now getting results from queue')
    print(q.get())
    print('now getting pickle dumps')
    print(pickle.loads(pickle.dumps(ar)).data)   