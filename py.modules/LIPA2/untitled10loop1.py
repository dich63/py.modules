# -*- coding: utf-8 -*-
"""
Created on Tue Apr 12 19:01:58 2022

@author: wwww
"""

import pickle
import multiprocessing
import multiprocessing as mp
import threading
import ctypes
from asyn.ipc_array import *

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
    q.put(ar)


def geval(iexpr,*_lp,**_kp): 
    try:
        global_ctx=globals()
        #print('lp=',_lp)
        global_ctx['_r']=None
        global_ctx['_lp']=_lp
        global_ctx['_kp']=_kp
        exec(iexpr,global_ctx);
        r=global_ctx['_r'];
        del global_ctx['_r'],global_ctx['_lp'],global_ctx['_kp']
        return r,0,''
    except Exception as e:
        return None,1,str(e)
        pass



def gelo(*_lp,**_kp): 
    return [_lp,_kp]

class response_t(object):
    def __init__(self,result=None,error=None):
        self.result,self.error=result,error
    

    
    
def eval_loop(n,begin_barrier,end_barrier,sh_icmd,queue):    
    
    invokers=();    
    lpf,kpf=(),{};
    
    icmd=sh_icmd.value;
    status=icmd[(2+n):]
    
    while(1):
        try:
            begin_barrier.wait();                               
            
            if icmd[0]==0:
                try:
                    status[0]=invokers[icmd[1]]();
                except:
                    status[0]=-1;
                
            else:            
                try:
                    if icmd[0]==1:
                        expr,lp,kp=queue.get()
                        r=geval(expr,lp,kp);
                        queue.put(response_t(r));                                             
                    else:
                        global_ctx=globals()
                        expr,lpf,kpf=queue.get()
                        func=eval(expr,global_ctx);                        
                        invokers+=(lambda : func(*lp,**kp) ,)
                        queue.put(response_t(len(invokers)-1));                    
                except BaseException as e:
                        queue.put(response_t(error=e));
                        
            end_barrier.wait();
        except threading.BrokenBarrierError:
            begin_barrier.reset();
            end_barrier.reset();
            pass




class workers_batch_t(object):  
    
    class worker_t(object):
        def __init__(self,owner):
            self.owner=owner;           
        
    
    
    def __init__(self,n):
        n=n;
        begin_barrier=mp.Barrier(n+1);
        end_barrier=mp.Barrier(n+1);
        queue=mp.Queue();#JoinableQueue()
        sh_icmd=ipc_array(np.zeros(2+n,dtype=np.int32));
        
        icmd=sh_icmd.value[:1];
        invok_id=sh_icmd.value[1:2];
        statuses=sh_icmd.value[2:];
        
        workers=[ mp.Process(target=eval_loop, args=(k,begin_barrier,end_barrier,sh_icmd,queue))  for k in range(n)]
        
        self._set_locals(locals());        
        pass
    
    def _set_locals(self,ctx):        
        self.__dict__.update(ctx)
        
    
    
        

    
    

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