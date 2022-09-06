import os
import multiprocessing  as mp
import affinity
import time
import inspect

class Testsw:
    def __init__(self,a=11,b=12,c=13,ipc_params=None):
        self.a=a;
        self.b=b;
        self.c=c;
        self.options=ipc_params;
    def dump(self,*args,**kwargs):
        return (self.a,self.b,self.c);
    def sh(self,*args,**kwargs):
        self.options[0].value[:]=os.getpid();
        return os.getpid();
    def m1(self,*args,**kwargs):
        return (args,kwargs);

class Tic:
    def __init__(self):
        self.start();
    def start(self):
        self.t=time.time()
    def sec(self):
        return time.time()-self.t;

class imap:
    def __init__(self):
        self._key=0;
        self._map={}
    def __getitem__(self,k):
        return self._map[k];

    def __len__(self):
        return self._key;

    def push(self,v):
        k=self._key;
        self._map[k]=v;
        self._key+=1;
        return k;


class Response:
    def __init__(self,result=None,error=None,tp=-1):
        self.error=error;
        self.result=result;
        self.tp=tp;


def worker(qin,qout,params):
    obj=None
    lastop=None
    _tic_=Tic();
    try:
        while(True):
            r=qin.get()
            _tic_.start()
            cmd=r[0]
            op=r[1]
            args=r[2][0];
            kw=r[2][1];
            if cmd==0:
                r=lastop(*args,**kw);
                qout.put(Response(result=r,tp=_tic_.sec()))
            elif cmd==1:
                lastop=getattr(obj,op);
                r=lastop(*args,**kw);
                qout.put(Response(result=r,tp=_tic_.sec()))
            elif cmd==2:
                if (params!=None) and ('ipc_params' in inspect.getargspec(op.__init__).args):
                    kw['ipc_params']=params;
                obj=op(*args,**kw);
                qout.put(Response(result=op.__name__,tp=_tic_.sec()))
            else:
                break
            qin.task_done()
            #qout.join()
            None
        None
        qout.put(Response(tp=_tic_.sec()));
    except Exception as e:
        qout.put(Response(error=e,tp=_tic_.sec()));
    None
    qin.task_done()
    #qout.join()


class asyn_object:
    def __init__(self,owner):
        self.owner=owner;
    def join(self):
        if not hasattr(self,'_response'):
            owner=self.owner;
            if not owner.process.is_alive():
                raise Exception('process is dead');
            owner.qin.join();
            self._response=owner.qout.get();
        return self
    @property
    def response(self):
        self.join();
        return self._response
    @property
    def result(self):
        r=self.response;
        if r.error:
            raise r.error;
        return r.result;






class shared_worker:
    def __init__(self,params=None,affinity_mask=0,daemon=False):
        self.qin=qin=mp.JoinableQueue();
        self.qout=qout=mp.JoinableQueue();
        self.process=process=mp.Process(target=worker,args=(qin,qout,params))
        process.daemon=daemon;
        process.start();
        n=mp.cpu_count();
        maxaff=(2**n)-1;
        if affinity_mask:
            affinity_mask&=maxaff;
            affinity.set_process_affinity_mask(process.pid,affinity_mask);
    def __del__(self):
        self.process.terminate();
        
    def __ops__(self,n,op,args,kwargs):
        if self.process.is_alive():
            self.qin.put((n,op,(args,kwargs)))
            return asyn_object(self);
        else:
            raise Exception('process is dead');

    def new(self,classtype,*args,**kwargs):
        return self.__ops__(2,classtype,args,kwargs);

    def call(self,name,*args,**kwargs):
        return self.__ops__(1,name,args,kwargs);

    def recall(self,*args,**kwargs):
        return self.__ops__(0,None,args,kwargs);


    def close(self):
        return self.__ops__(-1,None,(),{});





