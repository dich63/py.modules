import os
import multiprocessing  as mp
import time
import inspect
import numpy as np

#from . import pp_affinity as affinity
#import asyn.pp_affinity as affinity
#f_affinity=True


f_affinity=False
try:
    #from .affinity import affinity
    from . import pp_affinity as affinity
    f_affinity=True
    
except ImportError:
    pass


def get_affinity(htdisable=True,even=True):
    if f_affinity:
        affinity_mask=affinity.get_process_affinity_mask(0)
    else:
        affinity_mask=(2**mp.cpu_count()-1)
    if htdisable:
        affinity_mask=ht_disable(affinity_mask,even)
    return affinity_mask

def ht_disable(affinity_mask,even=True):
    m=0x5555555555555555&((2**mp.cpu_count())-1)
    m*=2**even;
    affinity_mask&=m;
    return affinity_mask


def affinity2cpus(affinity_mask=get_affinity()):
    cpus=[];
    affinity_mask=np.int64(affinity_mask)
    n=0;
    while affinity_mask:
        if affinity_mask&1:
            cpus.append(n)
        n+=1;
        affinity_mask>>=1;
    return cpus;

def cpus2affinity(cpus):
    affinity_mask=np.int64(0)
    if f_affinity:
        l=len(cpus);
        for n in range(l):
            affinity_mask+=2**cpus[n]
    return affinity_mask;


def thread_packets(datas,cpus=affinity2cpus(),fpack=True):
    datas=[d  for d in datas ]
    cn=len(cpus)
    packet=[([],2**(cpus[n]-0))  for n in range(cn) ]
    while(datas):
        for n in range(cn):
            if datas:
                d=datas.pop();
                packet[n][0].append(d);
            else:
                break
    return  [p for p in packet if p[0]] if fpack else packet;


def thread_packets_2(datas,np,cpus=affinity2cpus(),fpack=True):
    datas=[d  for d in datas ]
    cn=len(cpus)
    packet=[([],2**(cpus[n]-0))  for n in range(cn) ]
    while(datas):
        for n in range(cn):
            if datas:
                d=datas.pop();
                packet[n][0].append(d);
            else:
                break
    return  [p for p in packet if p[0]] if fpack else packet;


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
    def m1(self,l=0,r=2):
        return self.a.value[l:r];

class Tic:
    def __init__(self):
        self.start();
    def start(self):
        self.t=time.perf_counter();
    def sec(self):
        return time.perf_counter();-self.t;

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


def worker_proc(qin,qout,class_type,args_c,kw_c):
    obj=None
    lastop=None
    _tic_=Tic();
    
    
    try:        
        qin.get()
        obj=class_type(*args_c,**kw_c);
        qout.put(Response(result=class_type.__name__,tp=_tic_.sec()))
        qin.task_done()
        
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
                lastop=getattr(obj,op);
                qout.put(Response(result=True,tp=_tic_.sec()))
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

def register_methods(names,this,obj):
    for n in names:
        setattr(this,n,lambda *k,**o : obj.call(n,*k,**o))
    return this;


class AsynObject:
    def __init__(self,owner):
        self.worker=owner;
    def register_methods(self,names):
        return register_methods(names,self,self.worker);
    def join(self):
        if not hasattr(self,'_response'):
            w=self.worker;
            if not w.process.is_alive():
                raise Exception('process is dead');
            w.qin.join();
            self._response=w.qout.get();
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
    def call(self,name,*args,**kwargs):
        return self.worker.call(name,*args,**kwargs);

    def fake_call(self,name,*args,**kwargs):
        return self.worker.fake_call(name);
        
    def set_call(self,name,*args,**kwargs):
        return self.worker.fake_call(name);


    def recall(self,*args,**kwargs):
        return self.worker.recall(*args,**kwargs);

    def close(self):
        return self.worker.close();


class AsynGroup:
    def clear(self):
        (self.aos_s,self.aos,self.rsp)=([],[],[])
    def reset(self):
        (self.aos,self.rsp)=([],[])
    def __init__(self):
        self.clear();
    def __iter__(self):
        return iter(self.aos_s);
        
    def __lshift__(self,ao):
        if self.rsp:
            self.rsp=[]
            #raise Exception('AsynList already Responses');
        self.aos.append(ao)
        return self
    def join(self):
        if not self.rsp:
            self.rsp=[o.response for o in self.aos]
            self.aos_s,self.aos=self.aos,[]
        return self
    @property
    def responses(self):
        return self.join().rsp
    @property
    def results(self):
        return [r.result for r in self.responses]
    @property
    def errors(self):
        return [r.error for r in self.responses]




class SharedWorker:
    
    def new(self,*args,**kwargs):
        if  hasattr(self,'process'):
            raise Exception('object already exists');
            
        daemon=self.daemon;
        affinity_mask=self.affinity_mask;
        
        self.qin=qin=mp.JoinableQueue();
        self.qout=qout=mp.JoinableQueue();
        process=mp.Process(target=worker_proc,args=(qin,qout,self.classtype,args,kwargs))
        self.process=process;
        process.daemon=daemon;
        process.start();
        n=mp.cpu_count();
        maxaff=(2**n)-1;
        job=self.job;
        

        if f_affinity:
            if not job is None:
                affinity.assign_process_to_job(job,process.pid)
                #affinity.assign_process_to_job(job,0)
            if affinity_mask:
                affinity_mask&=maxaff;
                affinity.set_process_affinity_mask(process.pid,affinity_mask);
        
        return self.__ops__();

    def __call__(self,*args,**kwargs):
        return self.new(*args,**kwargs);
        
        
    def __init__(self,classtype,affinity_mask=0,daemon=False,job=None):
        self.daemon=daemon;
        self.affinity_mask=affinity_mask;
        self.classtype=classtype;
        self.job=job
    def __del__(self):
        self.terminate()
    
    def terminate(self):
        if  hasattr(self,'process'):
            self.process.terminate();
    
    def __ops__(self,n=None,op=None,args=None,kwargs=None):
        if self.process.is_alive():
            if n==None:
                self.qin.put(True)
            else:
                self.qin.put((n,op,(args,kwargs)))
            return AsynObject(self);
        else:
            raise Exception('process is dead');

    def call(self,name,*args,**kwargs):
        return self.__ops__(1,name,args,kwargs);
        
    def fake_call(self,name):
        return self.__ops__(2,name,[],{});


    def recall(self,*args,**kwargs):
        return self.__ops__(0,None,args,kwargs);


    def close(self):
        return self.__ops__(-1,None,(),{});

