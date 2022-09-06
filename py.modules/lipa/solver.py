import multiprocessing  as mp
import time
import affinity

class Tic:
    def __init__(self):
        self.start();
    def start(self):
        self.t=time.time()
    def sec(self):
        return time.time()-self.t;



def worker(qin,qout=None):
    tupl=[None,(None,),None]
    while(True):
        r=qin.get()
        cmd=r[0]
        op=r[1]
        if cmd==1:
            op=getattr(tupl[0],op);
            r=op(*tupl[1],**tupl[2])
            if qout!=None:
                qout.put(r)

        elif cmd==2:
            args=r[2][0];
            kw=r[2][1];
            #print('start...')
            tupl[0]=op(*args,**kw);
            #print('complete...')
            tupl[1]=r[3][0]
            tupl[2]=r[3][1]
        else:
            break
        qin.task_done()

    qin.task_done()

class asyn_object:
    def __init__(self,qin,quot):
        self.qin=qin;
        self.qout=qout;
    def wait(self):
        self.qin.join();
        self.result=self.qout.get();


class process_ipc:
    def __init__(self,affinity_mask=0):
        self.qin=qin=mp.JoinableQueue();
        self.qout=qout=mp.JoinableQueue();
        self.process=process=mp.Process(target=ls.worker,args=(qin,qout))
        process.daemon=True
        process.start();
        n=mp.cpu_count();
        maxaff=(2**n)-1;
        if affinity_mask:
            affinity_mask&=maxaff;
            affinity.set_process_affinity_mask(process.pid,affinity_mask);

        

