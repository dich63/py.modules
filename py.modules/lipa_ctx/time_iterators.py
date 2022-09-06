import numpy as np
from jsonrpc.jsonclass import *
import copy

class steps(object):
    def __init__(self,lst):
        self.lst=lst=copy.deepcopy(lst);        
        n=0; 
        self.ts=ts=[];
        self.tt=tt=[];
        self.to=to=[];
        self.dts=dts=[];
        self._times=None;
        ti=[0.0];
        t0=0;
        n=0;
        for i in lst:
            dt=i[0];
            tr=i[1];
            ti+=[tr];
            dts.append(dt);
            nt=int(np.ceil((tr-t0)/dt));            
            ts.append((dt,nt));
            t=[ t0+k*dt for k in range(1,nt+1)];            
            to.append(jsobject({"dt":dt,"nt":nt,"tt":t,"n":n}));
            n+=1;
            tt+=t;
            t0=tt[-1];
        self.ti=ti[0:-1];

    @property
    def times(self):
        if self._times is None:
            self._times=np.array(self.tt,dtype = np.double);
        return self._times
            
    def __iter__(self):
        #print('iter...')
        return iter(self.to)
    


