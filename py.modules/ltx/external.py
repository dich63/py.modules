#

import ltx.mm  as ltx
import jsonrpc.jsonclass as jc
#import jsonrpc.sparse_marshal
import jsonrpc.ndarray_marshal
from ltx.js import jscript;


class externObj(object):
    def __init__(this,f_ipc=True):        
        this._exobj=exobj=ltx.external();
        
        this.f_ipc=f_ipc;
        
        this._fjson=fjson=exobj.fjson
        
        objref="script" if fjson else "srv:script"         
        
        this._js=js=jscript(objref=objref,f_ipc=f_ipc);  
        js.jo('external=$$[0]',exobj)
        if not fjson:               
            js.jo("__stub__=bindObject('stub');__stub__(externRef());''")          
        #js.jo("require('register').register_vm('ee');''")       
    @property
    def count(this):
        return this._exobj.len
    @property
    def pid(this):
        return this._exobj.pid;
    
    @property
    def fjson(this):
        return this._fjson;
    
    
    @property
    def result(this):
        #return this._js("result");
        return None
    
    
    @result.setter
    def result(this,value):
        if this.fjson:
            s=jc.encode_ipc(value,this.f_ipc);
            this._exobj.result=s;
        else:
            js=this._js
            js("result=$$[0];null",value);               
            this._exobj.result=js.jo("result");
        
    def __getitem__(this,i):
          return this._js("external[$$[0]]",i);
    def raw_args(this,i):
        return this._exobj[i];
       
if __name__=='__main__':
    
    
    external=externObj()
    
    import os
    from  js_vm_call import vm_functor
    vm=external.raw_args(0)    
    #import code;code.interact(local=globals())  
    cprintf=vm_functor(vm,'cprintf')    
    cprintf(0xF4,"Ups from python [pid=%d]\n",os.getpid())
    
    shape=external[1]
    echo=external[2]
    
    import numpy as np
    f=np.random.rand(*(shape))+1j*np.random.rand(*(shape));
    
    external.result=[f,echo];
    