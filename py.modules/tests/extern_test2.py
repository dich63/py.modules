#
from utils import *
import jsonrpc.jsonclass as jc
from ltx.external import externObj


external=externObj(1)
    

shape=external[0]
echo=external[1]
fipc=external[3]
import numpy as np


tic()
f=np.random.rand(*(shape))+1j*np.random.rand(*(shape));
lf=f.reshape(-1);
t=toc()


# callback optional
if 1:
    try:
        import os
        from  ltx.js_vm_call import vm_functor,vm_context
        jc.ipc_mode(fipc)    
        js_this=vm_context(external.raw_args(2))   
        
        #    import code;code.interact(local=globals())  
        js_this.rf=lf
        js_this.cprintf(0xF4,"Ups from python [pid=%d]\n",os.getpid())
        import code;code.interact(local=globals())  
    except:
        pass
    

external.result=[f,echo,t];

