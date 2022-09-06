#
from ltx.external import externObj


external=externObj()
    

shape=external[0]
echo=external[1]

import numpy as np
f=np.random.rand(*(shape))+1j*np.random.rand(*(shape));

# callback optional
try:
    import os
    from  ltx.js_vm_call import vm_functor,vm_context
    vm=external.raw_args(2)    
    vmthis=vm_context(vm)
    
    #external.result=[1j*np.arange(111),echo];
    import code;code.interact(local=globals())  
    cprintf=vm_functor(vm,'cprintf')    
    cprintf(0xF4,"Ups from python [pid=%d]\n",os.getpid())
except:
    pass
    

external.result=[f,echo];

