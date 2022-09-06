# -*- coding: utf-8 -*-
"""
Created on Wed Sep 29 18:03:00 2021

@author: wwww
"""

def host_vm():
    import ltx.mm  as ltx
    return ltx.external()(0);

def host_vm_call(*d):    
    from ltx.js_vm_call import vm_call
    return vm_call(host_vm(),*d);
   
def host_vm_functor(sf):
    return lambda *d : host_vm_call(*((sf,)+d ))

def host_vm_detach(r):
    import ltx.mm  as ltx
    ltx.external().result=r;
    

def array2com(a):
    from ltx.mm import mm_buffer_from_array
    return mm_buffer_from_array(f).obj
    
    

if __name__=='__main__':
    
    printf=host_vm_functor('printf')
    jseval=host_vm_functor('globalEval')
    
    printf('modify  pyargs from python\n')
        
    
    import numpy as np    
    f=np.random.randn(10,20,30);
    
    jseval('pyargs=arguments[0]',[f,{"a":11,"f":f}])
    
    from ltx.mm import mm_buffer_from_array
    #
    import code;code.interact(local=globals())    
    host_vm_detach(array2com(f));
    
    
    
    