# -*- coding: utf-8 -*-
"""
Created on Thu Jul 14 14:24:05 2022

@author: DICH
"""

import ltx.mm  as ltx

#_external=ltx.external();
__vmref__=None #lambda *l: raise Exception('host vm not found')

def set_vm_ref(v):
    global __vmref__
    __vmref__=v

def call(*d):
    from ltx.js_vm_call import vm_call
    return vm_call(__vmref__(),*d);
   
def functor(sf):
    return lambda *d : host_vm_call(*((sf,)+d ))

def context():
    from ltx.js_vm_call import vm_context
    return vm_context(__vmref__());

def methods():
    from ltx.js_vm_call import vm_methods_context
    return vm_methods_context(__vmref__());
