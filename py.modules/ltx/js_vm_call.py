# -*- coding: utf-8 -*-

import jsonrpc.jsonclass as jc
import jsonrpc.marshals
from jsobj import js_getter_setter

def vm_call(vm,*d):    
    
    f=jc.ipc_mode();    
    if f:
        jc.ipc_cache_clear();    
        
    json_str=jc.encode(d);
    json_str=vm("require('jsoncall').jsoncall($$[0],$$[1])",json_str,not f);
    if not (json_str is None):
        r=jc.decode(json_str);        
        return r;

def vm_functor(vm,sf):
    return lambda *d : vm_call(vm,*((sf,)+d ))

    
def vm_context(vm):
    def getter(k):
                
        f,v=vm_call(vm,\
        'function(k){\
           var v=globalEval("["+k+"][0]");\
           return (typeof v=="function")?[true,null]:[false,v]\
         }'\
        ,k);
                    
        return vm_functor(vm,'['+str(k)+'][0]') if f else v;         
        #return vm_functor(vm,'global["'+str(k)+'"]') if f else v;         
        
    def setter(k,v):
        vm_call(vm,'globalEval','global[$$[0]]=$$[1];""',k,v);
        
    def caller(*d):
        d=("globalEval",)+d;
        return vm_call(vm,*d);
        
    
    return js_getter_setter(getter,setter,caller);
        
def vm_methods_context(vm):
    def getter(k):
        return vm_functor(vm,str(k));                         
        #return vm_functor(vm,'global["'+str(k)+'"]');                         
    def caller(*d):
        d=("globalEval",)+d;
        return vm_call(vm,*d);        
    
    return js_getter_setter(getter=getter,caller=caller);

