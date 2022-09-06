# -*- coding: utf-8 -*-
"""
Created on Thu Jan 17 05:06:14 2019

@author: wwww
"""
try:
    import ltx.mm as mm;
    from jsonrpc.jsonclass import config;
    import numpy as np;
    from jsonrpc.ndarray_marshal import js2np
    
    
    def objref_unmarshal(data,flink):
        if data[0:7]=='objref:':
            o=mm.bindObject(data);
            b=mm.mm_buffer(o);
            r=b.toarray(copy=not flink);
            if not flink:
                b.detach();
            return r;
        else:
            return None;
    
    
    config.classes.ipc_marshal['unmarshal']= objref_unmarshal; 
    try:
        stub_cache=mm.bindObject('stub.cache');    
        '''
        def objref_marshal(a): 
            global stubcache
            ma=mm.mm_buffer_from_array(a) 
            ma.detach();
            return stubcache(ma.obj);
        
        config.classes.ipc_marshal['marshal']= objref_marshal;       
        '''
        config.classes.ipc_marshal['marshal'] = lambda a: stub_cache(mm.mm_buffer_from_array(a).detach())
        config.classes.ipc_marshal['cache']   = lambda  : stub_cache(-1,0);
        
    except:
        pass
except:
    pass



 


    