import jsonrpc.sparse_marshal
import jsonrpc.jsonclass as jc
import sys,os
sys.path[0]=os.path.abspath(sys.path[0]+'/../');


import ltx.ipc_marshal
#print(dir(ltx))

global_ctx=globals();
local_ctx={};

__fjsobject__=False
__argv_base__=1
__vmref__=None

def globals_update_from(o):
    if not (type(o) in (dict,)):
        o=o.__dict__;    
    globals().update(o)

def host_vm_eval(*d):
    return __vmref__()(*d);


def host_vm_call(*d):
    from ltx.js_vm_call import vm_call
    return vm_call(__vmref__(),*d);
   
def host_vm_functor(sf):
    return lambda *d : host_vm_call(*((sf,)+d ))

def host_vm_context():
    from ltx.js_vm_call import vm_context
    return vm_context(__vmref__());

def host_vm_methods():
    from ltx.js_vm_call import vm_methods_context
    return vm_methods_context(__vmref__());



def set_module_path(p,ftop=True):    
    
    if type(p)==str:
        p=p.split(';');    
    
    if ftop:
        for i in reversed(p):
            sys.path.insert(0,i);           
    else:
        for i in p:
            sys.path.append(i);           
        
    return sys.path;
    

def exit(r=None):
    
    if not r is None:
        global_ctx['_r']=r;
        
    raise SystemExit();
    
    

def set_json_mode(m):
    global __fjsobject__
    __fjsobject__=m
    return __fjsobject__;
    
def plot2svg(plt):
    import sys
    if sys.version_info.major==3:
        import io
    else:
        import cStringIO as io
        
    buffer = io.StringIO()
    plt.savefig(buffer,format='svg')
    s=buffer.getvalue();    
    buffer.close()
    return s
    

def eval_callback(arguments):
    s=arguments.shift;
    if type(s)!=str:
        s=s.encode('ascii')
#    print(type(s))
#    print(s,'\n::len=',arguments.length)
    local_ctx['arguments']=arguments;
    local_ctx['result']=None;
    exec(s,globals(),local_ctx)
    #execfile(s)
    r=local_ctx['result'];
    del local_ctx['arguments'];
    del local_ctx['result'];
    return r

"""
def eval_json_class_callback(arguments):
    #print('eval_json_class_callback...')
    global __fjsobject__
    #print(__fjsobject__)
    s=arguments.shift;
    #print('s=',s,arguments[0])
    local_ctx['arguments']=jc.decode(arguments[0],__fjsobject__);
    local_ctx['result']=None;
    exec(s,global_ctx,local_ctx)
    r=local_ctx['result'];
    #print(r)
    r=jc.encode(r);
    #print(r)
    del local_ctx['arguments'];
    del local_ctx['result'];
    return r
"""
def __reset_arguments__(aa,f=False):
    def erase_global(*nn):
        for n in nn:
            if n in global_ctx:
                del global_ctx[n]    
    
    erase_global('_r','result');
    
    if f:
        l=len(aa)
        #global_ctx['_r']=None;        
        #global_ctx['result']=None;
        
        global_ctx['arguments']=aa;
        for k in range(l):
            global_ctx['_'+str(k+__argv_base__)]=aa[k]
        return l
    else:
        l=aa
        
        #[ erase_global('_'+str(k+__argv_base__)) for k in range(l)]
        erase_global(*[ '_'+str(k+__argv_base__) for k in range(l)]);
        
        '''
        if '_r' in global_ctx:
            del global_ctx['_r']
        if 'result' in global_ctx:
            del global_ctx['result'];
        
        global_ctx['arguments']=aa;
        for k in range(l):
            del global_ctx['_'+str(k+__argv_base__)]
            
        '''
        
def eval_json_class_callback(arguments):
    #print('eval_json_class_callback...')
    global __fjsobject__
    #print(__fjsobject__)
    _s=arguments.shift;    
    #print(arguments[0])
    #print('s=',s,arguments[0])    
#    global_ctx['arguments']=jc.decode(arguments[0],__fjsobject__);
#    global_ctx['result']=None;

    _l=__reset_arguments__(jc.decode(arguments[0],__fjsobject__),True)
    
    try:
        exec(_s,global_ctx)
    except SystemExit:
        print('exit');
    
    r=global_ctx.get('result',global_ctx.get('_r',None));
    #r=global_ctx.get('_r',global_ctx['result']);
    #print(r)
    r=jc.encode(r);
    
    __reset_arguments__(_l)
    
    #print(r)
#    del global_ctx['arguments'];
#    del global_ctx['result'];
    return r

def eval_json_class_callback2(*arguments):
    return eval_json_class_callback(arguments);

#os.system("pause")

if __name__=='__main__':
    from multiprocessing import freeze_support
    freeze_support()
    import os
    import sys
    import ltx.mm  as ltx
    #print(' ltx_py start...')
    pcallback=eval_callback
    
    external=ltx.external();
    #l=external.length
    
    fc=external(0,False)
    __fjsobject__=external(1,True)
    __argv_base__=external(2,__argv_base__)
    __vmref__=external(3,__vmref__)
    
    fc=1
    #print(fc)
    if fc:
            
        #import ltx.ipc_marshal
        import jsonrpc.sparse_marshal
        import jsonrpc.jsonclass as jc
        
        pcallback=eval_json_class_callback
    else:
        pcallback=eval_callback
        
#    pcallback=eval_compile_callback if fc else eval_callback
    hr=ltx.ltx_create_process_callback_loop(pcallback)

"""
    print 'finite...'
    quit()
    print 'Ok finite'
"""


#print(hr)
#os.system("pause")
"""
dtest=ltx.ltx_create_callback(eval_callback);
r=dtest('a=11\nb=22\nresult=a+b')
r=dtest('result=arguments[1]*arguments[2]',4,5)
r=dtest('result=a/b')
print(r)
"""
