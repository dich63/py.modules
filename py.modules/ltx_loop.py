"""
import os
import sys
import ltx.mm  as ltx
"""

import ltx.ipc_marshal
#print(dir(ltx))

global_ctx=globals();
local_ctx={};

__fjsobject__=False
__argv_base__=1


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
    
    if f:
        l=len(aa)
        global_ctx['result']=None;
        global_ctx['arguments']=aa;
        for k in range(l):
            global_ctx['_'+str(k+__argv_base__)]=aa[k]
        return l
    else:
        l=aa
        if '_r' in global_ctx:
            del global_ctx['_r']
            
        del global_ctx['result'];
        
        global_ctx['arguments']=aa;
        for k in range(l):
            del global_ctx['_'+str(k+__argv_base__)]
            
        
def eval_json_class_callback(arguments):
    #print('eval_json_class_callback...')
    global __fjsobject__
    #print(__fjsobject__)
    _s=arguments.shift;
    #print('s=',s,arguments[0])    
#    global_ctx['arguments']=jc.decode(arguments[0],__fjsobject__);
#    global_ctx['result']=None;
    _l=__reset_arguments__(jc.decode(arguments[0],__fjsobject__),True)
    try:
        exec(_s,global_ctx)
    except SystemExit:
        print('exit');
    
    r=global_ctx.get('_r',global_ctx['result']);
    #print(r)
    r=jc.encode(r);
    __reset_arguments__(_l)
    
    #print(r)
#    del global_ctx['arguments'];
#    del global_ctx['result'];
    return r

#os.system("pause")

if __name__=='__main__':
    from multiprocessing import freeze_support
    freeze_support()
    import os
    import sys
    import ltx.mm  as ltx
    
    pcallback=eval_callback
    
    external=ltx.external();
    #l=external.length
    
    fc=external(0,False)
    __fjsobject__=external(1,True)
    __argv_base__=external(2,__argv_base__)
    
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
