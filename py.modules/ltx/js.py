# -*- coding: utf-8 -*-
"""
Created on Sat Jun 03 18:06:31 2017

@author: dich
"""
import jsonrpc.jsonclass as jc

import ltx.mm as mm
import ltx.ipc_marshal
#import jsonrpc.sparse_marshal
import jsonrpc.marshals

def fun00(lp):
    
    return 'aaaaaaaaaaaaaaa'

class func_holder_t(object):
    def __init__(self,func):
        self.func=func;
    def __call__(self,*la):
        return self.func(*la);


def _functor_cache_op(idf,prefix='',postfix=''):
    s='if(typeof __functor_cache__=="undefined" ){ __functor_cache__={} };'
    n='__functor_cache__['+str(idf)+']'
    s+=prefix+n+postfix
    return s,n;
    
    

def _reduce_name(jo,name,idf=0):
    
    
    #jo=self.jo;            
    name=name.lstrip();
    i=name.find(':=')                
    if i>=0:
        head,tail=name[:i],name[(i+2):]            
        head=head.rstrip();
        name=head            
        
        if i==0:
            s,name=_functor_cache_op(idf,postfix='='+tail)
        else:
            s=name+'='+tail;                        
            
        
        jo(s)  
        
    return name;


class jscript(object):    
    
    def __init__(self,opts='',objref=None,rpc={},f_ipc=True,jo=None):
        
        class _jfs_(object):
            pass
        
        if jo is None:
            if objref is None:
                jo=mm.jsvm(opts);           
            else:
                jo=mm.bindObject(objref);
            
        self.jo=jo
        
        jo('require("jsonrpc",global);');
        self.jeval=jo('json_eval_ipc');
        self.jcall=jo('json_call_ipc');
        self.sleep=jo('require("tools").sleep')      
        
        self.JsonClassEncodeIPC=jo('JsonClassEncodeIPC');
        self.JsonClassDecode=jo('JsonClassDecode');
        
        
        self.f=_jfs_();
        self.f_ipc=f_ipc;
        self._cbl=[];
        
    def __call__(self,*la):
        f_ipc=self.f_ipc;
        s=jc.encode_ipc(la,f_ipc);
        #print(s)
        
        s=self.jeval(s,f_ipc);
        return jc.decode(s) if not s is None else s;
    
    def call(self,*la):
        f_ipc=self.f_ipc;
        s=jc.encode_ipc(la,f_ipc);
        s=self.jcall(s,f_ipc);
        #print(type(s))        
        return jc.decode(s) if not s is None else s;
    
    def wrap_func(self,func,name=None): 
        
        f_ipc=self.f_ipc;
        decoder=self.JsonClassDecode;
        encoder=self.JsonClassEncodeIPC;
        
        
        def marshaller(args):
            
            len=args.len;
            #print('len=',len)
            
            jeval=self.jeval
            params=();
            
            #print('for')
            
            for k in range(len):
                s=encoder(args[k],f_ipc)
                params+=( jc.decode(s),);
                
            result=func(*params);
            #print('forend')
            result=jc.encode_ipc(result,f_ipc);
            return decoder(result);
            
            #return '++111!!!'
        
        cm=func_holder_t(marshaller);
        self._cbl+=[cm];
        #disp=mm.ltx_create_callback(marshaller);
        #disp=mm.ltx_create_callback(fun00);
        disp=mm.ltx_create_callback(cm);
        
        
        if type(name)==str:
            self.jo('global.'+name+'=$$[0];""',disp);
        
        return disp;
    
    
    
    
    def functor(self,name,fname=None,init=None): 
        
        class jsfunctor(object):           
            
            
            def __init__(self,js,name):
                
                self.name=_reduce_name(js.jo,name,id(self));
                self.js=js;
                #self._holds=hold_list(js.jo,name);
                
            def __del__(self):
                try:
                    s=_functor_cache_op(id(self),prefix='delete ');
                    self.js.jo(s)
                except:
                    pass;
                
                
            def __call__(self,*la):
                
                la=(self.name,)+la;
                return self.js.call(*la);
            
            def __repr__(self):
                d=self.__dict__
                if len(d):                    
                    tn=str(type(self))                    
                    s='functor:{\n'
                    for k,v in d.items():
                        if str(type(v))==tn:
                            s+=str(k)+': functor\n'
                            #s+=str(k)+':'+str(type(v))+'\n'                    
                    s+='}'            
                else:
                    s='functor:{}';            
                return s;

        
        if init:
            self.jo(init);
            
        #name=self.__reduce_name__(name);
        jf=jsfunctor(self,name);
        
        name=jf.name;
            
        if(fname=='*'):
            fname=name;
        
        
        if(fname):
            self.f.__dict__[fname]=jf;
            
        return jf;
    
    def functors_tree(self,name,fname=None,init=None,call=None):        
        
        
        
        jo=self.jo
        name=_reduce_name(jo,name,0);
        nc=':='+name
        if call is None:
            f=self.functor(nc,fname,init);
        else:
            f=self.functor(nc+'.'+call,fname,init);
            
        ff=self('__function_list__('+name+')');
        
        for n in ff:
            f.__dict__[n]=self.functor(nc+'.'+n);
        
        f._com=jo(name);   
        
        jo(_functor_cache_op(0,prefix='delete ')[0])    
        return f;
    rfunctor=r_functor=functors_tree
    

__global_script__=False;

def global_js():
    global __global_script__;
    if not __global_script__:
        __global_script__=jscript();
    return __global_script__;

def ShowFigures(f):
    global_js().call('require("matlab-engine2").ShowFigures',f)


def jslink(name):
    return jscript(objref=name,f_ipc=1);

def jsbind(parsestr,f_ipc=1):
    jo=mm.bindObject(parsestr);
    return jscript(jo=jo,f_ipc=f_ipc);

def jshost(job_nested=0):
    j=jscript();
    jo=j.jo("require('utils').global_vm($$[0])",job_nested);
    return jscript(jo=jo,f_ipc=1);
    
        