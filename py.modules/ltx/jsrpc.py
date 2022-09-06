# -*- coding: utf-8 -*-
"""
Created on Sat Jun 03 18:06:31 2017

@author: dich
"""
import jsonrpc.jsonclass as jc
import jsonrpc.sparse_marshal
import ltx.mm as mm
from ltx.js import jscript,ShowFigures;


class jscript_rpc(object):    
    
    def __init__(self,**rpc):
        
        class _jfs_(object):
            pass
        self.id=0;
        self.f=_jfs_(); 
        
        if '__parent__' in rpc:
            cp=rpc['__parent__']
            self.parent=cp;
            self.js=cp.js;
            self.jo=jo=cp.jo;
            self.jsmode=rpc.get('jsmode',cp.jsmode)
            
            r=rpc.get('pipe',{})
            srpc=jc.encode(r);            
            self.__channel__=jo('$$[0].pipe(JSON.parse($$[1]))',cp.__channel__,srpc);
            
        else:
            self.js=js=jscript(f_ipc=False);
            self.jo=jo=js.jo;
            self.jsmode=rpc.get('jsmode',False)
            
            
            s="\n".join(
            [
            'require("jsonrpc",global)',
            'new RPC_channel($$[0])'
            #'rpc_channel=RPC_channel($$[0])'
            ]
            );                     
            
            srpc=jc.encode(rpc);
            
            self.__channel__=jo(s,srpc);
        
        
        
    def __send_recv__(self,params=[],name='globalEval'):
         self.id+=1;
         msg={
         'id':self.id,
         'method':name,
         'params':params
         };
         
         smsg=jc.encode(msg);
         smsg=self.jo('$$[0].send_recv_raw($$[1])',self.__channel__,smsg);        
         msg=jc.decode(smsg);
         if msg.get('error',None):             
             raise Exception(msg['error']);
             
         r=msg.get('result',None)
         
         if self.jsmode:
             r=jc.jsobject_decode(r);
         
         return r;  
         
         
    def __call__(self,*params):              
        return self.__send_recv__(params);        
        
    def call(self,*params):
        params=[p for p in params];
        method=params.pop(0);
        return self.__send_recv__(params,method);
    
    def close(self):
        self.jo('$$[0].close()',self.__channel__);
        
    def pipe(self,**opts):
        return jscript_rpc(__parent__=self,pipe=opts)
    
    
    def functor(self,name,fname=None,init=None,prop_get='',prop_put=''):
        
        
        
        class jsfunctor(object):
            def __init__(self,jsrpc,name,fcache):
                self.name=name;
                self.__fcache__=fcache;
                self.__gcache__=fcache+'.'+prop_get;
                self.__pcache__=fcache+'.'+prop_put;
                self.js=jsrpc.js;
                self.rpc=self.jsrpc=jsrpc
            def __call__(self,*la):
                #la=(self.name,)+la;
                la=(self.__fcache__,)+la;                
                return self.jsrpc.call(*la);
            def __getitem__(self,key):                
                return self.jsrpc.call(self.__gcache__,key);
            def __setitem__(self,key,val):                
                return self.jsrpc.call(self.__pcache__,key,val);

        
        i=self('__functors__.push(['+name+'][0])')-1
        
        fcache='(__functors__['+str(i)+'])';
        
        if init:
            self(init);
        if(fname=='*'):
            fname=name;
        jf=jsfunctor(self,name,fcache);
        
        if(fname):
            self.f.__dict__[fname]=jf;
            
        return jf;
    def functors_tree(self,name,fname=None,init=None,defmethod=None,prop_get='',prop_put=''):
        f=self.functor(name,fname,init,prop_get,prop_put);        
        name=f.__fcache__;
        ff=self('__function_list__('+name+')');
        if defmethod and (defmethod in ff):
            #f=self.functor(name+'.'+defmethod,fname);
            f.__fcache__=f.__fcache__+'.'+defmethod;
        for n in ff:
            f.__dict__[n]=self.functors_tree(name+'.'+n);
        return f;
    rfunctor=r_functor=functors_tree
    
    def folder2sandbox(self,folder):
        return self.jo('$$[0].folder2sandbox($$[1])',self.__channel__,folder); 
    
    def push_sandbox(self,folder):
        #self.jo('$$[0].exec("global.__stack__.push($$[0])",dir2buffer($$[1]))',self.__channel__,folder);
        return self.jo('$$[0].push_sandbox($$[1])',self.__channel__,folder);
        
        
    
    
    def show(self,s):
        from IPython.display import SVG
        SVG(data=s)
        #import tempfile,shutil
        #from IPython.display import SVG
        #fd=tempfile.gettempdir()+'/svg000.svg';
        #self.jo('StringToFile($$[0],$$[1])',s,fd);        
        #SVG(fd);
        #fn=self.jo('StringToTmpFile($$[0])',s);        
        #shutil.copyfile(fn,fd)
        
        
    
    
jsrpc=jscript_rpc;
        