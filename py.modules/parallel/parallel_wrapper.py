# -*- coding: utf-8 -*-
"""
Created on Mon Dec 05 05:11:11 2016

@author: dich
"""

from ctypes import *
import os
import numpy as np
import env
import context_wrapper.context_wrapper as ctx
import jsonrpc.jsonclass as jc
# low level

long_pvoid_int64_pvoid_p=CFUNCTYPE(c_long,c_void_p,c_longlong,c_void_p);
long_pvoid_p=CFUNCTYPE(c_long,c_void_p);
p_void_p=POINTER(c_void_p);

class context_list_t(Structure):
    pass

p_context_list_t  = POINTER(context_list_t)
pp_context_list_t  = POINTER(p_context_list_t)

class invoke_proc_u_t(Union):
    _pack_=8
    _fields_ = [
       ("pinvoke",long_pvoid_int64_pvoid_p),
       ("pinvoke0",long_pvoid_p),
       ("ptr_invoke",c_void_p)
       ]

class cmd_u_t(Union):
    _pack_=8
    _fields_ = [
       ("ccmd",c_char*16),
       ("wcmd",c_wchar*8),
       ("icmd",c_longlong),
       ("pcmd",c_void_p)
       ]
       

       
    
class base_context_t(Structure):
    _anonymous_ = ("u")
    _pack_=8
    _fields_ = [    
    ("weak_ref_handle",c_void_p),
    ("uuid",c_char*16),
    ("u",invoke_proc_u_t),  
    ("context",c_void_p)  
    ];

p_base_context_t  = POINTER(base_context_t)
pp_base_context_t  = POINTER(p_base_context_t)


class params_u_t(Union):
    _pack_=8
    _fields_ = [    
       ("pcontext",p_base_context_t),
       ("params",c_void_p),       
       ("iparams",c_longlong),
       ("fvalue",c_double),
       ("pclist",p_context_list_t),
       ]
 

class invoke_context_t(base_context_t):
    _anonymous_ = ("uc","up")
    _pack_=8
    _fields_ =[
    ("uc",cmd_u_t),
    ("up",params_u_t),
    ("status",c_longlong)
    ]

p_invoke_context_t  = POINTER(invoke_context_t)
pp_invoke_context_t  = POINTER(p_invoke_context_t)


context_list_t._pack_=8
context_list_t._fields_=[
("pp",pp_invoke_context_t),
("count",c_longlong),
("rep",c_longlong),
("once",c_longlong)
]



class buffer_flat_params_t(Structure):    
    _pack_=8
    _fields_ = [
    ("op",c_longlong),
    ("N",c_longlong),
    ("count",c_longlong),          
    ("pp",c_void_p)  
    ];





def get_lib():
    #
    p=os.environ['parallel_wrapper_lib'];
    #p="V:/Projects/tbb_wrapper/x64/Debug/tbb_wrapper.dll"    
    #p="V:/Projects/tbb_wrapper/x64/Release/tbb_wrapper.dll"
    '''
    try:
        pass
        #cdll.LoadLibrary(p+'/../tbb.dll');
    except Exception:
        pass
    '''
    #
    return ctx.load_library_with_dir(p)
    #return cdll.LoadLibrary(p);
    
pwl=get_lib();

create_lin_comb=pwl.create_lin_comb;
create_lin_comb.rstype=c_int32
create_lin_comb.argtypes = (pp_base_context_t,)

create_LA_op=pwl.create_LA_op;
create_LA_op.rstype=c_int32
create_LA_op.argtypes = (pp_base_context_t,)


create_buffer_context=pwl.create_buffer_context
create_buffer_context.rstype=c_int32
create_buffer_context.argtypes = (c_longlong,pp_base_context_t)

#_buffer_from_memory = pythonapi.PyBuffer_FromMemory
#_buffer_from_memory.restype=py_object

create_ref_context=pwl.create_ref_context;
create_ref_context.rstype=c_int32
create_ref_context.argtypes = (pp_base_context_t,pp_base_context_t)



create_invoke_context=pwl.create_invoke_context;
create_invoke_context.rstype=c_int32
create_invoke_context.argtypes = (p_base_context_t,pp_invoke_context_t)





create_invoke_context_indirect=pwl.create_invoke_context_indirect;
create_invoke_context_indirect.rstype=c_int32
create_invoke_context_indirect.argtypes = (c_void_p,long_pvoid_int64_pvoid_p,long_pvoid_p,pp_invoke_context_t)


_create_invoke_context_indirect=pwl.create_invoke_context_indirect;
_create_invoke_context_indirect.rstype=c_int32
_create_invoke_context_indirect.argtypes = (c_void_p,c_void_p,c_void_p,pp_base_context_t)

    
create_invoke_context_batch=pwl.create_invoke_context_batch;
create_invoke_context_batch.rstype=c_int32
#create_invoke_context_batch.argtypes = (c_int32,pp_invoke_context_t,c_int32,c_int32,pp_invoke_context_t)
create_invoke_context_batch.argtypes = (c_int64,c_void_p,c_int64,c_int64,c_int64,pp_invoke_context_t)

invoke_context=pwl.invoke_context;
invoke_context.rstype=c_int32
invoke_context.argtypes = (p_base_context_t,)


invoke_context_batch=pwl.invoke_context_batch;
invoke_context_batch.rstype=c_int32
invoke_context_batch.argtypes = (c_int32,p_base_context_t,c_int32)


release_context=pwl.release_context;
release_context.rstype=c_int32
release_context.argtypes = (c_void_p,)

addref_context=pwl.addref_context;
addref_context.rstype=c_int32
addref_context.argtypes = (c_void_p,)


create_parallel_group=pwl.create_parallel_group;
create_parallel_group.rstype=c_int32
create_parallel_group.argtypes = (pp_base_context_t,)


create_parallel_group_ex=pwl.create_parallel_group_ex;
create_parallel_group_ex.rstype=c_int32
create_parallel_group_ex.argtypes = (c_longlong,pp_base_context_t)


call_context=pwl.call_context;
call_context.rstype=c_int32
#call_context.argtypes = (p_base_context_t,c_longlong,c_void_p)
call_context.argtypes = (c_void_p,c_longlong,c_void_p)


tls_buffer=pwl.tls_buffer;
tls_buffer.rstype=c_int32
#call_context.argtypes = (p_base_context_t,c_longlong,c_void_p)
tls_buffer.argtypes = (c_longlong,c_void_p)

tls_buffer_clear=pwl.tls_clear;
tls_buffer_clear.rstype=c_int32
tls_buffer_clear.argtypes = ()

_byref=byref;


def ptr_array(pp,offset=0):    
    l=[p.ctypes.data for p in pp]
    for k in range(offset):
        l.insert(0,0);
    return np.array(l,dtype=np.uint64);
        

def struct_def(struct,opts):
    for p in struct._fields_:
        k=p[0]
        #if opts.has_key(p[0]):            
        if k in opts.keys():
            struct.__setattr__(k,opts[k] );
    return struct;
        

def refcount_context(p):
    addref_context(p);
    return release_context(p);
    
def icmd(s):
    u=cmd_u_t() 
    u.icmd=0
    #u.ccmd=bytes(s[0:15]);
    if not s is None:
        t=type(s)
        if t in jc.string_types:        
            u.ccmd=s[0:15].encode('utf8');
        else:
            u.icmd=s;        
    
    return u.icmd;
    
    
    
def ctx_check_err(hr):
    if hr:
        raise Exception('base_context: error:'+hex(hr))
    return 0;
    
def cast_context(p,byref=False,faddress=False):
    if  issubclass(type(p),hcontext):
        p= p.byref() if byref  else p.ptr;
        return addressof(p.contents) if faddress else p;
    else:
        return p;
        
def cast_params(p,byref=False):
    if  issubclass(type(p),hcontext):
        p= p.byref() if byref  else p.ptr;
        return addressof(p.contents);
    elif issubclass(type(p),np.ndarray):        
        return p.ctypes.data;
    return addressof(p) if byref  else p;


class hcontext(object):
    def __init__(self,ptr=None,ptype=None, constructor=None,addref=False):
        ptr=self._ptr=c_void_p() if ptr is None else ptr;        
        self.ptype=type(ptr) if ptype is None else ptype;
        if not constructor is None:
            ctx_check_err(constructor(self.byref())); 
        if addref:
            addref_context(self._ptr);
                    
    
    def __del__(self):        
        try:
            #if bool(i):
            
            i=release_context(self._ptr);

            ctx.dbg_print(str(type(self))+' ~ ptr_context: {0}'.format(int(i)))
            #print(str(type(self))+' ~ ptr_context: {0}'.format(int(i)))
            if(i<0):
                raise Exception(str(type(self))+' ~ ptr_context: {0}'.format(int(i)))

            
        except Exception:
            pass
        
    def __nonzero__(self):
        return bool(self._ptr);
        
    def __bool__(self):
        return bool(self._ptr);

            
    @property
    def refcount(self):
        return refcount_context(self._ptr)
    
    @property
    def ptr(self):
        return cast(self._ptr,self.ptype);
    @property
    def address(self):
        return addressof(self.ptr.contents) if self._ptr else 0;
    @property
    def value(self):
        return self.ptr.contents;
    
    @property
    def handle(self):
        return self.ptr.contents.weak_ref_handle;
        
    def cast_ptr(self,ptype=None):        
        if ptype is None:
            ptype=self.ptype;
        return cast(self._ptr,ptype);
    
    def byref(self,ptype=None):
        if ptype is None:
            ptype=self.ptype;
        return cast(addressof(self._ptr),POINTER(ptype));
        
    def ref_context(self):
        pp=self.byref();
        pref=p_base_context_t();
        ctx_check_err(create_ref_context(pp,_byref(pref)))        
        h=hcontext(pref);
        h.link(self);
        return h;
        
    def __call__(self,cmd=None,params=None,byref=False):
        i=icmd(cmd)
        p=cast_params(params,byref);
        ptr=self.ptr
        return call_context(ptr,i,p);
    
    def create_invoke_context(self,cmd=None,params=None,byref=False,links=None):
        return hinvoke_context(self,cmd,params,byref,links);
        
    def link(self,*links):
        hls=ctx.wrap_context(links);
        ctx.link_context(self.handle,hls)
        ctx.release_context(hls);
        return self;
        


class hinvoke_context(hcontext):
    
    def __init__(self,context=None,cmd=None,params=None,byref=False,links=None,raw_context=None):
        #if issubclass(type(context),hcontext):
         #   context=context.ptr;
        
        if not raw_context is None:
            super(hinvoke_context,self).__init__(raw_context);
            return
            
        if type(context) is pp_base_context_t:
            ref,context=context,p_base_context_t();
            ctx_check_err(create_ref_context(ref,_byref(context)))            
        else:
            context=cast_context(context)
        
        p=p_invoke_context_t()
        ctx_check_err(create_invoke_context(context,_byref(p)))
        super(hinvoke_context,self).__init__(p);
        
        self.link(params,links)
        
        self.value.icmd=icmd(cmd)              
        self.value.params=cast_params(params,byref=byref);
        self.value.status=-111222333444;        
        

        
    def __del__(self):
        try:
            hcontext.__del__(self)
        except Exception:
            pass
        #super(hinvoke_context,self).__del__();
        
    @property    
    def status(self):
        return self.value.status;
    def __call__(self):
        return invoke_context(self._ptr);
           

class hinvoke_batch(hinvoke_context):
    def __init__(self,contexts=None,asyn=False,links=None,rep=1,once=0):
        l=len(contexts)
        
        cc=np.empty(shape=(l,),dtype=np.int64);
        for i in range(l):
            cc[i]=cast_params(contexts[i]);
        pc=p_invoke_context_t();
        ctx_check_err(create_invoke_context_batch(l,cc.ctypes.data,asyn,rep,once,_byref(pc)));
        
        super(hinvoke_batch,self).__init__(raw_context=pc);        
        self.link(cc,links);        
        
        self.value.icmd=asyn;
    @property
    def asyn(self):
        return bool(self.value.icmd);    
    @asyn.setter
    def asyn(self,asyn):
        self.value.icmd=bool(asyn);    

    def __del__(self):
        try:
            hinvoke_context.__del__(self)
        except Exception:
            pass

        

def create_empty_context():
    p=p_base_context_t();
    _create_invoke_context_indirect(None,None,None,byref(p));
    return p;

class pp_group_ex(object):
    def __init__(self,mode=1):

        self.mode=mode;        
        g=p_base_context_t();
        if create_parallel_group_ex(mode,byref(g)):
            raise Exception('_create_parallel_group error')
        self.group=g;
        self.iwait=hinvoke_context(g,'join');       

        
    def __del__(self):
        
        try:
        
            i=release_context(self.group);
            #print('~ ~pp_group:',i)
            dbg_print('~ ~pp_group: {0}'.format(i))
        except Exception:
            pass

    def __call__(self,invoke_context):                
        invoke_context=cast_context(invoke_context)
        return call_context(self.group,ord('r'),invoke_context);

    def join(self):
        return call_context(self.group,ord('j'),p_invoke_context_t());

    def link(self,*links):
        if not self.mode==1:
            raise Exception('pp_group_ex link support only mode==1')
        h=ctx.wrap_context(links);
        ctx.link_context(self.group.contents.weak_ref_handle,h)
        ctx.release_context(h);
        return self;


    
class pp_group(object):
    def __init__(self,sync=False):
        self.sync=sync;
        if not sync:       
            g=p_base_context_t();
            if create_parallel_group(byref(g)):
                raise Exception('_create_parallel_group error')
            self.group=g;
            self.iwait=hinvoke_context(g,'join');
        else:            
            self.iwait=hinvoke_context();

        #self.join_invoker()
    def __del__(self):
        
        try:
            if not self.sync:
                i=release_context(self.group);
                #print('~ ~pp_group:',i)
                dbg_print('~ ~pp_group: {0}'.format(i))
        except Exception:
            pass

    def __call__(self,invoke_context):        
        if self.sync:
            return invoke_context();
        else:
            invoke_context=cast_context(invoke_context)
            return call_context(self.group,ord('r'),invoke_context);
    def join(self):
         if not self.sync:      
             return call_context(self.group,ord('j'),p_invoke_context_t());
         else:
             return 0;
        
    def join_invoker(self):        
        
        '''
        pi=p_invoke_context_t()        
        call_context(self.group,ord('j'),byref(pi));
        return hinvoke_context(raw_context=pi);
        '''
        
class buffer_flat_invoker(hcontext):
    def __init__(self,xx,xx_unpack=None): 
        if xx_unpack is None:
            xx_unpack=xx;
        self.xx_pack=self.xx=xx;
        self.xx_unpack=xx_unpack;
        count = len(xx);
        if not count:
            raise Exception(' hbuffer_scatter_context error: must be count>0')
        N=xx[0].size*xx[0].itemsize;  
        
        sz=N*count;
        
        bfp=buffer_flat_params_t();
        bfpu=buffer_flat_params_t();
        bfpu.count=bfp.count=count;
        bfpu.N=bfp.N=N;
        bfp.op=0;
        bfpu.op=1;
        self.hbuffer=hbuffer=hcontext(ptype=p_base_context_t,constructor =lambda p : create_buffer_context(sz,p) );   
        
        lpk=ptr_array(xx);
        luk=ptr_array(xx_unpack);
        
        bfp.pp=lpk.ctypes.data;
        bfpu.pp=luk.ctypes.data;
        
        
        ipack=hbuffer.create_invoke_context('scatter',bfp,byref=True,links=(xx,lpk));
        iunpack=hbuffer.create_invoke_context('scatter',bfpu,byref=True,links=(xx_unpack,luk));
        
        self.ipack=ipack;
        self.iunpack=iunpack;
        
    def __del__(self):
        try:
            hcontext.__del__(self)
        except Exception:
            pass
        
        
        
        #copy.
    


def create_buf_pw(N,fzero=True, dtype=np.complex128):        
    
    dtype=np.dtype(dtype)
    sz=dtype.itemsize*np.prod(N);
    hb=hcontext(ptype=p_base_context_t,constructor =lambda p : create_buffer_context(sz,p) );   
    #hb('realloc',sz);
    if(fzero):
        hb('zero');
        
    #buf=_buffer_from_memory(hb.address,sz);
    p=c_void_p();
    hb('get',byref(p));
    buf=(c_char*(sz)).from_address(p.value);
    buf.__holder__=hb
    r=np.frombuffer(buf,dtype=dtype)
    r=r.reshape(N);    
    return r;
    

def create_buf(N,fzero=True, dtype=np.complex128):
            
    #return create_buf_pw(N,fzero=fzero,dtype=dtype);
    #
    return np.zeros(N,dtype=dtype);

def create_buf_like(xx,fzero=True, dtype=np.complex128):
    if type(xx) in [list,tuple]:
        return [create_buf_like(x,fzero,dtype) for x in xx];
    else:
        create= np.zeros_like if fzero else np.empty_like;
    return create(xx,dtype=dtype);
