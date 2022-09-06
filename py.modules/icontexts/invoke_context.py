#import env
import os,sys;
from ctypes import *
#



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




'''
from icontexts import context as ctx
iclib=ctx.load_library_with_dir(os.environ['invoke_context_lib'])
create_context=iclib.create_context
create_context.rstype=c_int32
create_context.argtypes = (c_void_p,long_pvoid_int64_pvoid_p,c_void_p)
'''



