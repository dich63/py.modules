# -*- coding: utf-8 -*-
"""
Created on Mon Dec 05 05:11:11 2016

@author: dich
"""
from ctypes import *
from p23 import *
import os
import sys

# low level

class external_callbacks_t(Structure):
    _fields_ = [
    ("name",c_char_p),
    ("proc",c_void_p),
    ("attr",c_int)
    ];

long_void_p=CFUNCTYPE(c_long);
long_pvoid_p=CFUNCTYPE(c_long,c_void_p);
#long_pvoid_p=CFUNCTYPE(c_long,c_void_p,c_int);
long_pvoid_pvoid_p=CFUNCTYPE(c_long,c_void_p,c_void_p);
void_pvoid_p=CFUNCTYPE(None,c_void_p);
long_pvoid_p_proc=CFUNCTYPE(c_long,c_void_p,void_pvoid_p);
long_pvoid_pvoid_p_proc=CFUNCTYPE(c_long,c_void_p,c_void_p,void_pvoid_p);
long_void_proc=CFUNCTYPE(c_int32)

#long_pvoid_p_proc_pvoid_p=CFUNCTYPE(c_long,c_void_p,void_pvoid_p,c_void_p);
long_pvoid_p_proc_pvoid_p=CFUNCTYPE(c_long,c_void_p,c_void_p,c_void_p);
def from_param(cls, obj):
    if obj is None:
        return None # return a NULL pointer
    from ctypes import _CFuncPtr
    return _CFuncPtr.from_param(obj)

void_pvoid_p.from_param=from_param;

class context_holder_t(Structure):
    _fields_ = [    
    ("addref",long_pvoid_p),
    ("release",long_pvoid_p),
    ("unwrap_context",long_pvoid_pvoid_p),
    ("wrap_context",long_pvoid_p_proc_pvoid_p),
    ("link",long_pvoid_pvoid_p),
    ("clear_links",long_pvoid_p),
    ("global_ref_count",long_void_proc),      
    ("unlink",long_pvoid_pvoid_p),
    ("wrap_context_ex",long_pvoid_pvoid_p_proc),
##--------------------------------------------------------------------
    ("create_handle_ex",long_pvoid_pvoid_p_proc),
    ("get_context",long_pvoid_pvoid_p),
    ("create_handle",long_pvoid_p_proc_pvoid_p),
    
    ("tss_link_context",long_pvoid_p),
    ("tss_unlink_context",long_pvoid_p),
    ("tss_clear_links",long_void_p),
    ("tss_onexit",long_pvoid_p_proc),
    ("load_lib_path",long_pvoid_pvoid_p)
    ];

pcontext_holder_t=POINTER(context_holder_t);

def get_lib():
    p=os.environ['context_wrapper'];    
    return cdll.LoadLibrary(p);
    
        

_get_context_utils=get_lib().get_context_utils
_get_context_utils.restype =c_void_p
#_get_context_utils();
context_utils=context_holder_t.from_address(_get_context_utils())
NullFunct=void_pvoid_p.from_address(0)

ginit=context_utils.global_ref_count()
_g_dict__={};

def get_gc():
    global _g_dict__;
    return  _g_dict__;

# high level    
class context_handle_t(Structure):
    None
   
p_context_handle_t=POINTER(context_handle_t);        


class test:
    def __del__(self):
        print('del id=',id(self) );


def __onexit__(idh):
    #global _g_dict__
    #print('release:'+hex(idh))
    r=context_utils.global_ref_count()
    dbg_print('[grefc={0} ]:release: {1}'.format(r,hex(idh)))
    #print('[grefc={0} ]:release: {1}'.format(r,hex(idh)))
    #print('[grefc={0} ]:release: BEGIN {1}'.format(r,hex(idh)))
    
    get_gc().pop(idh);
    #print('[grefc={0} ]:release END: {1}'.format(r,hex(idh)))

p__onexit__=void_pvoid_p( __onexit__);




def _CHECK_ERR(hr=-1,msg='Error:'):
    if hr:
        raise Exception(msg+hex(hr));


def wrap_context(o):
    pp=c_void_p(0);
    o=(o,);
    i=id(o);
    get_gc()[i]=o;
    #hr=context_utils.wrap_context(i,p__onexit__,byref(pp))
    #if hr:
    #   raise Exception('Error wrap_context code='+hex(hr));
    _CHECK_ERR(context_utils.wrap_context(i,p__onexit__,byref(pp)),
               'Error wrap_context code=')
        
    r=context_utils.global_ref_count()
    dbg_print('[grefc={0} ]:addref: {1}'.format(r,hex(i)))
    
    return pp;
        
def unwrap_context(h,p=c_void_p(0)):
    '''
    hr=context_utils.unwrap_context(h,byref(p))
    if hr:
        raise Exception('Error unwrap_context code='+hex(hr));
    '''
    _CHECK_ERR(context_utils.unwrap_context(h,byref(p)),
               'Error unwrap_context code=');
    return p;  
    
def release_context(h):
    c=context_utils.release(h);
    if c<0:
        raise Exception('Error release_context code='+hex(c));
    return c;
    
def link_context(h,h2):
    hr=context_utils.link(h,h2)
    if hr:
        raise Exception('Error link_context code='+hex(hr));
    
    

def clear_context_links(h):
    return context_utils.clear_links(h);      

add_lib_path= lambda p: 0
rem_lib_path= lambda p: 0

def __dbg__print(s):
    print(s);
    
dbg_print=__dbg__print
if sys.platform=='win32':
    
    #windll.kernel32.SetDefaultDllDirectories(0x00001e00);
    _AddDllDirectory=windll.kernel32.AddDllDirectory
    _AddDllDirectory.restype=c_longlong;
    _RemoveDllDirectory=windll.kernel32.RemoveDllDirectory
    _RemoveDllDirectory.argtypes=(c_longlong,);
    
    _OutputDebugString=windll.kernel32.OutputDebugStringW
    _OutputDebugString.argtypes=(c_wchar_p,);
    
    #
    dbg_print=lambda s :_OutputDebugString(unicode(s+'\n'));
    
    add_lib_path =lambda p: _AddDllDirectory(unicode(p))
    rem_lib_path = _RemoveDllDirectory

def load_library_with_dir_2(p):
    cookie=add_lib_path(p+'/../');
    lib=cdll.LoadLibrary(p);
    rem_lib_path(cookie)
    return lib

def load_library_with_dir_0(p):
    pw=unicode(p+'/../');
    windll.kernel32.SetDllDirectoryW(pw);
    '''
    sp=os.getenv('path');
    pw+=unicode(';'+sp)
    windll.kernel32.SetEnvironmentVariableW(unicode(path),unicode(pw))
    '''
    try:
        lib=cdll.LoadLibrary(p);
    except:
        windll.kernel32.SetDllDirectoryW(0)
        #windll.kernel32.SetEnvironmentVariableW(unicode(path),unicode(sp))
        raise
        
    #windll.kernel32.SetEnvironmentVariableW(unicode(path),unicode(sp))
    windll.kernel32.SetDllDirectoryW(0);
    return lib

def load_library_with_dir_1(p):
    ps=os.environ['PATH'];
    try:
        pn=p+'/../;'+ps;
        os.environ['PATH']=pn
        lib=cdll.LoadLibrary(p);
    except:
        os.environ['PATH']=ps
        raise
    os.environ['PATH']=ps    
    return lib

load_library_with_dir=load_library_with_dir_0
#'''
