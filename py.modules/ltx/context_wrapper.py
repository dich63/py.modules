# -*- coding: utf-8 -*-
"""
Created on Mon Dec 05 05:11:11 2016

@author: dich
"""
from ctypes import *
from ctypes.wintypes import BOOL
import os
#from .mm import *;



class external_callbacks_t(Structure):
    _fields_ = [
    ("name",c_char_p),
    ("proc",c_void_p),
    ("attr",c_int)
    ];


long_pvoid_p=CFUNCTYPE(c_long,c_void_p);
long_pvoid_pvoid_p=CFUNCTYPE(c_long,c_void_p,c_void_p);
void_pvoid_p=CFUNCTYPE(None,c_void_p);
long_pvoid_p_proc=CFUNCTYPE(c_long,c_void_p,void_pvoid_p);

long_pvoid_p_proc_pvoid_p=CFUNCTYPE(c_long,c_void_p,void_pvoid_p,c_void_p);


class context_holder_t(Structure):
    _fields_ = [    
    ("addref",long_pvoid_p),
    ("release",long_pvoid_p),
    ("unwrap_context",long_pvoid_pvoid_p),
    ("wrap_context",long_pvoid_p_proc_pvoid_p),
    ("link",long_pvoid_pvoid_p),
    ("clear_links",long_pvoid_p)    
    ];

pcontext_holder_t=POINTER(context_holder_t);

def get_lib():
    p=os.environ['context_wrapper'];
    return cdll.LoadLibrary(p);
    
        

_get_context_utils=get_lib().get_context_utils
_get_context_utils.restype =c_void_p

context_utils=context_holder_t.from_address(_get_context_utils())
NullFunct=void_pvoid_p(0)
