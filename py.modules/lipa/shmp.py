# -*- coding: utf-8 -*-
"""
Created on Sat Jan  2 21:25:59 2016

@author: dich6_000
"""
import ctypes
import numpy as np
import multiprocessing  as mp
from multiprocessing import heap
from multiprocessing.reduction import ForkingPickler
from multiprocessing.context import assert_spawning


def reduce_ctype(obj):
    assert_spawning(obj)
    return rebuild_ctype, (type(obj), obj._wrapper, None)

"""
    if isinstance(obj, ctypes.Array):
        return rebuild_ctype, (obj._type_, obj._wrapper, obj._length_)
    else:
        return rebuild_ctype, (type(obj), obj._wrapper, None)
"""

def rebuild_ctype(type_,wrapper,length=None):
    ForkingPickler.register(type_,reduce_ctype)
    buf = wrapper.create_memoryview()
    obj = type_.from_buffer(buf)
    obj._wrapper = wrapper
    return obj


def sh_array(count):
    """
    if count<0x7FFFFFF:
        return mp.RawArray('b',count);
    
    count=np.int64(count);
    ir=count and 0x0FFF;
    cp=count>>12;
    if ir!=0:
        cp+=1;
    """
    ct=ctypes.c_byte*count;
    size=ctypes.sizeof(ct)
    wrapper=heap.BufferWrapper(size)
    return rebuild_ctype(ct,wrapper);



