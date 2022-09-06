# -*- coding: utf-8 -*-
"""
Created on Wed Jan 26 16:05:22 2022

@author: wwww
"""
import sys,os
from ctypes import *
import numpy as np


flib=cdll.LoadLibrary('E:\\projects\\repos\\Dll3\\x64\\Debug\\Dll3.dll')

ptr_at = lambda x: x.ctypes._data

my_fortran_function=flib.my_fortran_function
my_fortran_function.restype = c_int32
my_fortran_function.argtypes = (c_int32,c_double,c_void_p,c_void_p)

N=50
c=-10.00e0;
x=np.arange(N,dtype=np.float64)+1
y=np.zeros(N,dtype=x.dtype);
print('pid=',os.getpid())
err=my_fortran_function(N,c,ptr_at(x),ptr_at(y));

print('err=',err)
print('c=',c)
print('x=',x)
print('y=',y)