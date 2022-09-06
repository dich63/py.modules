# -*- coding: utf-8 -*-
"""
Created on Wed Jun 22 17:14:45 2016

@author: dich
"""
import os
import numpy.ctypeslib as npct

midlib=None

lib_path=os.getcwd();

try:
    midlib = npct.load_library('MID2D.dll',lib_path)
except Exception:
    lib_path=__file__+'/../bin';
    midlib = npct.load_library('MID2D.dll',lib_path)


    