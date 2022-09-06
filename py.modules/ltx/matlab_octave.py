# -*- coding: utf-8 -*-
"""
Created on Mon Jan 21 01:57:09 2019

@author: wwww
"""
import ltx.js as ltxjs
import ltx.ipc_marshal

def octave(f_ipc=1):
    return ltxjs.jscript(f_ipc=f_ipc).rfunctor('oc=require("oe").OctaveEngine()',call='feval');
def matlab(f_ipc=1):
    return ltxjs.jscript(f_ipc=f_ipc).rfunctor('oc=require("me").MatlabEngine()',call='feval');