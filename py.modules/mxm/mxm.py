# -*- coding: utf-8 -*-
"""
Created on Thu Oct 15 22:33:35 2020

@author: wwww
"""




from ltx.js import jslink

from jsobj import jso

def communicator(name='QtApp'):
    j=jslink(name)
    return jso({
    'readScalarLogData':j.functor('communicator.readScalarLogData'),
    'readDepthData':j.functor('communicator.readDepthData'),
    'writeLogData':j.functor('communicator.writeLogData'),
    'readMatrixLogData':j.functor('communicator.readMatrixLogData')    
    })

maxim=communicator;