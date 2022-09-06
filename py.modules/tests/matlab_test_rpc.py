# -*- coding: utf-8 -*-
"""
Created on Mon Jun 05 10:53:25 2017

@author: Administrator
"""

  

#from IPython import get_ipython
#get_ipython().magic('reset -sf') 

import sys,os
sys.path.append('v:/ipc/py.modules')





import numpy as np
from ltx.jsrpc import jsrpc#,ShowFigures;


#from IPython.display import SVG

def ShowFigures(ff):   
    import IPython.display as dd
    for f in ff:
        dd.display_svg(f[1],raw=True)
    

host='dg'
#
host='78.138.137.166'
#ml=jsrpc(host='78.138.137.166').functors_tree('ml','*',"ml=require('matlab_engine2').MatlabEngine()","exec")
ml=jsrpc(host=host).r_functor('ml','*',"ml=require('matlab_engine2').MatlabEngine()","exec")



#

print ml('testLIPA_ctx2')

ShowFigures(ml.figures())



ml('close all')

print ml('testLIPA_ctx')
#q=js.f.ml.show_figures()
ShowFigures(ml.figures())

ml.feval('close', 'all')
print ml.execute('testLIPA_z')


ShowFigures(ml.figures())

am=ml.feval('sin',np.array([1,2,3j]))
print am
ml.rpc.close()

#ml.terminate() 