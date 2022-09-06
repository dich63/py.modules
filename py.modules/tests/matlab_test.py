# -*- coding: utf-8 -*-
"""
Created on Mon Jun 05 10:53:25 2017

@author: Administrator
"""

  

#from IPython import get_ipython;get_ipython().magic('reset -sf') 

import sys,os
sys.path.append('v:/ipc/py.modules')






from ltx.js import jscript;

js=jscript()

        
ml=js.functors_tree('ml','*',"ml=new (require('matlab_engine2').MatlabEngine)({rpc:{host:'dg'}})")
#ml=js.functors_tree('ml','*',"ml=require('matlab_engine2').MatlabEngine({rpc:{host:'dg'}})")
#ml=js.functors_tree('ml','*',"ml=require('matlab_engine2').MatlabEngine({rpc:{host:'78.138.137.166'}})")
#ml=js.functors_tree('ml','*',"ml=require('matlab_engine2').MatlabEngine()")

#ml=js.functors_tree('ml','*',"ml=require('matlab_engine2').MatlabEngine({rpc:{host:'localhost'}})")
#js.f.ml.execute('testLIPA_z')

#

print ml.execute('testLIPA_ctx2')
q=ml.show_figures()


ml.execute('close all')

print ml.execute('testLIPA_ctx')
#q=js.f.ml.show_figures()
q=ml.show_figures()

ml.feval('close', 'all')
'''
print ml.execute('testLIPA_z')
q=ml.show_figures()
'''
am=ml.feval('sin',np.array([1,2,3j]))
am

#ml.terminate() 