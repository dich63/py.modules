# -*- coding: utf-8 -*-
"""
Created on Tue Jun 19 00:10:03 2018

@author: dich
"""

import os,sys;sys.path.append(os.getenv('ltx_js_root')+'py.modules')

import numpy as np


from  ltx.jsrpc import jsrpc as rpc

carat_host='172.20.20.106'

carat=rpc(host=carat_host)



carat.folder2sandbox('T:/.felix/nout-2/*')

matlab=carat.rfunctor('require("matlab_engine2").MatlabEngine()',defmethod='exec',prop_get='getvar',prop_put='put_var')

print(matlab("cd(getenv('sandbox_dir'));cd"))

print(matlab("mainSolve8sTime2"))

tJ,Jf=matlab['tJ'],matlab['Jf']


import matplotlib.pyplot as plt

plt.plot(np.abs(tJ),np.abs(Jf[:,1]),np.abs(tJ),np.abs(Jf[:,2]))
