# -*- coding: utf-8 -*-
"""
Created on Wed Jan 19 09:32:01 2022

@author: wwww
"""
import matplotlib.pyplot as plt
from ltx.mm import ltx_create_callback
from ltx.ltx_srv import eval_json_class_callback,eval_callback
from ltx.js import jscript
import jsonrpc.jsonclass as jc

g=jscript(objref='host.script')
a=jscript(objref='script:')
wsleep=g.functor('require("tools").sleep')



class ccr(object):
    def __init__(self,fun):
        self.fun=fun;
    def __call__(self,*la):
        return self.fun(*la);
    

def jsonwraper(site,com_args):
    print(site,com_args)
    lp=site('$$[0]',com_args)
    print(lp)
    r=fun(*lp)
    site('__result__=$$[0];""',r);
    cr=site.jo('__result__')
    return cr;

def fun(lp):
    
    return lp.len

def fun2(a,b):
    return {'a':a,'b':b}
#o=ltx_create_callback(eval_json_class_callback)

cc=ccr(lambda x : 'aaaaaaaaaaaaaa');

o=ltx_create_callback(fun )

mplot=g.functor('matlab.plot')
msurf=g.functor('matlab.surf')
mfeval=g.functor('matlab.feval')

o=ltx_create_callback(cc )
oo=g.wrap_func(fun2)
g.sleep(100)
#o=ltx_create_callback(lambda args: jsonwraper(g,fun,args) )
#o=ltx_create_callback(lambda args: args )
g.jo('fun2=$$[0]',oo)
g.jo('py=$$[0]',o)

