# -*- coding: utf-8 -*-
"""
Created on Sun Aug  8 03:50:43 2021

@author: wwww
"""
from parallel.sparse import *

class sps_invoker(sparse_solver_invoker):

    def __init__(self,sA,**opts):
        super(sps_invoker,self).__init__(sA,opts);
        self._ifs=hinvoke_batch((self.factorize,self.solve));
    def __call__(self,y):
        self.y=y;
        self._ifs();
        return self.y.copy();
    @property
    def handle(self):
        return self._ifs;

