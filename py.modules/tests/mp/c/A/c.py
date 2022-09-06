# -*- coding: utf-8 -*-
"""
Created on Tue Feb 15 18:15:40 2022

@author: public_user
"""
import os ;print('[ pid='+str(os.getpid())+'] module <<'+__name__+'>> locate: ',__file__)


def makeFilter(depth_step, depth_units, window_size_m,responses):
    return  'FilterExampleCall(depth_step, depth_units, window_size_m,responses)'