# -*- coding: utf-8 -*-
"""
Created on Sat Jun 25 01:37:04 2022

@author: DICH
"""
import jsonrpc.marshals
import jsonrpc.jsonclass as jc
import numpy as np

a=np.arange(11)
s=jc.encode(a)