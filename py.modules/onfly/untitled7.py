# -*- coding: utf-8 -*-
"""
Created on Fri Feb 18 02:21:42 2022

@author: wwww
"""
import logging
import httpimport
from httpimport import add_remote_repo, remove_remote_repo;
httpimport.INSECURE=1
httpimport.RELOAD=1
logger=httpimport.logger
httpimport.logger.setLevel(logging.INFO)
url='http://127.0.0.1:8000/'
#mm=httpimport.load('mf',url);
ii=add_remote_repo('mf',url)

import mf.c
#ii=add_remote_repo('mf.c.A',url)

import mf.c.A.B.c