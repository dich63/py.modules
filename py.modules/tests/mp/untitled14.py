# -*- coding: utf-8 -*-
"""
Created on Fri Apr 29 03:56:17 2022

@author: wwww
"""

import onfly.repo;
onfly.repo.pkg('mp');
import mp.test_mp as tt;


from multiprocessing import Process
print('Ok')
tt.start()
#onfly.repo.urlset('http://192.168.0.110:8000/')


