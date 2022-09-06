# -*- coding: utf-8 -*-
"""
Created on Fri Apr 29 02:19:31 2022

@author: wwww
"""



import onfly.repo;

from multiprocessing import Process


import os ;print('[ pid='+str(os.getpid())+'] module <<'+__name__+'>> locate: ',__file__)


onfly.repo.pkg('mp')

#import mp




def print_func(continent='Asia'):
    print('start ['+str(os.getpid())+'] The name of continent is : ', continent)
    
    #import mp.c.A.B.c
    print('stop ['+str(os.getpid())+'] The name of continent is : ', continent)

def start():
    #import mp.pp
    from multiprocessing import freeze_support
    freeze_support()
    names = ['America', 'Europe', 'Africa']
    procs = []
    proc = Process(target=print_func)  # instantiating without any argument
    procs.append(proc)
    proc.start()

    # instantiating process with arguments
    for name in names:
        # print(name)
        proc = Process(target=print_func, args=(name,))
        procs.append(proc)
        proc.start()

    # complete the processes
    for proc in procs:
        proc.join()
