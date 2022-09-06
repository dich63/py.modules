# -*- coding: utf-8 -*-
"""
Created on Fri Apr 29 02:19:31 2022

@author: wwww
"""

import onfly.repo


from multiprocessing import Process


onfly.repo.urlset('http://192.168.0.110:8000/')

onfly.repo.pkg('mp')



import os
import mp.pp


def print_func(continent='Asia'):
    import mp.c.A.B.c
    print('['+str(os.getpid())+'] The name of continent is : ', continent)

if __name__ == "__main__":  # confirms that the code is under main function
    import mp.c.A.c
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