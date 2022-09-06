# -*- coding: utf-8 -*-
"""
Created on Tue Apr 12 19:01:58 2022

@author: wwww
"""

import pickle
import multiprocessing

class Tricky:
    def __init__(self,x):
        self.data=x

    def __setstate__(self, d):
        print('setstate happening')
        self.data = 10

    def __getstate__(self):
        return self.data
        print('getstate happening')

def report(ar,q):
    q.put(ar.data)

if __name__ == '__main__':
    
    ar = Tricky(5)
    
    q = multiprocessing.JoinableQueue()
    p = multiprocessing.Process(target=report, args=(ar, q))
    print('now starting process')
    p.start()
    print('now joining process')
    p.join()
    print('now getting results from queue')
    print(q.get())
    print('now getting pickle dumps')
    print(pickle.loads(pickle.dumps(ar)).data)   