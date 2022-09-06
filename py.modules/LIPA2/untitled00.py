import pickle
import multiprocessing
from jsobj import *
from tric import *
'''
class Tricky:
    def __init__(self,x):
        self.data=x

    def __setstate__(self, d):
        print('setstate happening')
        self.data = 10

    def __getstate__(self):
        return self.data
        print('getstate happening')
'''
def report(ar,q):
    #ar=Tricky(7)
    #q.put(arg2jso(a=11,b=333))
    ar.data=-7777
    q.put(ar)

if __name__ == '__main__':
    ar = Tricky(5)
    q = multiprocessing.Queue()
    
    q.put(ar)
    print(q.get().data)
    
    p = multiprocessing.Process(target=report, args=(ar, q))
    print('now starting process')
    p.start()
    print('now joining process')
    p.join()
    print('now getting results from queue')
    print(q.get().data)
    print('now getting pickle dumps')
    print(pickle.loads(pickle.dumps(ar)).data)   
