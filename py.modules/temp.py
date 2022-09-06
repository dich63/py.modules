# -*- coding: utf-8 -*-
"""
Редактор Spyder

Это временный скриптовый файл.
"""


import asyn.SharedWorker as sw
tic=sw.Tic()
import ltx.ipc_marshal
import jsonrpc.ndarray_marshal
import jsonrpc.jsonclass as jc


import ltx.mm as mm
mm.create_mm_buffer
mm.create_mm_buffer(100)
import numpy as np


#a=np.random.rand(int(1024/8),1024,1024)
N=1
#a=np.random.rand(int(1024/8),1024,1024)
a=np.random.rand(int(1024/8)*1024*1024)
#jc.ipc_cache_clear()
jc.ipc_mode(0)

tic.start()
for k in range(N):
    s=jc.encode(a);
    ac=jc.decode(s)
    #jc.ipc_cache_clear()
    
print('mode:',jc.ipc_mode() ,  'tic:',tic.sec()*1.e3/N,' Ms len(s)=',len(s)/1e6,'MB')
    
jc.ipc_mode(1)
for k in range(N):
    tic.start()
    s=jc.encode(a);
    ac=jc.decode(s)
    #jc.ipc_cache_clear()
    t=tic.sec()
print('mode:',jc.ipc_mode() ,  'tic:',t*1.e3/N,' Ms len(s)=',len(s)/1e3,'KB')
    
clipbrd=mm.bindObject('module: lib=**;proc=ltx_callback_list::clipbrd;flags.ftm=1')
clipbrd(s)
#tic.start();jc.ipc_cache_clear();t=tic.sec();print('ipc_cache_clear tic:',t*1.e3/N,' Ms')
'''
tic.start()
N=4
for k in range(N):
    b=mm.mm_buffer_from_array(a)
    s=ltx.ipc_marshal.stubcache(b.obj)
    
print('tic:',tic.sec()*1.e3/N,' Ms')
#print('mon:',s)

a=np.random.rand((4,3))

r=np.ravel(a)
dtype=r.dtype;
mr=mm.create_mm_buffer(r.size,type=np.dtype(dtype).name)
rd=mr.toarray(copy=False)
np.copyto(rd,r)
'''



