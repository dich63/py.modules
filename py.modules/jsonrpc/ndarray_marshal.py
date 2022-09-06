# -*- coding: utf-8 -*-
"""
Created on Sat Jan 23 21:39:52 2016

@author: dich
"""
import numpy as np
import base64
from .jsonclass import config,class_decode,class_encode
#from gzip import zlib

"""
np2js={np.int8:'Int8Array',np.int16:'Int16Array',np.int32:'Int32Array'
,np.uint8:'Uint8Array',np.uint16:'Uint16Array',np.uint32:'Uint32Array'
,np.float32:'Float32Array',np.float64:'Float64Array'}

js2np={'Int8Array':np.int8,'Int16Array':np.int16,'Int32Array':np.uint32,
'Uint8Array':np.uint8,'Uint16Array':np.int16,'Uint32Array':np.uint32,
'Float32Array':np.float32,'Float64Array':np.float64}
"""

np2js={'int8':'Int8Array','int16':'Int16Array','int32':'Int32Array'
,'uint8':'Uint8Array','uint16':'Uint16Array','uint32':'Uint32Array'
,'float32':'Float32Array','float64':'Float64Array',
'int64':'Int64Array','uint64':'Uint64Array'}

js2np={'Int8Array':'int8','Int16Array':'int16','Int32Array':'uint32',
'Uint8Array':'uint8','Uint16Array':'int16','Uint32Array':'uint32',
'Float32Array':'float32','Float64Array':'float64'}

def _reverse(l):
    l=list(l);
    l.reverse();
    return l
    
def base64_decompress(s):
    (prfx,s)=s.split(',',1)
    s=base64.b64decode(s)
    fz=prfx.split('base64');
    if fz[1]==':lz4':
        s=config.compressor.decompress(s)
    return s;
         
     
     
#    data=base64.b64decode(o[dn].split(',',1)[1])
    
def ndarray_base64_marshal(a,dn='data',fshape=True):
    atype=np2js[a.dtype.name];
    nc=np.int(np.prod(a.shape))
    ipc=config.classes.ipc_marshal;
    if ipc['enabled'] and (not ipc['marshal'] is None):
        data=ipc['marshal'](a);
    else:
        prfx='data:'+atype.lower()+'_le['+str(nc)+']:base64';
        if config.compress_mode['mode']:
            prfx+=':lz4,'
            ca=config.compressor.compress(a.tobytes(),config.compress_mode['level'])
        else:
            prfx+=','
            ca=np.ascontiguousarray(a);
        data=prfx+base64.b64encode(ca).decode()
        
    r={dn:data,"type":atype,"size":_reverse(a.shape)} if fshape else {dn:data,"type":atype}
    return r



def ndarray_base64_unmarshal(o,dn='data'):
    #data=base64.b64decode(o[dn].split(',',1)[1])      
    data=o[dn];
    #
    r=config.classes.ipc_marshal['unmarshal'](data,False);
    #r=config.classes.ipc_marshal['unmarshal'](data,True);
    if  r is None:
        data=base64_decompress(data);
        dtype=js2np[o['type']];
        r=np.frombuffer(data,dtype=dtype)
        
    shape=o.get('size',None)
    if shape:
        r=r.reshape(_reverse(shape));
    return r



def ndarray_json_marshal(a,json_class_name):
    t=a.dtype.name
    if (t=='complex128') or (t=='complex64'):
        r=complex_json_marshal(a)
    else:
        if t=='bool':
            a=np.array(a,dtype=np.int16);

        r={'__jsonclass__':[json_class_name,[ndarray_base64_marshal(a)]]}
    return r

def ndarray_json_unmarshal(d):
    #o=d['__jsonclass__'][1][0];
    return ndarray_base64_unmarshal(d[1][0])
    

def matlab_cell_json_unmarshal(d):
    return class_decode(d[1][0]['data'])
    

class matlab_CellTensor:
    pass

"""
class json_complex:
    pass
"""

def complex_json_marshal(a,json_class_name='complex'):
    r={'__jsonclass__':[json_class_name,[{'re':class_encode(a.real),'im':class_encode(a.imag)}]]}
    return r

def complex_json_unmarshal(d):
    #print(d)
    o=d[1][0]
    re=class_decode(o['re'])
    im=class_decode(o['im'])
    return re+1j*im



config.classes.register(np.ndarray,'NumericTensor',ndarray_json_marshal,ndarray_json_unmarshal)
config.classes.register(complex,'complex',complex_json_marshal,complex_json_unmarshal)
config.classes.register(np.complex128,'complex',complex_json_marshal,complex_json_unmarshal)
config.classes.register(np.complex64,'complex',complex_json_marshal,complex_json_unmarshal)
config.classes.register(matlab_CellTensor,'CellTensor',None,matlab_cell_json_unmarshal)
