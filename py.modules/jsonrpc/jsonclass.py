import json
import sys
import copy
import numbers
from gzip import zlib
from jsobj import *


class local_marshal_classes:
    def __init__(self):
#        self.compress_mode={'mode':0,'level':6}
#        self.compressor=zlib
        self.ipc_marshal={'enabled':False,'marshal':None,'unmarshal':lambda data,flink :None,'cache':lambda:None}
        self.marshalers={};
        self.unmarshalers={};
    
        
    def register(self,cls,jsonclass=None,marshaller=None,unmarshaller=None):
        if jsonclass is None:
            jsonclass=cls.__name__
        self.marshalers[cls]=(jsonclass,marshaller)
        self.unmarshalers[jsonclass]=(cls,unmarshaller)
    

class marshal_classes:
    classes=local_marshal_classes()
    compress_mode={'mode':0,'level':6}
    compressor=zlib
    marshal_error_handler=None;
    unmarshal_error_handler=None;
    
    _instance = None
    @classmethod
    def instance(cls):
        if not cls._instance:
            cls._instance = cls()
        return cls._instance

config=marshal_classes.instance()


if sys.version_info.major==3:
    xrange=range
    unicode=str
    long=int


iter_types = [
    dict,
    list,
    tuple
]

string_types = [
    bytes,
    unicode
]

numeric_types = [
    int,
    long,
    float
]

value_types = [
    bool,
    type(None)
]

def ipc_mode(mode=None):
    m=config.classes.ipc_marshal['enabled'];
    if not mode is None:
        config.classes.ipc_marshal['enabled']=bool(mode);        
    return m;

def ipc_cache_clear():
    return config.classes.ipc_marshal['cache']();    

def compress_mode(mode=-1,level=-1):
    if mode>=0:
        config.compress_mode['mode']=mode
    if level>=0:
        config.compress_mode['level']=level
    return copy.copy(config.compress_mode);

        

def _encode(obj):
    
    obj_type = type(obj)
    
    if isinstance(obj, numbers.Number):
        
        if isinstance(obj, numbers.Integral):
            return int(obj);
        elif isinstance(obj, numbers.Real):
            return float(obj);
        elif isinstance(obj, numbers.Complex):
            obj=complex(obj);            
    
    else:
        
        if obj_type==jsobject:
            return _encode(obj.__dict__)
        
        #if obj_type in numeric_types+string_types+value_types:
        if obj_type in string_types+value_types:   
            return obj
        if obj_type in iter_types:
            if obj_type in (list, tuple):
                new_obj = []
                for item in obj:
                    new_obj.append(_encode(item))
                if obj_type is tuple:
                    new_obj = tuple(new_obj)
                return new_obj
            # It's a dict...
            else:
                new_obj = {}
                for key, value in list(obj.items()):
                    new_obj[key] = _encode(value)
                return new_obj
        
    try:
        (jsonclass,marshaller)=config.classes.marshalers[obj_type];
    except Exception:
        return obj
    return marshaller(obj,jsonclass);
    



def _decode(obj):
    if type(obj) in string_types+numeric_types+value_types:
        return obj
    if type(obj) is list:
        return_list = []
        for entry in obj:
            return_list.append(_decode(entry))
        return return_list
    # Othewise, it's a dict type
    if '__jsonclass__' not in list(obj.keys()):
        return_dict = {}
        for key, value in list(obj.items()):
            new_value = _decode(value)
            return_dict[key] = new_value
        return return_dict
    # It's a dict, and it's a __jsonclass__
    p=obj['__jsonclass__']
    try:
        (cls,unmarshaller)=config.classes.unmarshalers[p[0]];
    except:
        error_handler=config.unmarshal_error_handler;
        if error_handler:
           unmarshaller=error_handler;
        else:
            raise Exception('__jsonclass__: '+p[0]+' is not registered')
            
        
        
    return unmarshaller(p);

    

def jsobject_decode(obj={}):
    if type(obj) in (list, tuple):
        return_list = []
        for entry in obj:
            return_list.append(jsobject_decode(entry))
        return return_list
    if not type(obj) ==dict:
        return obj;
    return_dict = {}
    for key, value in list(obj.items()):
        new_value = jsobject_decode(value)
        return_dict[key] = new_value
    return jsobject(return_dict)

def _jsclass_decode(obj):
    return jsobject_decode(_decode(obj))


def encode_ipc(obj,f_ipc=None):
    old=ipc_mode(f_ipc);
    try:
        t=_encode(obj)
        return json.dumps(t,indent=4, separators=(',', ':'));
    finally:
        ipc_mode(old);
    
def encode(obj):             
    t=_encode(obj)
    return json.dumps(t,indent=4, separators=(',', ':'),default=config.marshal_error_handler)

def decode(s,jslike=False):
    t=json.loads(s);
    return _jsclass_decode(t) if jslike else _decode(t);

def encode_to_file_raw(obj,fn):
    if sys.version_info.major<3:
        open(fn,'wb').write(encode(obj))
    else:
        open(fn,'wb').write(encode(obj).encode('utf8'))
    
def encode_to_file(obj,fn):
    i=ipc_mode(False);
    try:
        encode_to_file_raw(obj,fn);
    finally:
        ipc_mode(i);
    
def decode_from_file(fn,jslike=False):
    if sys.version_info.major<3:
        return decode(open(fn,'rb').read(),jslike)
    else:
        return decode(open(fn,'rb').read().decode('utf8'),jslike)

jslike=jsobject_decode
class_decode=_decode
class_encode=_encode
jso=jsobject




"""
class marshal:
    def __init__(self):
        self.marshalers={}
        self.unmarshalers={}
    

    
    def register(self,cls,jsonclass=cls.__name__,marshaller=None,unmarshaller=None):
        self.marshalers[cls]=(jsonclass,marshaller)
        self.unmarshalers[jsonclass]=(jsonclass,unmarshaller)
    
"""