# -*- coding: utf-8 -*-

from jsonrpc.jsonclass import config,class_decode,class_encode,json

class JSON_marshaling_error(Exception):
    def __init__(self,message=""):
        self.message=message;
    def __str__(self):
        return "JSON_marshaling_error: "+self.message;
    def __repr__(self):
        return self.__str__();
    
def unmarshal_error(d):
    msg=d[1][0].get('message','?unknown?');
    return JSON_marshaling_error(msg);

def marshal_error(v,json_class_name='JSON_marshaling_error'):    
    s='Object of type '+str(type(v))+' is not JSON serializable';
    r={'__jsonclass__':[json_class_name,[{'message':s}]]}
    return r;
    #return json.dumps(r);
    


config.unmarshal_error_handler=lambda p : JSON_marshaling_error('__jsonclass__: '+p[0]+' is not registered')
config.marshal_error_handler=lambda v : marshal_error(v);
config.classes.register(JSON_marshaling_error,'JSON_marshaling_error',marshal_error,unmarshal_error)
   