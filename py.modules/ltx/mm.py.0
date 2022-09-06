#
import sys
import os     
import win32com.client
import win32com.client.dynamic
import pythoncom
from ctypes import *
from ctypes.wintypes import BOOL
#from comtypes import IUnknown
#from comtypes.automation import IDispatch
from win32com.client import GetObject
import numpy as np

class  OLE_VAR_U(Union):
    _fields_=[("llVal",c_longlong),("lVal",c_long),
              ("dblVal",c_double),("fltVal",c_float),
              ("bstrVal",c_wchar_p),("byref",c_void_p)]


class OLE_VARIANT(Structure):
    _fields_ = [("vt", c_ushort),("rsvd__",c_ushort*3),
    ("V",OLE_VAR_U)]

_ltx_error_raiser_PT=CFUNCTYPE(c_int,c_wchar_p,c_wchar_p)

class callback_context_arguments_t(Structure):
    _fields_ = [("args", c_void_p),("error",c_void_p)]


#_ltx_on_callback_PT=CFUNCTYPE(c_int,POINTER(OLE_VARIANT),py_object,POINTER(callback_context_arguments_t))
_ltx_on_callback_PT=CFUNCTYPE(c_int,POINTER(OLE_VARIANT),py_object,c_void_p)
_ltx_on_exit_PT=CFUNCTYPE(c_int,py_object)

_PyComLib = PyDLL(pythoncom.__file__)

_PyCom_VariantFromPyObject =_PyComLib.PyCom_VariantFromPyObject
_PyCom_VariantFromPyObject.restype = c_int
_PyCom_VariantFromPyObject.argtypes = (py_object,POINTER(OLE_VARIANT))


_PyCom_PyObjectFromIUnknown = _PyComLib.PyCom_PyObjectFromIUnknown
_PyCom_PyObjectFromIUnknown.restype = py_object
_PyCom_PyObjectFromIUnknown.argtypes = (c_void_p, c_void_p, BOOL)

_PyCom_InterfaceFromPyObject = _PyComLib.PyCom_InterfaceFromPyObject
_PyCom_InterfaceFromPyObject.restype = c_int32
_PyCom_InterfaceFromPyObject.argtypes = (py_object, c_void_p,c_void_p, BOOL)







#eo=GetObject('ltx.bind:external')
bindObject = GetObject('ltx.bind:')

ltx_js = cdll.ltx_js




com_apartment_type=ltx_js.ltx_apartment_type
com_apartment_type.restype = c_int32
com_apartment_type.argtypes = ()


get_tid=windll.kernel32.GetCurrentThreadId
get_tid.restype = c_int32
get_tid.argtypes = ()


_ltx_raise_exception=ltx_js.ltx_raise_exception
_ltx_raise_exception.restype = c_int32
_ltx_raise_exception.argtypes = (c_wchar_p,c_wchar_p)

_ltx_fill_exception=ltx_js.ltx_fill_exception
_ltx_fill_exception.restype = c_int32
_ltx_fill_exception.argtypes = (c_void_p,c_wchar_p,c_wchar_p)


_ltx_callback=ltx_js.ltx_callback
_ltx_callback.restype = c_int32
_ltx_callback.argtypes = (c_void_p,_ltx_on_callback_PT,py_object,c_void_p)

_ltx_process_callback_loop=ltx_js.ltx_process_callback_loop
_ltx_process_callback_loop.restype = c_int32
_ltx_process_callback_loop.argtypes = (_ltx_on_callback_PT,py_object,c_void_p)



_ltx_mm_DataInfo = ltx_js.ltx_mm_DataInfo
_ltx_mm_DataInfo.restype = c_int32
_ltx_mm_DataInfo.argtypes = (c_void_p, c_void_p,c_void_p,c_void_p)
_ltx_mm_CommitRegionPtr = ltx_js.ltx_mm_CommitRegionPtr
_ltx_mm_CommitRegionPtr.restype = c_int32
_ltx_mm_CommitRegionPtr.argtypes = (c_void_p, c_int64,c_int64,c_void_p)
_ltx_mm_DecommitRegionPtr = ltx_js.ltx_mm_DecommitRegionPtr
_ltx_mm_DecommitRegionPtr.restype = c_int32
_ltx_mm_DecommitRegionPtr.argtypes = (c_void_p,c_void_p)



_ltx_mm_Lock = ltx_js.ltx_mm_Lock
_ltx_mm_Lock.restype = c_int32
_ltx_mm_Lock.argtypes = (c_void_p,)


_ltx_mm_Unlock = ltx_js.ltx_mm_Unlock
_ltx_mm_Unlock.restype = c_int32
_ltx_mm_Unlock.argtypes = (c_void_p,)


def __t(s):
    return s;

if sys.version_info.major==3:
    xrange=range
    unicode=__t

"""
"""
class external_callbacks_t(Structure):
    _fields_ = [
    ("name",c_char_p),
    ("proc",c_void_p),
    ("attr",c_int)
    ];


long_pvoid_p=CFUNCTYPE(c_long,c_void_p);
long_pvoid_pvoid_p=CFUNCTYPE(c_long,c_void_p,c_void_p);
void_pvoid_p=CFUNCTYPE(None,c_void_p);
long_pvoid_p_proc=CFUNCTYPE(c_long,c_void_p,void_pvoid_p);

long_pvoid_p_proc_pvoid_p=CFUNCTYPE(c_long,c_void_p,void_pvoid_p,c_void_p);


class context_holder_t(Structure):
    _fields_ = [    
    ("addref",long_pvoid_p),
    ("release",long_pvoid_p),
    ("unwrap_context",long_pvoid_pvoid_p),
    ("wrap_context",long_pvoid_p_proc_pvoid_p),
    ("link",long_pvoid_pvoid_p),
    ("clear_links",long_pvoid_p)    
    ];

pcontext_holder_t=POINTER(context_holder_t);
_get_context_utils=ltx_js.get_context_utils
_get_context_utils.restype =c_void_p

context_utils=context_holder_t.from_address(_get_context_utils())
NullFunct=void_pvoid_p(0)
"""
"""
class GUID(Structure):
    _fields_ = [("c",c_char * 16)]

pc = pythoncom
vtypes = {pc.VT_R8:c_double,pc.VT_R4:c_float,pc.VT_I1:c_int8,pc.VT_UI1:c_uint8,pc.VT_I2:c_int16,pc.VT_UI2:c_uint16,pc.VT_BOOL:c_uint16,pc.VT_I4:c_int32,pc.VT_UI4:c_uint32,pc.VT_I8:c_int64,pc.VT_UI8:c_uint64,pc.VT_I8:c_int64}
del pc

def com_check(hr):
    if hr != 0:
        raise Exception('COM Error:',hr)
    return hr

def test_mm_buffer(m):
    for i in xrange(m.count):
        m[i] = i

def test_mm_buffer2(m):
    v = m.value
    for i in xrange(m.count):
        v[i] = i

def GUID_from_string(s='{00000000-0000-0000-C000-000000000046}'):
    iid = GUID();
    s=unicode(s);
    com_check(hr=windll.ole32.CLSIDFromString(s,byref(iid)))
    return iid;

"""
def comtypes2pywin(ptr, interface=None):
    
    if interface is None:
        interface = IUnknown
    return _PyCom_PyObjectFromIUnknown(ptr, byref(interface._iid_), True)
    """

def comtypes2pywin(ptr, interface='{00000000-0000-0000-C000-000000000046}'):
    #print('1.1.0')
    iid=GUID_from_string(interface)
    #print('1.1.1')
    return _PyCom_PyObjectFromIUnknown(ptr, byref(iid), True)

#'{00020400-0000-0000-C000-000000000046}'

def comtypes2Dispatch(ptr):
    pd=comtypes2pywin(ptr,'{00020400-0000-0000-C000-000000000046}')
    #print('1.1')
    disp=win32com.client.dynamic.Dispatch(pd)
    return disp

def _ltx_disp_method(result,funct,callback_context):
    
    #print('start')
    
    try:
        
        #print('fun=',funct)
        cca=(callback_context_arguments_t).from_address(callback_context)
        #print('1')
        arguments=comtypes2Dispatch(cca.args);
        #print('2')
        r=funct(arguments)
        _PyCom_VariantFromPyObject(r,result)
    except Exception as e:
        _ltx_fill_exception(callback_context,repr(e),sys.version)
    
    #print('tid=',get_tid(),'CC:: ',funct,callback_context)
    return 0

_ltx_wrap_dc=_ltx_on_callback_PT(_ltx_disp_method)

def ltx_create_callback(funct):
     pp=c_void_p(0);
     np=c_void_p(0);
     #print('main-tid=',get_tid())
     hr=_ltx_callback(byref(pp),_ltx_wrap_dc,funct,np)
     if hr==0:
         #pd=comtypes2pywin(pp,IDispatch)
         #disp=win32com.client.dynamic.Dispatch(pd)
         return comtypes2Dispatch(pp)
     else:
        return None


def _ltx_disp_method2(result,funct,callback_context):
    
    try:
        #print('start')
        cca=(callback_context_arguments_t).from_address(callback_context)
        #print('1')
        arguments=comtypes2Dispatch(cca.args);
        #print('2')
        r=funct(*arguments)
        _PyCom_VariantFromPyObject(r,result)
    except Exception as e:
        _ltx_fill_exception(callback_context,repr(e),sys.version)
    
    #print('tid=',get_tid(),'CC:: ',funct,callback_context)
    return 0

_ltx_wrap_dc=_ltx_on_callback_PT(_ltx_disp_method)

def ltx_create_callback2(funct):
     pp=c_void_p(0);
     np=c_void_p(0);
     #print('main-tid=',get_tid())
     hr=_ltx_callback(byref(pp),_ltx_wrap_dc,funct,np)
     if hr==0:
         #pd=comtypes2pywin(pp,IDispatch)
         #disp=win32com.client.dynamic.Dispatch(pd)
         return comtypes2Dispatch(pp)
     else:
        return None


def ltx_create_process_callback_loop(funct):
     np=c_void_p(0);
     hr=_ltx_process_callback_loop(_ltx_wrap_dc,funct,np)
     return hr




class mm_buffer:
    def __init__(self,o):
        self.obj = o
        u = unicode('{C6D197F2-609B-4514-8315-BF6810CA6C1C}')
        iid = GUID()
        com_check(hr=windll.ole32.CLSIDFromString(u,byref(iid)))
        pd = c_void_p()
        hr = _PyCom_InterfaceFromPyObject(o._oleobj_,byref(iid),byref(pd),0)
        if hr == 0:
            com_check(0x80070057)
        self.pobj = pd
        count = c_int64()
        elsize = c_int32()
        vt = c_uint16()
        hr = _ltx_mm_DataInfo(pd,byref(vt),byref(count),byref(elsize))
        com_check(hr)
        self.count = count.value
        self.elsize = elsize.value
        ll = c_int64(count.value * elsize.value)
        ptr = c_void_p()
        hr = _ltx_mm_CommitRegionPtr(pd,0,ll,pointer(ptr))
        com_check(hr)
        self.ptr = ptr
        cty = vtypes[vt.value]
        self.type = cty
        NN = count.value
        self.__value = (cty * NN).from_address(ptr.value)
        self.__value.__holder__ =o;
    
    def toarray(self,shape=None,dtype=None,copy=True,count=-1,offset=0):
        import numpy as np
        dtype = dtype if dtype else self.type;
        v=np.frombuffer(self.value,dtype,count,offset);
        if copy:
            v=v.copy();
        else:
            #v.base.__holder__=self.__value.__holder__;
            v.base.__holder__=self;
        if shape:
            v=v.reshape(shape);
        return v
    
    @property
    def value(self):
        return  self.__value

    def detach(self):        
        _ltx_mm_DecommitRegionPtr(self.pobj,self.ptr)
        #print('_ltx_mm_DecommitRegionPtr')
        self.__value = None
        return self.obj
    def lock(self):
        _ltx_mm_Lock(self.pobj);
        
    def unlock(self):
        _ltx_mm_Unlock(self.pobj);

    def objref_detach(self):
        o = self.detach()
        self.stub = bindObject('stub')
        return self.stub(o)

    def __del__(self):
        self.detach()

    def __getitem__(self,k):
        return self.__value[k]

    def __setitem__(self,k,v):
         self.__value[k] = v

    def __len__(self):
        return int(self.count)


def objref2array(objref,flink=False):
    o=bindObject(objref);
    b=mm_buffer(o);
    return b.toarray(copy= not flink);

def lock_array(f):
    if '__holder__' in f.base.__dict__:
        f.base.lock();
        return True;
    return False;
        
def unlock_array(f):
    if '__holder__' in f.base.__dict__:
        f.base.unlock();
        return True;
    return False;
    
    

def create_mm_buffer(count,type='double',options=''):
    b = bindObject('mm_buffer:length=#1;type=#2;'+options,count,type)
    return mm_buffer(b)

def mm_buffer_from_array(a,locked=True):
    r=np.ravel(a)
    dtype=r.dtype;
    if dtype==bool:
        r=np.array(r,dtype=np.int16);
    buffer=create_mm_buffer(r.size,type=np.dtype(dtype).name)
    rd=buffer.toarray(copy=False)
    np.copyto(rd,r)
    if not locked:
        buffer.unlock();
    return buffer

    
    
    
def  jsvm(opts=''):
    s='script:'+opts;
    return bindObject(s);

def external():
    return bindObject('external');

_jsg_=[False];

 


def wrap_call_method(o,name):
#     if gg:
#        _jsg_=jsvm() 
#     gg=_jsg_
     if not _jsg_[0] :
         _jsg_[0]=jsvm('debug=0')
         _jsg_[0]('require("utils",global)')

     return  _jsg_[0]('(function(o,n){ var obj=o,name=n;\n ; return function (){ debugger; return apply_function(obj,arguments,name) }; })(arguments[0],arguments[1])',o,name)

