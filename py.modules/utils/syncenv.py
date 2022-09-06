#
import sys,os
if sys.platform=='win32':
    from ctypes import *
    import ctypes
    _setenv=ctypes.windll.kernel32.SetEnvironmentVariableW;
    _setenv.rstype=ctypes.c_uint32
    _setenv.argtypes = (ctypes.c_wchar_p,ctypes.c_wchar_p)
    
    
    
    _getenv=ctypes.windll.kernel32.GetEnvironmentVariableW;
    _getenv.rstype=ctypes.c_uint32
    _getenv.argtypes = (ctypes.c_wchar_p,ctypes.c_wchar_p,ctypes.c_uint32)
    
    def getenv(name,dflt=''):
        buf=(ctypes.c_wchar*4096)()
        if _getenv(name,buf,4095):
            return buf.value;
        else:
            return dflt;
    
    def setenv(name,value):
        return _setenv(name,value);
    
else:    
    getenv=os.getenv
    setenv=os.setenv
    
pid=os.getpid
    
