#
import sys,os
if sys.platform=='win32':
    
    import ctypes
    _setenv=ctypes.windll.kernel32.SetEnvironmentVariableW;
    _setenv.rstype=ctypes.c_uint32
    _setenv.argtypes = (ctypes.c_wchar_p,ctypes.c_wchar_p)
    
    
    
    _getenv=ctypes.windll.kernel32.GetEnvironmentVariableW;
    _getenv.rstype=ctypes.c_uint32
    _getenv.argtypes = (ctypes.c_wchar_p,ctypes.c_wchar_p,ctypes.c_uint32)
    
    def isenv(name):
        return not not _getenv(name,(ctypes.c_wchar*8)(),0);
    
    
        
        
    def igetenv(name,dflt=''):
        N=_getenv(name,(ctypes.c_wchar*8)(),0)
        if N>0:
            buf=(ctypes.c_wchar*N)()
            if _getenv(name,buf,2*N):
                return (N,buf.value);
        else:
            return (0,dflt);
        
    def getenv(name,dflt=''):        
        return igetenv(name,dflt=dflt)[1];
    
    
            
    
    def setenv(name,value):
        os.environ[name]=value;
        return _setenv(name,value);
    
else:    
    getenv=os.getenv
    setenv=os.setenv
    
    def isenv(name):
        try:
            q=os.environ[name]            
        except Exception:
            return False;
        return True;
    
    
pid=os.getpid
    
