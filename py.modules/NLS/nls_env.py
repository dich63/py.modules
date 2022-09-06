# 
import os;
import utils.process as ps

'''
import ctypes

_setenv=ctypes.windll.kernel32.SetEnvironmentVariableW;
_setenv.rstype=ctypes.c_uint32
_setenv.argtypes = (ctypes.c_wchar_p,ctypes.c_wchar_p)
'''


__binpath__=os.path.abspath(__file__+'/../bin')+'/'

ps.setenv('context_wrapper',__binpath__+'ictx/ictxwrr.dll')
ps.setenv('nls_pade_lib',__binpath__+'nls_pade.dll')
ps.setenv('nls_pade_lib.path',__binpath__+'tbb/')
ps.setenv('ssf_lib',__binpath__+'ssf.dll')

'''
os.environ['context_wrapper']=__binpath__+'ictx/ictxwrr.dll'

os.environ['nls_pade_lib']=__binpath__+'nls_pade.dll'
os.environ['nls_pade_lib.path']=__binpath__+'tbb/'
os.environ['ssf_lib']=__binpath__+'ssf.dll'
'''
#_setenv('context_wrapper',os.environ['context_wrapper'])



_dbg=0
_dump=1
if _dbg: 
    #os.environ['nls_pade_lib']='V:/HW/nls_pade/x64/Debug/nls_pade.dll'
    ps.setenv('nls_pade_lib','V:/HW/nls_pade_op/x64/Debug/nls_pade.dll')
    ps.setenv('nls_pade_lib.path','D:/Programs/Intel/2020_09_17_02-23-24/')
    
    '''
    os.environ['nls_pade_lib']='V:/HW/nls_pade_op/x64/Debug/nls_pade.dll'
    os.environ['nls_pade_lib.path']='D:/Programs/Intel/2020_09_17_02-23-24/'
    '''



if _dump: 
    print ('pid=',os.getpid())
    print ('context_wrapper='+os.environ['context_wrapper'])
    print ('nls_pade_lib='+os.environ['nls_pade_lib'])
    print ('ssf_lib='+os.environ['ssf_lib'])


def load_lib_with_path1(name):
    
    import ctypes 
    
    path=spath=ps.getenv('path');
    
    #htbb=ctypes.cdll.LoadLibrary(__binpath__+'tbb/tbb.dll')
    fm=False;
    try:
        
        [ok,fn]=ps.igetenv(name)
        if ok:
            path=fn+'/../;'+path;
            fm=True;       
            
        
        [ok,np]=ps.igetenv(name+'.path')
        if ok:
            path=np+';'+path;
            fm=True;
            
        ps.setenv('path',path)
        ps.setenv('lib_path',path)       
        
        return ctypes.cdll.LoadLibrary(fn);
                    
    finally:
        if fm:
            ps.setenv('path',spath);            
            

         
def load_lib_with_path2(name):  
    
    import ctypes     
    
    [ok,fn]=ps.igetenv(name)
    if ok:
        cookie=os.add_dll_directory(fn+'\\..\\');       
        
    [ok,np]=ps.igetenv(name+'.path')    
    if ok:
        cookie2=os.add_dll_directory(np);
    
    return ctypes.cdll.LoadLibrary(fn);
    

def load_lib_with_path(name):
    try:
        return load_lib_with_path2(name)
    except Exception:
        return load_lib_with_path1(name)

    