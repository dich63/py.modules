import os;
print ('pid='+str(os.getpid()))

#os.environ['MID_lib']="V:\\pycpp\\MID2D\\x64\\release\\MID2D.dll"
'''
os.environ['context_wrapper']="V:/Projects/tbb_wrapper/x64/Release/context_holder.dll"
os.environ['parallel_wrapper_lib']="V:/Projects/tbb_wrapper/x64/Release/tbb_wrapper.dll"
os.environ['sparse_context_lib']='O:\\__ss\\suitesparse-metis-for-windows\\bin\\x64\\release\\tpt_klu.dll'
'''

#
#os.environ['context_wrapper']="V:/Projects/tbb_wrapper/x64/release/ictxwrr.dll"  
#os.environ['context_wrapper']="V:/Projects/tbb_wrapper/x64/DEBUG/ictxwrr.dll"  
#os.environ['context_wrapper']="V:/ipc/bin/x64/ltx_js.dll"
#os.environ['context_wrapper']="V:/ipc/bin_d/x64/ltx_js.dll"

''''''
__binpath__=os.path.abspath(__file__+'/../bin')+'/'
print (__binpath__)

def __set_lib_env__(en,ln):
    global __binpath__;
    os.environ[en]=__binpath__+ln;
    
__set_lib_env__('context_wrapper',"ictxwrr.dll")
__set_lib_env__('parallel_wrapper_lib',"tbb_wrapper.dll")
__set_lib_env__('sparse_context_lib','tpt_klu.dll')
    

os.environ['nls_wrapper_lib']='??????'
    
'''
os.environ['context_wrapper']=bpath+"ictxwrr.dll"
os.environ['parallel_wrapper_lib']=bpath+"tbb_wrapper.dll"
os.environ['sparse_context_lib']=bpath+'tpt_klu.dll'
'''



if 0:
    '''
    os.environ['context_wrapper']="V:/Projects/tbb_wrapper/x64/release/ictxwrr.dll"
    os.environ['parallel_wrapper_lib']="V:/Projects/tbb_wrapper/x64/release/tbb_wrapper.dll"
    os.environ['sparse_context_lib']='O:\\__ss\\suitesparse-metis-for-windows\\bin\\x64\\release\\tpt_klu.dll'
    '''
    os.environ['context_wrapper']="V:/Projects/tbb_wrapper/x64/Debug/ictxwrr.dll"     
    os.environ['parallel_wrapper_lib']="V:/Projects/tbb_wrapper/x64/Debug/tbb_wrapper.dll"
    os.environ['sparse_context_lib']='O:\\__ss\\suitesparse-metis-for-windows\\bin\\x64\\Debug\\tpt_klu.dll'
    os.environ['nls_wrapper_lib']="V:\\HW\\nls\\x64\\Debug\\nls.dll"     
    
    #'''
    os.environ['context_wrapper']="V:/Projects/tbb_wrapper/x64/release/ictxwrr.dll"     
    os.environ['parallel_wrapper_lib']="V:/Projects/tbb_wrapper/x64/release/tbb_wrapper.dll"
    os.environ['nls_wrapper_lib']="V:\\HW\\nls\\x64\\release\\nls.dll"
    #'''
     
     




#os.environ['context_wrapper']="V:/Projects/tbb_wrapper/x64/release/ictxwrr.dll"

#os.environ['context_wrapper']="V:/Projects/tbb_wrapper/x64/Release/context_holder.dll"
#os.environ['context_wrapper']="V:/Projects/tbb_wrapper/x64/Debug/ictxwrr.dll"
#os.environ['context_wrapper']="c:/temp/0/ictxwrr.dll"  
#os.environ['context_wrapper']="c:/temp/0/xz.a"  
#os.environ['context_wrapper']="V:/Projects/tbb_wrapper/x64/release/ictxwrr.dll"  
#os.environ['context_wrapper']="V:/Projects/tbb_wrapper/x64/debug/context_holder.dll"    
print ('context_wrapper='+os.environ['context_wrapper'])
print ('parallel_wrapper_lib='+os.environ['parallel_wrapper_lib'])
print ('sparse_context_lib='+os.environ['sparse_context_lib'])
print ('nls_wrapper_lib='+os.environ['nls_wrapper_lib'])
#print ('MID_lib='+os.environ['MID_lib'])