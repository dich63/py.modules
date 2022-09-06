import os;
print ('pid='+str(os.getpid()))

__binpath__=os.path.abspath(__file__+'/../../bin')+'/'
print (__binpath__)

def __set_lib_env__(en,ln):
    global __binpath__;
    os.environ[en]=__binpath__+ln;
    
__set_lib_env__('context_wrapper',"ictxwrr.dll")
__set_lib_env__('parallel_wrapper_lib',"tbb_wrapper.dll")
__set_lib_env__('sparse_context_lib','tpt_klu.dll')
    
    



if 1:
     os.environ['context_lib']="V:/Projects/tbb_wrapper/x64/Debug/ictxwrr.dll"     
     os.environ['invoke_context_lib']="V:/Projects/tbb_wrapper/x64/Debug/invoke_context.dll"
     #os.environ['sparse_context_lib']='O:\\__ss\\suitesparse-metis-for-windows\\bin\\x64\\Debug\\tpt_klu.dll'



#os.environ['context_wrapper']="V:/Projects/tbb_wrapper/x64/release/ictxwrr.dll"

#os.environ['context_wrapper']="V:/Projects/tbb_wrapper/x64/Release/context_holder.dll"
#os.environ['context_wrapper']="V:/Projects/tbb_wrapper/x64/Debug/ictxwrr.dll"
#os.environ['context_wrapper']="c:/temp/0/ictxwrr.dll"  
#os.environ['context_wrapper']="c:/temp/0/xz.a"  
#os.environ['context_wrapper']="V:/Projects/tbb_wrapper/x64/release/ictxwrr.dll"  
#os.environ['context_wrapper']="V:/Projects/tbb_wrapper/x64/debug/context_holder.dll"    
print ('context_wrapper='+os.environ['context_wrapper'])
print ('invoke_context_lib='+os.environ['invoke_context_lib'])
#print ('sparse_context_lib='+os.environ['sparse_context_lib'])
#print ('MID_lib='+os.environ['MID_lib'])