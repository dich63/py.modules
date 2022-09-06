
import sys    

try:
   from colorama import init
   init()
   #from colorama import Fore, Back, Style
except Exception:
    pass

def sprintf(*args):
    
    f=args[0];
    a=args[1:];
    s=f% tuple(a);
    return s;

def printf(*args):      
    sys.stdout.write(sprintf(*args));

