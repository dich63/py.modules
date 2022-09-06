#
from ltx.external import externObj
import numpy as np    

external=externObj()

a=external[0]
b=external[1]   
c=np.array(a)+np.array(b);
external.result=c

