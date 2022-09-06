import sys

from utils.p23 import *

'''
to_bytes=bytes

if sys.version_info.major==3: #moronic python3
    xrange=range
    
    def unicode(s):
        if type(s)==bytes:
            return s.decode('utf8');
        else:
            return s
        
    def to_bytes(s):
        if type(s)==bytes:
            return s;
        else:
            return bytes(str(s),'utf8')


'''

