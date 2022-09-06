# -*- coding: utf-8 -*-
import sys
import pydoc

def help2file(filepath, request):
    f = open(filepath, 'w')
    try:
        sys.stdout = f
        sys.stderr = f
        pydoc.help(request)
    except Exception:
        pass
        
    f.close()
    sys.stdout = sys.__stdout__
    sys.stderr = sys.__stderr__
    return
