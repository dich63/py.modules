#
from utils import set_cb_log

def ltx_event_set_ref(key,ref=None):
    
    if not ref is None:
        from win32com.client import GetObject
        ref=GetObject(ref);
        
    set_cb_log(key,ref);    

    