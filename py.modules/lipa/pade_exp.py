#from  ltx.mm  import ltx_js
#import win32api
import os
from ctypes import *
import numpy as np
import pickle

#_tmp=win32api.GetModuleFileName(ltx_js._handle)+'\\..\\..\\..\\tp\\mkl\\bin\\exp_pade_data.dll'
_tmp=os.environ['ltx_js_root']+'tp/mkl/bin/exp_pade_data'
_epdlib=cdll.LoadLibrary(_tmp)

_max_poles_num=_epdlib.max_poles_num
_max_poles_num.restype = c_int

_pade_poles_res=_epdlib.pade_poles_res
_pade_poles_res.restype = c_int
_pade_poles_res.argtypes = (c_int32,c_int32,c_void_p,c_void_p)

_pade_poles_res_half=_epdlib.pade_poles_res_half
_pade_poles_res_half.restype = c_int
_pade_poles_res_half.argtypes = (c_int32,c_int32,c_void_p,c_void_p)

def load_epd(fhalf=False):
    s='/../pade_exp_poles_res_half.pickle' if fhalf else '/../pade_exp_poles_res.pickle'
    with open(__file__+s, 'rb') as f:
        pr=pickle.load(f)
    return pr;

def get_poles_res(n,m,fhalf=False):
    return load_epd(fhalf)[m][n]


class exp_pr:
    def __init__(self,n,m,fhalf=False,dt=1):
        mc=2*m
        p=(c_double*mc)();
        r=(c_double*mc)();
        if fhalf:
            c=_pade_poles_res_half(c_int32(n),c_int32(m),byref(p),byref(r));
        else:
            c=_pade_poles_res(c_int32(n),c_int32(m),byref(p),byref(r));

        pr=();
        poles = np.zeros(c,dtype=np.complex);
        res   = np.zeros(c,dtype=np.complex);
        bt=1.0/dt;
        for k in range(c):
            poles[k]=bt*complex(p[2*k],p[2*k+1]);
            res[k]=bt*complex(r[2*k],r[2*k+1]);
            pr+=((poles[k],res[k]),);
        self.poles=poles;
        self.res=res;
        self.count=c;
        self.dt=dt;
        self.poles_res=self.pr=pr;
        self.fhalf=fhalf;



class poles_res:
    def __init__(self,n,m,fhalf=False,dt=1):
        mc=2*m
        p=(c_double*mc)();
        r=(c_double*mc)();
        if fhalf:
            c=_pade_poles_res_half(c_int32(n),c_int32(m),byref(p),byref(r));
        else:
            c=_pade_poles_res(c_int32(n),c_int32(m),byref(p),byref(r));

        pr=();
        poles = np.zeros(c,dtype=np.complex);
        res   = np.zeros(c,dtype=np.complex);
        bt=1.0/dt;
        for k in range(c):
            poles[k]=bt*complex(p[2*k],p[2*k+1]);
            res[k]=bt*complex(r[2*k],r[2*k+1]);
            pr+=((poles[k],res[k]),);
        self.poles=poles;
        self.res=res;
        self.count=c;
        self.dt=dt;
        self.poles_res=self.pr=pr;
        self.fhalf=fhalf;
#    def pr(self):
#        return (self.poles,self.res)

