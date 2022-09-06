import sys,os
import pickle

__pade_exp_poles_res_half__,__pade_exp_poles_res__=None,None;

def load_epd(fhalf=False):
    
    global __pade_exp_poles_res_half__,__pade_exp_poles_res__
    
    if fhalf:
        if not __pade_exp_poles_res_half__ is None:
            return __pade_exp_poles_res_half__
    else:
        if not __pade_exp_poles_res__ is None:
            return __pade_exp_poles_res__
    

    s='/../pade_exp_poles_res_half.pickle' if fhalf else '/../pade_exp_poles_res.pickle'
    
    with open(os.path.abspath(__file__+s), 'rb') as f:
        pr=pickle.load(f,encoding='latin1') if sys.version_info.major==3 else pickle.load(f)

    if fhalf:
        __pade_exp_poles_res_half__ =pr
    else:
        __pade_exp_poles_res__ =pr

    return pr;
    
def poles_res(n,m,fhalf=False):
    return load_epd(fhalf)[m][n]
    
def get(nm,fhalf=False):
    return load_epd(fhalf)[nm[1]][nm[0]]

def get_poles_res_array(nm,fhalf=False):
    import numpy as np
    pr=get(nm,fhalf);
    M=len(pr);
    poles=np.empty((M,),dtype=complex);
    res=np.empty((M,),dtype=complex);
    
    for m in range(M):
        poles[m]=pr[m][0]
        res[m]=pr[m][1]
        
    return (poles,res)