from MID.clib import *
from ctypes import *
import numpy as np
import copy
class fbs_params_ids_t(Structure):
    _pack_=8
    _fields_ = [
    ("R_in",c_double),
    ("R_out",c_double),
    ("R_obs",c_double),
    ("R_cut",c_double),
    ("gamma",c_double),
    ("D",c_longlong),  
    ("nvxs",c_longlong),
    ("ntrs",c_longlong),
    ("p_ntrs_cut",c_void_p),
    ("p_vxs",c_void_p),
    ("p_trs",c_void_p),
    ("p_mask_trs",c_void_p),
    ("scale",c_double),
    ("p_morph_scale",c_void_p),
    ("p_mDxD_deform",c_void_p),
    ("p_morph_proc",c_void_p),
    ];

p_fbs_params_ids_t=POINTER(fbs_params_ids_t)
'''
_fuzzy_boundaries_condition_ids=midlib.fuzzy_boundaries_condition_ids
_fuzzy_boundaries_condition_ids.restype=c_int32    
_fuzzy_boundaries_condition_ids.argtype=[p_fbs_params_ids_t]
'''
class fuzzy_bc_t(object):
    def __init__(self,vxs,trs,mask_trs=None,morph_scale=None,mDxD_deform=None,D=2):
        

        szv=np.size(vxs);
        szt=np.size(trs);

        nT=D+1;
        

        nvxs=int(szv/D);
        ntrs=int(szt/nT);
        vxs=np.reshape(vxs,(nvxs,D))    
        trs=np.reshape(trs,(ntrs,nT))    
        
        self.ntrs_cut=np.array([ntrs],dtype=np.int64)
        
        if mask_trs is None:
            mask_trs=np.zeros(ntrs,dtype=np.uint16);
            
        if morph_scale is None:
            morph_scale=np.zeros(ntrs,dtype=np.float64);
            
        self.mask_trs=mask_trs
        self.morph_scale=morph_scale
        self.mDxD_deform=mDxD_deform;
        
        self.vxs,self.trs=vxs,trs
        
        self.D,self.nT,self.nvxs,self.ntrs=D,nT,nvxs,ntrs

    @property        
    def  trs_cut(self):
        return float(self.ntrs_cut[0])/float(self.ntrs);
        
    def make(self,R_in=1,R_out=1,R_obs=1,R_cut=1,gamma=1.0,scale=1.0,mDxD_deform=None):
        
        
        if mDxD_deform is None:
            mDxD_deform=self.mDxD_deform;
            
        self.params=p=fbs_params_ids_t()    
        
        
        
        p.R_in,p.R_out,p.R_obs,p.R_cut,p.gamma,p.scale=R_in,R_out,R_obs,R_cut,gamma,scale        
        
        
        p.D,p.nvxs,p.ntrs=self.D,self.nvxs,self.ntrs
        
        p.p_mask_trs=self.mask_trs.ctypes.data
        p.p_morph_scale=self.morph_scale.ctypes.data
        p.p_vxs=self.vxs.ctypes.data
        p.p_trs=self.trs.ctypes.data              
        p.p_ntrs_cut=self.ntrs_cut.ctypes.data              
        
        if not mDxD_deform is None:            
            p.p_mDxD_deform=np.array(self.mDxD_deform,dtype=np.float64).ctypes.data;
            
        _fuzzy_boundaries_condition_ids(byref(p));   
        return self;
        
        
            

    
    
    
            
         
