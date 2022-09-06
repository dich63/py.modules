# -*- coding: utf-8 -*-
"""
Created on Sat Mar 25 21:11:54 2017

@author: dich
"""

import jsonrpc.jsonclass as jsncls
import jsonrpc.sparse_marshal
import copy
from MID.rect_triangulation import *


def mesh_refine(mesh,expand_z=0,expand_r=0,centroid=0,unit_scale=1000,fcopy=False):
    
    if fcopy:
        mesh=copy.deepcopy(mesh);
        
    rl=mesh["init_opts"]["rl"];
    zl=mesh["init_opts"]["zl"];
    lbound=mesh.get("lbound",0);
    for k in range(expand_z):
        zl=expand1(zl);
    for k in range(expand_r):
        rl=expand1(rl);        
    (trs,vxs,fbc)=make_triangulation_c(rl,zl,lbound,unit_scale=unit_scale);
    
    (mesh['trs'],mesh['vxs'],mesh['fbc'],mesh['lbound'])=(trs,vxs,fbc,lbound);
    
    return mesh
    
#def make_triangulation_c(rl,zl,lbound=0,fbc_mask=(1,1,1,1),unit_scale=1.0,info={}):
    
