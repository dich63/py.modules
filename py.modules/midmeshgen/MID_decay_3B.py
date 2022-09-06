# -*- coding: utf-8 -*-
"""
Created on Fri Jan  8 04:18:53 2016

@author: dich6_000
"""


import jsonrpc.jsonclass as jc
import jsonrpc.sparse_marshal
from jsonrpc.jsonclass import jslike
import MID.mesh_regions as mr
import asyn.SharedWorker as sw
from scipy.spatial import Delaunay

from scipy.sparse import  linalg as sla
import MID.MIDdecay0 as mmd
from lipa.parallel import LIPA_solver_st,LIPA_solver_pp
#import affinity
import sys
import numpy as np
#from midmeshgen.gen2barMesh2 import gen2barMesh
from midmeshgen.gen2barMesh import  gen2barMesh as gen2barMesh_bag
from midmeshgen.gen2barMeshDich import  gen2barMesh as gen2barMesh
from  MID.rect_triangulation import make_triangulation_c

near_t = np.array([0.000133, 0.000178, 0.000237, 0.000316, 0.000365, 0.000422, 0.000488, 0.00056200000000000011,
                   0.00064900000000000005, 0.00075, 0.000866, 0.001, 0.00115, 0.00133, 0.0015400000000000001,
                   0.0017800000000000001, 0.00206, 0.00237, 0.0027400000000000002, 0.00316, 0.00365, 0.00422,
                   0.00487, 0.00562, 0.00649, 0.0075, 0.008660000000000001, 0.01, 0.0115, 0.013300000000000001,
                   0.0154, 0.0178, 0.0206, 0.0237, 0.0274, 0.0316, 0.0365, 0.0422, 0.048700000000000007,
                   0.056200000000000007, 0.064900000000000013, 0.075])

if sys.version_info.major==3:
    xrange=range

def delta(x,y):
    return np.linalg.norm(np.array(x)-np.array(y))

def derr2m(x,y):
    y=np.array(y,dtype='d');
    x=np.array(x,dtype='d');
    xn=np.abs(x)+np.abs(y)+1.0e-15;
    xd=x-y;
    return xd/xn;

def make_barriers_3B(**kw):
    barriers={
            "tbg":{
                   "d":140.,
                   "th":7.,
                   "th0":7.,
                   "sigma":6.7e6,
                   "mu":35.9
                },
            "csg":{
                   "d":240.,
                   "th":10.,
                   "th0":10.,
                   "sigma":6.7e6,
                   "mu":30.6
                },
            "csg2":{
                   "d":340.,
                   "th":10.,
                   "th0":10.,
                   "sigma":6.7e6,
                   "mu":30.6
                },
            "sensor":'LS',
            'zmax':2500
            }
    for k in barriers:
        v=kw.get(k,None);
        if not v is None:
            if (type(barriers[k]) is dict) and (type(v) is dict):
                barriers[k].update(v)
            else:
                barriers[k]=v;
    for k in kw:
        v=barriers.get(k,None);
        if  v is None:
            barriers[k]=kw[k]
        elif (type(kw[k]) is dict) and (type(v) is dict):
            v.update(kw[k])
    #barriers.update(kw)
    return barriers;

def make_bc(nxny,fbc_mask):
    fbc=np.zeros(nxny,dtype=np.uint16)
    fbc[ :, 0]=fbc_mask[0]
    fbc[ :,-1]=fbc_mask[1]
    fbc[ 0, :]=fbc_mask[2]
    fbc[-1, :]=fbc_mask[3]
    return fbc




def make_triangulation(rl,zl,fbc_mask,unit_scale):    
    
    [X, Y] = np.meshgrid(rl,zl)
    nxny = np.shape(X)
    fbc=make_bc(nxny,fbc_mask)
    t = np.zeros( (np.size(X), 2) )

    t[:,0] = X.reshape( np.size(X) )
    t[:,1] = Y.reshape( np.size(Y) )

    vxs = np.ascontiguousarray(t, dtype="d")
    vxs /= unit_scale

    tri = Delaunay(t)
    trs = np.ascontiguousarray(tri.simplices, dtype=np.uint32) #trs

    '''
    norm=np.linalg.norm
    ff = [];
    
    nmi0=np.min(vxs[:,0])
    nmi1=np.min(vxs[:,1])
    nma0=np.max(vxs[:,0])
    nma1=np.max(vxs[:,1])

    
    eps=1e-14;

    #xu=np.unique(vxs[:,0])
    #yu=np.unique(vxs[:,1])
    
    dist = lambda x: np.abs(x)<eps

    for one in vxs:
        #if (nmi0 in one) or (nmi1 in one) or (nma0 in one) or (nma1 in one):
        if dist(nmi0- one[0]) or dist(nmi1- one[1]) or dist(nma0- one[0]) or dist(nma1- one[1]):
            ff.append(1)
        else:
            ff.append(0)

    fbc0 = np.array(ff)
    nn=norm(np.ravel(fbc)-np.ravel(fbc0));
    ns=np.sum(np.abs((np.ravel(fbc))));
    ns0=np.sum(np.abs((np.ravel(fbc0))));
    nsd=np.sum(np.abs((np.ravel(fbc)-np.ravel(fbc0))));
    fbc=fbc0;
    '''


    return (trs,vxs,fbc)



def make_MID_3B_mesh(barriers=make_barriers_3B(),**opts):
    
    params=jslike(barriers);
    params.t = near_t
    params.supp = np.nonzero(np.logical_and(near_t > 0.5e-3, near_t<=21e-3))
    params.tsupp = near_t[params.supp[0]-1]
    #params.dt = 0.1e-3
    #params.sensor = 'LS'
    params.wCollar = 140.
    params.level = 0.01
    grd = {}
    unit_scale=1000;

    

    ftriangle=opts.get('triangulate',True);

    sy0 = gen2barMesh((params.sensor), params,frlzl=ftriangle)
    
    (rl,zl)=(sy0['rl'],sy0['zl'])

    
    if ftriangle:
        make_tri=make_triangulation_c if opts.get('centroid',True) else make_triangulation;
        fbc_mask=opts.get('fbc_mask',(1,1,1,1))
        (trs,vxs,fbc)=make_tri(rl,zl,fbc_mask=fbc_mask,unit_scale=unit_scale);
        grd["trs"] = trs
        grd["vxs"] = vxs
        grd["fbc"] = fbc
    
    grd["regions"] = sy0["regions"]
    grd["scheme"] = 1
    
    grd["units"] = {"regions": 0.001}
    grd["lbound"] = 0
    grd["pade"] = {"n": 2, "m": 4}
    grd["rlzl"] = [rl, zl]

    smbg=opts.get('sigma_mu_bg',None)
    if not smbg is None:
        grd['sigma_mu_bg'] = smbg;
    return grd



def get_MID_decay_3B(barriers=make_barriers_3B(),dt=0.001,count=100,gcoil='sgcoil',rcoil='srcoil',**opts):
    
    params=jslike(barriers);
    params.t = near_t
    params.supp = np.nonzero(np.logical_and(near_t > 0.5e-3, near_t<=21e-3))
    params.tsupp = near_t[params.supp[0]-1]
    #params.dt = 0.1e-3
    #params.sensor = 'LS'
    params.wCollar = 140.
    params.level = 0.01

    syba = gen2barMesh_bag((params.sensor), params)
    sydi=  gen2barMesh((params.sensor), params)

    errz=delta(syba['zl'],sydi['zl'])
    errr=delta(syba['rl'],sydi['rl'])

    if opts.get('bag',False):
        sy0 = gen2barMesh_bag((params.sensor), params)
    else:
        sy0 = gen2barMesh((params.sensor), params)
    
    (rl,zl)=(sy0['rl'],sy0['zl'])

    fbc_mask=opts.get('fbc_mask',(1,1,1,1))

    (trs,vxs,fbc)=make_triangulation(rl,zl,fbc_mask);
    
    grd = {}
    grd["fbc"] = fbc
    grd["pade"] = {"n": 2, "m": 4}
    grd["regions"] = sy0["regions"]
    grd["scheme"] = 1
    grd["trs"] = trs
    grd["vxs"] = vxs
    grd["units"] = {"regions": 0.001}
    grd["lbound"] = 0

    smbg=opts.get('sigma_mu_bg',None)
    if not smbg is None:
        grd['sigma_mu_bg'] = smbg;

    

    LIPA_solver_thread_model = LIPA_solver_pp if opts.get('mt',True) else LIPA_solver_st;
    tic=sw.Tic()
    tic.start()
    ms=mmd.MIDdecay0(grd)
    t=tic.sec()
    print(np.round(t),'sec MIDdecay0')
    tic.start()
    ms.init(dt=dt,LIPA_solver=LIPA_solver_thread_model,**opts)
    t=tic.sec()
    print(np.round(t),'sec :ms.init(',LIPA_solver_thread_model,')')
    tic.start()
    jst1=ms.decay(count)
    t=tic.sec()
    print(np.round(t),'sec :count=',count)
    if opts.get('binout',True):
        jst1=np.array(jst1);
    return jst1




def json_get_MID_decay_3B(s):
    params=jc.decode(s);
    def update(n,dflt):
        params[n]=params.get(n,dflt)
    update('barriers',make_barriers_3B());
    update('dt',0.001);
    update('count',100);
    update('gcoil','sgcoil');
    update('rcoil','srcoil');
    j=get_MID_decay_3B(**params);
    s=jc.encode(j);
    return s;

def json_make_barriers_3B(s):
    params=jc.decode(s);
    params=make_barriers_3B(**params);
    s=jc.encode(params);
    return s;

def json_make_MID_3B_mesh(s):
    params=jc.decode(s);
    params=make_MID_3B_mesh(**params);
    s=jc.encode(params);
    return s;


    #barriers=make_barriers_3B(),dt=0.001,count=100,gcoil='sgcoil',rcoil='srcoil',**opts):

