
import numpy as np
import scipy

from .mid_utils import mid_utils_ as mid_utils
from .minmorph import minmorph
from .hybrid_morph1D import hybrid_morph1D_

try:
    import matplotlib.pylab as plt
except ImportError:
    pass

def barrier_in_out(tbg):
    tub_in=tbg.d/2.-tbg.th0
    tub_out=tub_in+tbg.th
    return (tub_in,tub_out);



def renorm_region(barriers):
    None


def gen2barMesh(sensor, p,frlzl=True):

    midu = mid_utils()
    _switch_val=sensor
    senLen = 0
    if False:
        pass
    elif _switch_val == 'SS':
        senLen = 60.
    elif _switch_val == 'LS':
        senLen = 160.
    else:
        print('Unexpected sensor type.')
        return []
    zmax=p.zmax
    zmin=-zmax
    #%======== tool description ========
    core = midu.set_region_params(name = 'core', sigma = 0.000, mu = 5e4, rects = np.array([[0., 7.], [-senLen, senLen]]))
    #
    rcoil = midu.set_region_params('rcoil',4.6e2, 1, np.array([[11., 12.], [-senLen, senLen]]))
    gcoil = midu.set_region_params('gcoil',1.5e3, 1, np.array([[13.5, 14.91], [-senLen, senLen]]))
    hull = midu.set_region_params('hull',0.9e6, 1, np.array([[17., 21.], [-600., 600.]]))
    #
    #======== pipes description =======
    (tub_in,tub_out)=barrier_in_out(p.tbg);
    (cas_in,cas_out)=barrier_in_out(p.csg);
    (cas2_in,cas2_out)=barrier_in_out(p.csg2);
    """"
    tub = midu.set_region_params('tubing',p.tbg.sigma,p.tbg.mu,np.array([[tub_in,tub_out], [zmin,zmax]]))
    cas = midu.set_region_params('casing',p.csg.sigma,p.csg.mu,np.array([[cas_in,cas_out], [zmin,zmax]]))
    cas2 = midu.set_region_params('casing2',p.csg2.sigma,p.csg2.mu,np.array([[cas2_in,cas2_out], [zmin,zmax]]))
    
    meshes = [core, rcoil, gcoil, hull, tub, cas, cas2]
    """
    if frlzl:
        #%======== radial mesh ========
        rreg = []
        rreg.append( np.linspace(core["rects"][0][0], core["rects"][0][1], 12) )
        rreg.append( np.linspace(rcoil["rects"][0][0], rcoil["rects"][0][1], 2 ) )
        rreg.append( np.linspace(gcoil["rects"][0][0], gcoil["rects"][0][1], 2 ) )
        rreg.append( np.linspace(hull["rects"][0][0], hull["rects"][0][1], 12) )
        rreg.append( minmorph(np.array([tub_in,tub_out]), 0.25) )
        rreg.append( minmorph(np.array([cas_in,cas_out]), 0.25) )
        rreg.append( minmorph(np.array([cas2_in,cas2_out]), 0.25) )
        rreg.append( [5e3] )
        rl = hybrid_morph1D_(rreg, 1.2).conj().T
        #%======== z mesh ========
        zreg = []
        zreg.append( np.dot(np.array([-2.5, -2.3]), 1e3) )
        zreg.append( minmorph(np.array([-senLen, senLen]), 4.) )
        zreg.append( np.dot(np.array([2.3, 2.5]), 1e3) )

        zl = hybrid_morph1D_(zreg, 1.3).conj().T
    else:
        (rl,zl)=(None,None)


    core = midu.set_region_params(name = 'core', sigma = 0.000, mu = 5e4, rects = [[0., -senLen, 7., senLen]])
    rcoil = midu.set_region_params('srcoil',4.6e2,1, [[11., -senLen, 12., senLen]])
    gcoil = midu.set_region_params('sgcoil',1.5e3,1, [[13.5, -senLen, 14.91, senLen]])
    hull = midu.set_region_params('hull',0.9e6,1, [[17., -600, 21., 600.]] )
    
    tub = midu.set_region_params('tubing',p.tbg.sigma,p.tbg.mu,[[tub_in,zmin,tub_out,zmax]])
    cas = midu.set_region_params('casing',p.csg.sigma,p.csg.mu,[[cas_in, zmin,cas_out,zmax]])
    cas2 = midu.set_region_params('casing2',p.csg2.sigma,p.csg2.mu,[[cas2_in,zmin,cas2_out,zmax]])

    hull2tubing = midu.set_region_params('hull2tubing',0,1, [[21.,zmin,tub_in,zmax]] )

    tub2cas = midu.set_region_params('tubing2casing',0,1,[[tub_out,zmin,cas_in,zmax]])
    cas2cas2 = midu.set_region_params('casing2casing',0,1,[[cas_out,zmin,cas2_in,zmax]])
    cas2inf = midu.set_region_params('casing2inf',0,1,[[cas2_out,zmin,5000,zmax]])
    #meshes = [core, rcoil, gcoil, hull, tub, cas, cas2]
    
    
    meshes = [core, rcoil, gcoil, hull,hull2tubing, tub,tub2cas, cas,cas2cas2, cas2,cas2inf]

    strdat = {}

    strdat["rl"] = rl
    strdat["zl"] = zl
    strdat["regions"] = meshes
    return strdat