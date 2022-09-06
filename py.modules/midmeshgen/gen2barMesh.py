
import numpy as np
import scipy

from .mid_utils import mid_utils_ as mid_utils
from .minmorph import minmorph
from .hybrid_morph1D import hybrid_morph1D_

try:
    import matplotlib.pylab as plt
except ImportError:
    pass

def gen2barMesh(sensor, p):

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
    
    #%======== tool description ========
    core = midu.set_region_params(name = 'core', sigma = 0.000, mu = 5e4, rects = np.array([[0., 7.], [-senLen, senLen]]))
    #
    rcoil = midu.set_region_params('rcoil',4.6e2, 1, np.array([[11., 12.], [-senLen, senLen]]))
    gcoil = midu.set_region_params('gcoil',1.5e3, 1, np.array([[13.5, 14.91], [-senLen, senLen]]))
    hull = midu.set_region_params('hull',0.9e6, 1, np.array([[17., 21.], [-600., 600.]]))
    #
    #======== pipes description =======
    tub = midu.set_region_params('tubing',p.tbg.sigma,p.tbg.mu,np.array([[p.tbg.d/2.-p.tbg.th0, p.tbg.d/2.-p.tbg.th0+p.tbg.th], [-2000., 2000.]]))
    cas = midu.set_region_params('casing',p.csg.sigma,p.csg.mu,np.array([[p.csg.d/2.-p.csg.th0, p.csg.d/2.-p.csg.th0+p.csg.th], [-2000., 2000.]]))
    cas2 = midu.set_region_params('casing2',p.csg2.sigma,p.csg2.mu,np.array([[p.csg2.d/2.-p.csg2.th0, p.csg2.d/2.-p.csg2.th0+p.csg2.th], [-2000., 2000.]]))
    meshes = [core, rcoil, gcoil, hull, tub, cas, cas2]
    #%======== radial mesh ========
    rreg = []
    rreg.append( np.linspace(core["rects"][0][0], core["rects"][0][1], 12) )
    rreg.append( np.linspace(rcoil["rects"][0][0], rcoil["rects"][0][1], 2 ) )
    rreg.append( np.linspace(gcoil["rects"][0][0], gcoil["rects"][0][1], 2 ) )
    rreg.append( np.linspace(hull["rects"][0][0], hull["rects"][0][1], 12) )
    rreg.append( minmorph(np.array([p.tbg.d/2.-p.tbg.th0, p.tbg.d/2.-p.tbg.th0+p.tbg.th]), 0.25) )
    rreg.append( minmorph(np.array([p.csg.d/2.-p.csg.th0, p.csg.d/2.-p.csg.th0+p.csg.th]), 0.25) )
    rreg.append( minmorph(np.array([p.csg2.d/2.-p.csg2.th0, p.csg2.d/2.-p.csg2.th0+p.csg2.th]), 0.25) )
    rreg.append( [5e3] )
    rl = hybrid_morph1D_(rreg, 1.2).conj().T
    #%======== z mesh ========
    zreg = []
    zreg.append( np.dot(np.array([-2.5, -2.3]), 1e3) )
    zreg.append( minmorph(np.array([-senLen, senLen]), 4.) )
    zreg.append( np.dot(np.array([2.3, 2.5]), 1e3) )

    zl = hybrid_morph1D_(zreg, 1.3).conj().T



    core = midu.set_region_params(name = 'core', sigma = 0.000, mu = 5e4, rects = [[0., -senLen, 7., senLen]])
    rcoil = midu.set_region_params('srcoil',4.6e2,1, [[11., -senLen, 12., senLen]])
    gcoil = midu.set_region_params('sgcoil',1.5e3,1, [[13.5, -senLen, 14.91, senLen]])
    hull = midu.set_region_params('hull',0.9e6,1, [[17., -600, 21., 600.]] )
    tub = midu.set_region_params('tubing',p.tbg.sigma,p.tbg.mu, [[p.tbg.d/2.-p.tbg.th0, -2000, p.tbg.d/2.-p.tbg.th0+p.tbg.th, 2000.]])
    cas = midu.set_region_params('casing',p.csg.sigma,p.csg.mu, [[p.csg.d/2.-p.csg.th0, -2000, p.csg.d/2.-p.csg.th0+p.csg.th, 2000.]])
    cas2 = midu.set_region_params('casing2',p.csg2.sigma,p.csg2.mu, [[p.csg2.d/2.-p.csg2.th0, -2000, p.csg2.d/2.-p.csg2.th0+p.csg2.th, 2000.]])
    meshes = [core, rcoil, gcoil, hull, tub, cas, cas2]

    strdat = {}

    strdat["rl"] = rl
    strdat["zl"] = zl
    strdat["regions"] = meshes
    return strdat