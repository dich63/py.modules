import numpy as np
from MID.hallmesh.minmorph import minmorph
from MID.hallmesh.hybrid_morph1D import hybrid_morph1D_
from dataclasses import dataclass
from typing import List
from functools import reduce
from itertools import chain
from MID.hallmesh import well_construction as wc




@dataclass(frozen=True)
class Region:
    name: str
    sigma: float
    mu: float
    rects: List[List[float]]
    disabled: bool = False

    def __post_init__(self):
        if self.mu <= 0:
            raise Exception('Error: Mu must bu great zero')


@dataclass
class MeshData:
    rl: np.ndarray
    zl: np.ndarray
    regions: List[Region]


def gen_mesh(sensor, tubes: List, decentralization_mm: float = 0.) -> MeshData:
    senLen = 240.
    d = decentralization_mm

    tube = tubes[0]
    tube_inner_r, tube_outer_r = tube.d / 2. - tube.th0, tube.d / 2. - tube.th0 + tube.th

    # mesh description
    mmtube=senLen+20
    mmtube=300
    
    meshes = [
        Region(name='domain', sigma=0.,     mu=1.,      rects=[[-500.,   -1000.,      500.,     1000.      ]]),
        Region(name='tips',   sigma=0.,     mu=1.,    rects=[
            [-16. + d,    senLen - 20,  16. + d,      senLen + 20],
            [-16. + d,    -senLen - 20, 16. + d,     -senLen + 20],
        ]),
        Region(name='core',   sigma=0.,     mu=5e4,     rects=[[-10. + d,    -senLen,     10. + d,      senLen     ]]),
        Region(name='gcoil1', sigma=1.5e3,  mu=1.,      rects=[[13.736 + d,  -(senLen-20),     16. + d,      senLen - 20]]),
        Region(name='gcoil2', sigma=1.5e3,  mu=1.,      rects=[[-16. + d,    -(senLen-20),     -13.736 + d,  senLen - 20]]),
        Region(name='tube0',  sigma=tube.sigma, mu=tube.mu, rects=[
            [tube_inner_r,  -mmtube,  tube_outer_r, mmtube],
            [-tube_outer_r, -mmtube, -tube_inner_r, mmtube],
        ]),
    ]

    # ======== z mesh ========
    z_bounds = [
        [meshes[0].rects[0][1], meshes[0].rects[0][1] + 25.],
        minmorph(
            [
                [meshes[1].rects[0][1], meshes[1].rects[0][3]],
                [meshes[1].rects[1][1], meshes[1].rects[1][3]],
                [meshes[2].rects[0][1], meshes[2].rects[0][3]],
                [meshes[3].rects[0][1], meshes[3].rects[0][3]],
                [meshes[4].rects[0][1], meshes[4].rects[0][3]],
            ],
            2.5
        ),
        [meshes[0].rects[0][3], meshes[0].rects[0][3] + 25.],
    ]
    z_mesh = hybrid_morph1D_(z_bounds, 1.2)

    # ======== r mesh ========
    r_bounds = [
        [meshes[0].rects[0][0], meshes[0].rects[0][0] + 15.],
        minmorph(
            [
                meshes[5].rects[1][2],
                meshes[5].rects[1][0],
            ],
            .5
        ),
        minmorph(
            [
                meshes[3].rects[0][2],
                meshes[3].rects[0][0],
                meshes[1].rects[0][0],
                meshes[1].rects[0][2],
                meshes[2].rects[0][0],
                meshes[2].rects[0][2],
                meshes[4].rects[0][0],
                meshes[4].rects[0][2],
            ],
            .5
        ),
        minmorph(
            [
                meshes[5].rects[0][0],
                meshes[5].rects[0][2],
            ],
            .5
        ),
        [meshes[0].rects[0][2] - 15., meshes[0].rects[0][2]],
    ]
    r_mesh = hybrid_morph1D_(r_bounds, 1.2)
    return MeshData(r_mesh, z_mesh, meshes)


def _expand(grid):
    out = []
    for i in range(len(grid)-1):
        out.append(grid[i])
        out.append(0.5 * (grid[i]+grid[i+1]))
    out.append(grid[-1])
    return out

def hall_mesh_def(tube,decentralization_mm=0,Np=9,Ndeg=2):
    from MID.hall.hsolver import HCz_functor2
    from MID.hall.rectmesh import MeshData2mesh
    tubes = [
        wc.Tube(**tube),
    ]
    meshdata = gen_mesh(None, tubes, decentralization_mm)
    print('MeshData2mesh start...dc=',decentralization_mm)
    mesh,J=MeshData2mesh(meshdata);  
    print('HCz_functor2 start...')
    f=HCz_functor2(mesh,J,Np,Ndeg)
    print('HCz_functor2 end')
    return f

def hall_mesh_def0(tube,decentralization_mm=0,Np=9,Ndeg=2):
    from MID.hall.hsolver import HCz_functor2
    from MID.hall.rectmesh import MeshData2mesh
    tubes = [
        wc.Tube(**tube),
    ]
    meshdata = gen_mesh(None, tubes, decentralization_mm)
    
    mesh,J=MeshData2mesh(meshdata);    
    return mesh,J;



if __name__ == '__main__':
    from MID.hall.hsolver import *
    from MID.hall.rectmesh import *

    import matplotlib.pyplot as plt
    decentralizations = [0., 8.37, 11.84, 14.49, 16.73, 17.37]

    td={
        'd':88.9,
        'th':5.7,
        'mu':80
        };

    fu=hall_mesh_def(td,decentralizations[2]);







    tubes = [
        wc.Tube(d=88.9, th=5.7, sigma=3e6, mu=80.),
    ]

    tt=wc.Tube()


    
    meshdata = gen_mesh(None, tubes, decentralization_mm=decentralizations[0])
    m=meshdata

    mesh,J=MeshData2mesh(meshdata);

    hz=HCz_functor2(mesh,J,8,2);


    plt.figure()
    X, Y = np.meshgrid(meshdata.rl, meshdata.zl)
    ax = plt.subplot(1, 1, 1)
    ax.pcolor(X, Y, np.abs(Y*X), linewidths=.5,  edgecolors='k')
    plt.show()


pass
