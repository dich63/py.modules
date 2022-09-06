import MID.mesh_regions as mr
from lipa.solvers import *
import jsonrpc.jsonclass as jsncls
import jsonrpc.sparse_marshal
from lipa.qpl import extend
from lipa.parallel import LIPA_solver_st,LIPA_solver_pp,lu_solver_factory
import numpy as np
import sys
import copy

if sys.version_info.major==3:
    xrange=range
    unicode=str



def barriers3_prepare_regions(regions):
    rg=mr.MID_regions(regions)

    if 'hull-tubing' in rg.rgnmap :
        return rg.regions;

    zmax=2500;
    zmin=-zmax;
    hull=rg.rgnmap['hull']['rects'][0];
    tub=rg['tubing']['rects'][0];
    cas=rg['casing']['rects'][0];
    cas2=rg['casing2']['rects'][0];
    
    hr=hull[2];
    tl=tub[0];
    tr=tub[2];
    cl=cas[0];
    cr=cas[2];
    c2l=cas2[0];
    c2r=cas2[2];

    ri=[{'name':'hull-tubing','rects':[[hr,zmin,tl,zmax]]},
        {'name':'tubing-casing','rects':[[tr,zmin,cl,zmax]]},
        {'name':'casing-casing2','rects':[[cr,zmin,c2l,zmax]]},
        {'name':'casing2-infinity','rects':[[c2r,zmin,5000,zmax]]} ]

    r=rg.regions;
    r+=ri
    return r;

def reset_nominal_barriers3(regions,barriers={}):
    regions=mr.MID_regions(barriers3_prepare_regions(regions))
    
    def update_params(r,b):
        r['sigma']=b.get('sigma',r.get('sigma',0));
        r['mu']=b.get('mu',r.get('mu',1));
        r['disabled']=b.get('disabled',r.get('disabled',False));
    def barrier_in_out(tbg):
        th=tbg["th"];
        th0=tbg.get("th0",th);
        d=tbg["d"];
        tub_in=d/2.-th0
        tub_out=tub_in+th
        return (tub_in,tub_out);

    def reset_triplet(cn,ln,rn,b):
        if b:
            rg=regions[cn];
            update_params(rg,b);        
            if not b.get('nominal',False):
                rc=rg['rects'][0];
                rl=regions[ln]['rects'][0];
                rr=regions[rn]['rects'][0];
                [l,r]=barrier_in_out(b);
                
                rc[0],rc[2]=l,r;
                rl[2],rr[0]=l,r;               

    if barriers:
        reset_triplet('tubing','hull-tubing','tubing-casing',barriers.get('tbg',{}))
        reset_triplet('casing','tubing-casing','casing-casing2',barriers.get('csg',{}))
        reset_triplet('casing2','casing-casing2','casing2-infinity',barriers.get('csg2',{}))         

    return regions.regions;

def update_nominal_barriers3(regions,barriers={}):
    regions[:]=reset_nominal_barriers3(regions,barriers)

def reparse_mesh(mesh):
    #print('reparse_mesh:...',type(mesh))
    l=6;
    if (type(mesh)==str) or (type(mesh)==unicode) :
            mesh=mesh.strip()
            if mesh.find(':file:',0,l)==0:
                if mesh[l:l+1]==':':
                    l=7;
                mesh=open(mesh[l:],'r').read()
            mesh=jsncls.decode(mesh)
    #print('reparse_mesh:Ok')
    return mesh


class MIDdecay0(object):
    def __init__(self,mesh,gcoil='sgcoil',rcoil='srcoil'):
        print('mesh preparing...')
        mesh=reparse_mesh(mesh)
        self.mesh=mr.MID_mesh(mesh);
        self.gcoil=self.mesh.region_mask(gcoil)
        self.rcoil=self.mesh.region_mask(rcoil)
        print('end')

    def __getitem__(self,name):
        return self.mesh.regions[name];

    def init(self,**kws):
        print('init...')
        self.solver=None
        opts=extend(kws,{'nd':2})
        LIPA_solver=opts.pop('LIPA_solver',LIPA_solver_pp);
        lu_factory=opts.pop('lu_factory',lu_solver_factory())
        self.mesh.make_vxs_deformation();
        (A,C)=self.mesh.makeAC(scheme=opts.pop('scheme',1))
        self.dt=opts.get('dt',1.0);
        #A.data=-A.data;
        self.solver=LIPA_solver((A,C),SolverFactory=lu_factory,**opts);
        print('end')

    def reset_corr(self,n):
        self.solver.reset_corr(n);

    def decay(self,count=100,J=1.0,fnparray=False):
        solver=self.solver;
        solver.reset_J(J*self.gcoil);
        rcoil=self.rcoil;
        jrs=[0.0]
        count=int(count)
        for n in range(count):
            solver.step();
            jr=solver.form(rcoil,1);
            jrs.append(jr);
        return np.array(jrs) if fnparray else jrs;

    def fields(self,count=5,J=1.0,fnparray=False):
        solver=self.solver;
        solver.reset_J(J*self.gcoil);
        jrs=[0.0]
        count=int(count)
        for n in range(count):
            solver.step();
        xn=self.solver.xn;                        
        return xn if fnparray else   [ v for v in xn ];


class MID3barrier_decay(MIDdecay0):
    def __init__(self,mesh,gcoil='gcoil',rcoil='rcoil'):
        mesh=reparse_mesh(mesh)
        print('nominal_barriers3 mesh preparing...')
        update_nominal_barriers3(mesh['regions']);
        super(MID3barrier_decay,self).__init__(mesh,gcoil,rcoil);
        
    def init(self,**kws):
        barriers=kws.get('barriers',{});
        update_nominal_barriers3(self.mesh.regions.regions,barriers)
        super(MID3barrier_decay,self).init(**kws);
    def reset_field(self,x=0):
        self.solver.x=x;


