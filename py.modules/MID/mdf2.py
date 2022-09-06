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

def printf(*args):
    f=args[0];
    a=args[1:];
    s=f% tuple(a);
    sys.stdout.write(s);

def pipes_prepare_regions(regions):
        
    rg=mr.MID_regions(regions);
    rgs=rg.regions;

    def get_rgs(name):
        p=[]
        for r in rgs:
            n=r['name']
            if n.find(name)==0:
                p.append(r);
        return p
        
    pipes=get_rgs('pipe');        
    
    if not pipes:
        return ([],[],rgs);

    intfs=get_rgs('intf');      

    if  intfs:
        return (pipes,intfs,rgs);

    rb=rg.rgnmap['hull']['rects'][0];

    

    i=0;
    
    for pn in pipes:
        re=rg.rgnmap[pn['name']]['rects'][0];
        left,right=rb[2],re[0];
        zmin,zmax=re[1],re[3];
        region={'name':'intf'+str(i),'rects':[[left,zmin,right,zmax]]}
        intfs.append(region);
        i+=1
        rb=re;
    inf=5000;#np.inf
    left,right=rb[2],inf;
    zmin,zmax=rb[1],rb[3];
    region={'name':'intf'+str(i),'rects':[[left,zmin,right,zmax]]}
    intfs.append(region);
    
    
    rgs+=intfs;
    return (pipes,intfs,rgs);

def reset_nominal_barriers_pipes(regions,barriers={}):

    (pipes,intfs,rgs)=pipes_prepare_regions(regions);
    lenp=len(pipes);

    if not (lenp+1==len(intfs)) :
        return None;
        
    regions=mr.MID_regions(rgs);
    
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
        for i in range(lenp):
            pn=pipes[i]['name'];
            ln=intfs[i]['name'];
            rn=intfs[i+1]['name'];
            reset_triplet(pn,ln,rn,barriers.get(pn,{}));
    return regions.regions;

"""
        reset_triplet('tubing','hull-tubing','tubing-casing',barriers.get('tbg',{}))
        reset_triplet('casing','tubing-casing','casing-casing2',barriers.get('csg',{}))
        reset_triplet('casing2','casing-casing2','casing2-infinity',barriers.get('csg2',{}))         
"""
   


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
    regions[:]=reset_nominal_barriers3(regions,barriers);

def update_nominal_barriers(regions,barriers={}):
    rs=reset_nominal_barriers_pipes(regions,barriers);
    if not rs:
        rs=reset_nominal_barriers3(regions,barriers);
    regions[:]=rs;

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
    def __init__(self,mesh,gcoil='sgcoil',rcoil='srcoil',rcoil_inflate=1e-3):
        print('mesh preparing...')
        mesh=reparse_mesh(mesh)
        self.mesh=mr.MID_mesh(mesh);
        self.gcoil=self.mesh.region_mask(gcoil)
        self.rcoil=self.mesh.region_mask(rcoil,rcoil_inflate)
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
        (A,C,V)=self.mesh.makeACV(scheme=opts.pop('scheme',1))
        self.dt=opts.get('dt',1.0);
        self.C=C;
        self.A=A;
        self.V=V;
        #A.data=-A.data;
        self.solver=LIPA_solver((A,C),SolverFactory=lu_factory,**opts);
        print('end')

    def reset_corr(self,n):
        self.solver.reset_corr(n);

    def decay(self,count=100,J=1.0,fnparray=False):
        solver=self.solver;
        C=self.C.tocsc();

        gcoil=C*self.gcoil.astype('d');
        rcoil=C*self.rcoil.astype('d');

        solver.reset_J(J*gcoil);
        
        jrs=[0.0]
        count=int(count)
        for n in range(count):
            solver.step();
            jr=solver.form(rcoil,1);
            jrs.append(jr);
            printf('\r %d of %d',n,count);
        printf('\n end \n')   

        return np.array(jrs) if fnparray else jrs;

    def decay_test(self,count=100,J=1e10,fnparray=True):
        solver=self.solver;
        C=self.C.tocsc();
        A=self.A.tocsc();
        V=self.V.tocsc();
        dt=self.dt;
        gcoil=self.gcoil.astype('d');
        rcoil=self.rcoil.astype('d');

        jgcoil=J*(C*gcoil);
        jrcoil=C*rcoil;

        solver.reset_J(jgcoil);

        def jsol():
            jf=(A*solver.vxn[0]) +(1)*jgcoil;
            return np.dot(rcoil,jf);


        
        
        jrs,jrst,jrAs=[],[],[]
        
        count=int(count)

        for n in range(count):
            solver.step();
            jr=solver.form(jrcoil,1);       

            jrs.append(jr);

            jrA=solver.form(jrcoil,0);
            jrAs.append(jrA);

            jrt=-jsol();
            jrst.append(jrt);
            printf('\r %d of %d',n,count);
        printf('\n end \n')   
        jrs=jrs[0:-1];
        jrst=jrst[0:-1];

        jra=[ (jrAs[k+1]-jrAs[k])/dt for k in range(count-1) ]

        r=np.array(jrs) if fnparray else jrs;
        rt=np.array(jrst) if fnparray else jrst;
        ra=np.array(jra) if fnparray else jra;
        a=np.array(jrAs) if fnparray else jrAs;

        return (r,rt,ra,a);

    def fields(self,count=1,J=1.0,fnparray=False):
        solver=self.solver;
        solver.reset_J(J*self.gcoil);
        jrs=[0.0]
        count=int(count)
        for n in range(count):
            solver.step();

        xn=self.solver.xn;                        
        return xn if fnparray else   [ v for v in xn ];


class MID3barrier_decay(MIDdecay0):
    def __init__(self,mesh,gcoil='gcoil',rcoil='rcoil',rcoil_inflate=1e-3):
        mesh=reparse_mesh(mesh)
        print('nominal_barriers3 mesh preparing...')
        update_nominal_barriers(mesh['regions']);
        super(MID3barrier_decay,self).__init__(mesh,gcoil,rcoil,rcoil_inflate);
        
    def init(self,**kws):
        barriers=kws.get('barriers',{});
        update_nominal_barriers(self.mesh.regions.regions,barriers)
        super(MID3barrier_decay,self).init(**kws);
    def reset_field(self,x=0):
        self.solver.x=x;


