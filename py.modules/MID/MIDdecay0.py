
import MID.mesh_regions as mr
from lipa.solvers import *
from klu_lipa import *
import jsonrpc.jsonclass as jsncls
import jsonrpc.sparse_marshal
from lipa.qpl import extend
from lipa.parallel import LIPA_solver_st,LIPA_solver_pp,lu_solver_factory
import asyn.SharedWorker as sw
import numpy as np
import sys
import copy


if sys.version_info.major==3:
    xrange=range
    unicode=str

inf=float("inf")
norm=np.linalg.norm;





def diff_sh(a,dt):
    ash1=np.arcsinh(a[1:]);
    ash0=np.arcsinh(a[:-1]);
    alpha=(ash1-ash0)/dt;
    return np.cosh(ash0)*alpha;


def derr2m(x,y):
    y=np.array(y,dtype='d');
    x=np.array(x,dtype='d');
    xn=np.abs(x)+np.abs(y)+1.0e-15;
    xd=x-y;
    return xd/xn;

def derr2m_norm(x,y,p=inf):   
    xn=norm(x,p)+norm(y,p)+1.0e-17;
    xd=norm(x-y,p);
    return xd/xn;


def cast_str(o):
    if type(o) in jsncls.string_types:
        o=eval(o)
    return o;

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
        pipes=get_rgs('tube');
    
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

def mesh_region_disabled_off(mesh):
    mesh=copy.copy(mesh);
    rgsn,rgs=[],mesh['regions'];
    for r in rgs:        
        if r and (not r.get('disabled',False)):
            rgsn.append(r);
    mesh['regions']=rgsn;
    return mesh;


class MIDdecay0(object):
    def __init__(self,mesh,gcoil='sgcoil',rcoil='srcoil',rcoil_inflate=1e-3):
        print('mesh preparing...')
        mesh=reparse_mesh(mesh);
        mesh=mesh_region_disabled_off(mesh);
        self.mesh=mr.MID_mesh(mesh);

        self.gcoil_name=gcoil;
        self.rcoil_name=rcoil;

        self.gcoil=self.mesh.region_mask(gcoil)
        self.rcoil=self.mesh.region_mask(rcoil,rcoil_inflate)
        self.fmd=True
        self.perf={}
        print('end')

    def __getitem__(self,name):
        return self.mesh.regions[name];



    def init(self,**kws):
        print('init...')
        tic=sw.Tic();
        self.solver=None
        opts=extend(kws,{'nd':4,'tic_corr':0,'fuzzy':None})
        LIPA_solver=cast_str(opts.pop('LIPA_solver',LIPA_solver_pp));
        lu_factory=cast_str(opts.pop('lu_factory',lu_solver_factory()));        
        print(LIPA_solver)    

        if self.fmd:
            self.mesh.make_vxs_deformation();

        scheme=opts.pop('scheme',1)
        fuzzy=opts.pop('fuzzy',None)
        printf('FEM -scheme: 0x%X\n',scheme);
        (A,C,V)=self.mesh.makeACV(scheme=scheme,fuzzy=fuzzy)
        self.dt=opts.get('dt',1.0);
        self.C=C;
        self.A=A;
        self.V=V;
        #(A,C)=(A.tocsc(),C.tocsc())
        #A.data=-A.data;
        tic.start()
        solver=LIPA_solver((A,C),SolverFactory=lu_factory,**opts);        
        solver.reset_corr(opts['tic_corr']);
        self.solver=solver;
        t=tic.sec();
        self.perf['tloading']=tlu=solver.tloading;
        printf('\r loading... elapsed: %f sec\n',tlu);        
        self.perf['tLU']=tlu=solver.tinit;        
        printf('\r LU elapsed: %f sec\n',tlu);

    def reset_corr(self,n):
        self.solver.reset_corr(n);

    def decay(self,count=100,J=1.0,fnparray=False,fctx=0):
        solver=self.solver;
        C=self.C.tocsc();

        gcoil=C*self.gcoil.astype('d');
        rcoil=C*self.rcoil.astype('d');

        solver.reset_J(J*gcoil);
        
        jrs=[0.0]
        count=int(count)
        if fctx:
            (i,jsr)=solver.form_tic_rep(rcoil,rep=count+1,nd=1);
            i();
            return jsr;
        else:       
            for n in range(count):            
                solver.step();
                jr=solver.form(rcoil,1);            
                jrs.append(jr);
                #printf('\r %d of %d normxx=%f',n,count,norm(solver.x));
                printf('\r %d of %d ',n,count);
            printf('\n end \n')   
            return np.array(jrs) if fnparray else jrs;


    def getAC(self,raw=False):                
        print('raw:',raw);
        if raw:
            return (self.A.tocoo(),self.C.tocoo())
        else:
            (a,c)=(self.A.tocsc(),self.C.tocsc());
            a.eliminate_zeros();
            c.eliminate_zeros();
            return (a,c);
                    

    def decay_test(self,count=100,count_front=-1,J=1e10,fnparray=True,**others):
        solver=self.solver;
        C=self.C.tocsc()
        A=self.A.tocsc();
        C.eliminate_zeros();
        A.eliminate_zeros();
        V=self.V.tocsc();
        dt=self.dt;
        gcoil=self.gcoil.astype('d');
        rcoil=self.rcoil.astype('d');
        field_step=others.get('field_step',-1);
        field_step=int(field_step);


        jgcoil=J*(C*gcoil);
        jrcoil=C*rcoil;

        jgcoil=self.mesh.form_volume(self.gcoil_name,eps=1e-4,scheme=1,with_deform=0);
        jrcoil=self.mesh.form_volume(self.rcoil_name,eps=1e-4,scheme=1,with_deform=0);

        #jgcoil=J*(V*gcoil);
        #jrcoil=V*rcoil;

        solver.reset_J(jgcoil);

        jfactor=1.0;

        def jsol():
            jf=-(A*solver.xn[0]) +(jfactor)*jgcoil;
            return (np.dot(rcoil,jf),jf);


        
        
        jrs,jrst,jrAs,df1,t,xx,ssf,ssf2=[],[],[],[],[],[],[],[]
        
        count=int(count)
        count_front=int(count_front);
        if count_front<0:
            count_front=count+1

        tt=0;
        tic=sw.Tic();

        for n in range(count):
            if n>count_front:
                count_front,jfactor=count+1,0.0;
                solver.reset_J(0);

            tic.start();
            solver.step();
            tt+=tic.sec();
            f1=solver.vxn[1];            
            jr=solver.form(jrcoil,1);       
            jr2=solver.form(jrcoil,2);
            jr3=solver.form(jrcoil,3);
            jrs.append(jr);

            jrA=solver.form(jrcoil,0);
            jrAs.append(jrA);

            #jrt=-jsol();
            #(tmp,f1e)=jsol();
            #jrt=-tmp;        
            (jrt,f1e)=jsol();    
            jrst.append(jrt);
            df1.append(100*derr2m_norm(C*f1,f1e,inf))
            tc=n*dt;
            
            ss,ss2=0.0,0.0;
            #ss=-1.e-7
            if n:
                ss=tc*jr2/jr;
                #ss=-tc*abs(jr2);                
                ss2=tc*jr3/jr2+0;
                
            ssf.append(-ss); 
            ssf2.append(-ss2);                           
            t.append(tc)            
            if (field_step>0) and (n%field_step==0) :
                xn=[ v for v in self.solver.xn ]
                jx=C*xn[1];
                xn.append(jx);
                xn=[ v.astype('float32') for v in xn ];
                xx.append(xn);
            printf('\r %d of %d t=%f <tic-time>=%f',n,count,tc,tt/(n+1));

        printf('\n end \n')   

        self.perf['tTic']=tt/count;

        jrs=jrs[0:-1];
        jrst=jrst[0:-1];
        t=t[1:];
        jra=[ (jrAs[k+1]-jrAs[k])/dt for k in range(count-1) ]
        jrash=diff_sh(jrAs,dt);

        errd=derr2m(jrs,jra)*100;
        erre=derr2m(jrs,jrst)*100;

        errdsh=derr2m(jrs,jrash)*100;

        if not fnparray:
            errd=errd.tolist();
            errdsh=errdsh.tolist();
            erre=erre.tolist();

        errdm=norm(errd,inf);
        errdshm=norm(errdsh,inf);
        errem=norm(erre,inf);
        df1m=norm(df1,inf);

        r=np.array(jrs) if fnparray else jrs;
        rt=np.array(jrst) if fnparray else jrst;
        ra=np.array(jra) if fnparray else jra;
        ra=np.array(jra) if fnparray else jra;
        a=np.array(jrAs) if fnparray else jrAs;

        res={'t':t,'jr':r,'jre':rt,'jrd':ra,'a':a,'errdp':errd,'errdhp':errdsh,'errep':erre,'maxerrdp':errdm,'maxerrdhp':errdshm,'maxerrep':errem,'errp_jpe':df1,'max_errp_jpe':df1m,'xx':xx,'ssf':ssf,'ssf2':ssf2};
        res['perf']=self.perf;
        return res;

    def fields(self,count=1,J=1.0,fnparray=False,**others):
        solver=self.solver;
        solver.reset_J(J*self.gcoil);
        jrs=[0.0]
        count=int(count)
        for n in range(count):
            solver.step();
        xn=self.solver.xn;                        
        return xn if fnparray else   [ v for v in xn ];


from jsobj import *

class MID3barrier_decay(MIDdecay0):
    def __init__(self,mesh,gcoil='gcoil',rcoil='rcoil',rcoil_inflate=1e-3):
        mesh=reparse_mesh(mesh)
        mesh=mesh_region_disabled_off(mesh);
        print('nominal_barriers mesh preparing...')
        update_nominal_barriers(mesh['regions']);
        super(MID3barrier_decay,self).__init__(mesh,gcoil,rcoil,rcoil_inflate);
        

    def mesh_reconstruct(self,**kws):
        print('mesh_reconstruct...')
        barriers=kws.get('barriers',{});        
        print(kws)        
        update_nominal_barriers(self.mesh.regions.regions,barriers)
        self.mesh.make_vxs_deformation();
        self.fmd=False;
        regions_new=self.mesh.regions.regions;
        vxs_new=self.mesh.vxs
        vxs_new=vxs_new.astype('float32');
        mesh_new=copy.copy(self.mesh.mesh);
        mesh_new['vxs']=vxs_new;
        mesh_new['regions']=regions_new;
        self.new_mesh=mesh_new;
        return mesh_new;
    
        
    def init(self,**kws):
        #barriers=kws.get('barriers',{});
        #update_nominal_barriers(self.mesh.regions.regions,barriers)
        
        if kws.get('barriers',{}):
            self.mesh_reconstruct(**kws);
        super(MID3barrier_decay,self).init(**kws);
        return self;

    def reset_field(self,x=0):
        self.solver.x=x;
    
    def all_currents(self,count=100,Jg=1.0,
                     jj=jsc(core=[],rcoil=[],pipe1=[],pipe2=[],pipe3=[]),
                     cutoff_rect=[-np.inf,+np.inf,-2,+2],
                     rcoil_inflate=1e-8):
    
        def isiter(j):
            try:
                iter(j)
                return True;
            except:
                return False;
        
        #jlist=('rcoil','pipe1','pipe2','pipe3')
        
        jlist=[n for n in jj];
        
        solver=self.solver;
        C=self.C.tocsc();
        A=self.A.tocsc();
        V=self.V.tocsc();
        gname=self.gcoil_name
        #gcoil=self.mesh.form_volume(gname,eps=1e-4,scheme=1,with_deform=0);
        gcoil=C@self.mesh.region_mask(gname,rcoil_inflate)
        
        cq=Jg if isiter(Jg) else (Jg,)
        
        Jq=({'qp':({'c':cq},),'x':gcoil.flatten()},);   
        jgl=[] 
        
        jm=jso();
        for n in jlist:
            jj[n]=jj._def(n,[]);
            #jf[n]=self.mesh.form_volume(n,eps=1e-4,scheme=1,with_deform=0);
            mask=self.mesh.region_mask(n,rcoil_inflate,cutoff_rect=cutoff_rect);
            jm[n]=mask!=0
            
            
        jj[gname]=jj._def(gname,[]);
        
        
        jlist=[n for n in jj];
        
        solver.reset_J(Jq);
        
        for c in range(count):
            xx=solver.step();
            jxx=C@xx[1];
            for n in jlist:
                if n==gname:
                    jj[n]+=[solver.solvers[0].qpj.jl[0].j0()]
                else:                      
                    mask=jm[n];
                    smm=np.sum(mask);
                    v=np.sum(jxx[mask]);             
                    jj[n]+=[v]
                
            
            printf('\r %d of %d ',c,count);
        
        return jj;    



