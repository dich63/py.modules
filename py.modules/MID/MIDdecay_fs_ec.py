
import MID.mesh_regions as mr
import MID.e_circuit as ec
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
import scipy.constants as scs

if sys.version_info.major==3:
    xrange=range
    unicode=str

inf=float("inf")
norm=np.linalg.norm;


def region_by_name(mesh,name):
    rgs= mesh['regions'];
    for r in rgs:
        if r['name']==name:
            return r;
    return None;



def mulxx(A,xx,dtype = np.complex):
    yy=np.empty_like(xx,dtype = dtype)
    (N,t)=xx.shape;
    for n in range(N):
        yy[n,:]=A*xx[n,:]

    return yy;


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
        region={'name':'intf'+str(i),'rects':[[left,zmin,right,zmax]],'disabled':77}
        intfs.append(region);
        i+=1
        rb=re;
    inf=5000;#np.inf
    left,right=rb[2],inf;
    zmin,zmax=rb[1],rb[3];
    region={'name':'intf'+str(i),'rects':[[left,zmin,right,zmax]],'disabled':77}
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

        def update_param(name,dft):
            r[name]=b.get(name,r.get(name,dft));

        update_param('sigma',0);
        update_param('mu',1.0);
        update_param('z2r',1.0);
        update_param('dc',0);
        update_param('dcr',0);
        update_param('disabled',False);

    def barrier_in_out(tbg,rin,rext):

        d=tbg.get("d",2.0*rext);
        th=tbg.get("th",rext-rin);
        th0=tbg.get("th0",th);
        
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
                [l,r]=barrier_in_out(b,rc[0],rc[2]);
                
                rc[0],rc[2]=l,r;
                dc=rg['dc'];
                if dc>0:
                    rg['rects_c']=copy.deepcopy(rg['rects']);
                    l-=dc;
                    r+=dc;
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

def region_by_name_update(regions,name,urg={}):
    if type(regions) is dict:
        regions=regions['regions'];    
    if type(name) is dict:
        ns=name.keys();
        for n in ns:
            region_by_name_update(regions,n,name[n])
    else:
        for r in regions:
            if r['name']==name:
                r.update(urg)
                return r;


class MIDdecay0(object):
    def __init__(self,mesh,gcoil='sgcoil',rcoil='srcoil',rcoil_inflate=1e-3):
        print('mesh preparing...')
        mesh=reparse_mesh(mesh);
        #mesh=mesh_region_disabled_off(mesh);
        self.mesh=mr.MID_mesh(mesh);

        self.gcoil_name=gcoil;
        self.rcoil_name=rcoil;
        self.rcoil_inflate=rcoil_inflate;
        #self.gcoil=self.mesh.region_mask(gcoil)
        #self.rcoil=self.mesh.region_mask(rcoil,rcoil_inflate)
        self.fmd=True
        self.perf={}
        self.finitialize=False;
        self.xf=None;
        self.units={'s_ph':1.0};
        print('end')

    def __getitem__(self,name):
        return self.mesh.regions[name];

    def e_circuit(self,jcc,ACV,units={}):
        
        def etrace(A):
            d=np.abs(A.tocoo().diagonal())
            return np.mean(d);

        cCGS=(10**2)*scs.c; 
        CGS9=1e9;       
        c4pi=cCGS/(4*scs.pi);
        sm2m=100.;

        
        ru=CGS9/cCGS**2;

        '''
        #us={'a':1.0,'e':0.5,'r':c4pi}
        vds=2./scs.pi/100.;
        vds=1./100.;
        us={'a':1.0/c4pi,'e':0.5/c4pi,'r':ru*vds}
        '''

        (a,e,r)=(1.0/c4pi,0.5/c4pi,ru/sm2m)

        

        #esp=etrace(ACV[0])        
        #esp=1.0;
        #us={'a':a/r,'e':esp*e,'r':esp*1.0,'s_ph':esp }

        #
        us={'a':a,'e':e,'r':r,'s_ph':1.0 }
        us.update(units);
        self.units=us;
        
        
        jcc=copy.deepcopy(jcc);
        self.coils=jcc;		
        mesh=self.mesh;	

        J_ids={}

        i=0;
        for jc in jcc:
            _id=jc.get('name',jc.get('id','id'+str(i)));
            J_ids[_id]=i;            
            for j in jc['coils']:
                j['mask']=mesh.region_mask(j.get('name',j.get('id',_id)))

            i+=1;
        J_ids['Mz']=i;
        self.J_ids=J_ids
        return ec.coil_chain_t(ACV,us).make_cicuit(jcc);                

    def getBzBy(self):
        if not '_BzBy' in self.__dict__: 
            V=self.mesh.make_volume();
            #V=self.V
            self._BzBy=self.mesh.make_gradxy(V)
            pass
        n=self.NField;
        xx=self.solver.xx[:,0:n];
        yy=mulxx(self._BzBy,xx);
        return yy.conj();

    def pro_init(self,**kws):

        if  self.finitialize:
            return (self.opts,self.A,self.C,self.V);

        print('pro-init...')
        self.solver=None


        opts=extend(kws,{'nd':4,'tic_corr':0,'fmbc':0,'fuzzy':None,'parallel':True,'z_scale':1.0,'pade_nm':(2,4),'e_circuit':None,'e_units':{}})
        
        if self.fmd:
            self.mesh.make_vxs_deformation();

        scheme=opts.pop('scheme',1)
        fuzzy=opts.pop('fuzzy',None)
        printf('FEM -scheme: 0x%X\n',scheme);
        
        (A,C,V,mmz)=self.mesh.makeACVM(scheme=scheme,fuzzy=fuzzy)
        V=self.mesh.make_volume()

        e_c=opts['e_circuit']
        self.fec=fec=not not e_c
        if(fec):
            (A,C,C0,self.jcoil,self.NField)=self.e_circuit(e_c,(A,C,V),opts['e_units'])

            
        
        

        self.fmbc=fmbc=not not opts['fmbc']

        if fmbc:
            (N,tmp)=A.shape;
            mMz=np.zeros(N,dtype = np.float64);
            mMz[0:self.NField]=mmz[:];
            (A,C)=ec.addMz(A,C,C0,self.mesh.vxs,mMz);


        self.C=C;
        self.A=A;
        self.V=V;       


        self.opts=opts;
        self.finitialize=True;
        return (opts,A,C,V);
        
    def init_ms(self,**kws):
        
        printf('init_ms...')
        (opts,A,C,V)=self.pro_init(**kws);
        tt=opts['times'];
        pnm=opts['pade_nm'];
        pade_nm_s=[];
        for t in tt:
            p=pnm if len(t)<3 else t[2];
            pade_nm_s.append(p);
        
        self.times=ts=itimes.steps(tt);
        tic=sw.Tic();
        common=opts.get('common',{});
        self.solver=sg=LIPA_solver_group_ctx((A,C),dt=ts.dts,parallel=opts['parallel'],nd=opts['nd'],pade_nm_s=pade_nm_s,common=common);
        printf('\r LU elapsed %f sec\n',tic.sec());
        for t in ts:
            t.solver=sg.solvers[t.n];
            
    def j_parse(self,J,b):

        tJ=type(J);

        if tJ in (dict,):
            J=(J,);

        if tJ in (tuple,list):
            cJ=[];
            for j in J:
                c=j.copy();
                c['b']=b;
                cJ+=[c]
            return cJ;
        else:
            return J*b;
            
    def j_parse_ec(self,J,NField,Nf):

        
        sb=self.units.get('s_ph',1.0);

        def blank(i=0):
            b=np.zeros(Nf,dtype =np.double);
            b[NField+i]=sb;
            return b

        

        

        tJ=type(J);

        if tJ in (dict,):
            J=(J,);
            tJ=type(J);

        if tJ in (tuple,list):
            
            J_ids=self.J_ids
            cJ=[];
            for j in J:
                if j:
                    c=j.copy();
                    q=c.get('U',c.get('Q',1.0));
                    q=to_list(q);
                    c['Q']=q
                    c['b']=blank(J_ids[j['id']]);
                    cJ+=[c]

            return cJ;
        else:
            return J*blank();

    def decay_ms_ec_mix(self,EMF=[],xinit=None,rcoil=None,xzero=True,eps=1e-4):

        if rcoil is None:
            rcoil=self.rcoil_name;

        rfield=self.mesh.form_grad_xy(self.rcoil_name,xzero=xzero,eps=eps,with_deform=0)[0];


        solver=self.solver;
        xx=self.solver.xx;

        if not xinit is None:
            xx[:]=xinit;        
        

        ts=self.times;
        tt=ts.times;

        nd,Nf=xx.shape;        
        NField=self.NField
        NT=tt.size;
        NJ=Nf-NField;
        rcoil=np.zeros(Nf,dtype=np.float64);
        rcoil[0:NField]=rfield;
        
        #gcoil[NField+1]=-1;


        cEMF=self.j_parse_ec(EMF,NField,Nf);

        solver.reset_J(cEMF);
        
        jjec=np.zeros((NT,nd,NJ),dtype=np.double)

        
        jrs=np.zeros_like(tt);
        jrs1=np.zeros_like(tt);
        jrs2=np.zeros_like(tt);
        
        count=np.size(tt);
        nc=0;
        for t in ts:
            solver  = t.solver;
            Nt = t.nt;
            for nt in range(Nt):
                solver.step();                
                jjet=self.jcoil(xx);
                jjec[nc]=jjet;

                jrs[nc]=solver.form(rcoil,1+0); 
                jrs1[nc]=solver.form(rcoil,1+1);
                jrs2[nc]=solver.form(rcoil,1+2);                                 

                nc+=1;
                printf('\rec: %d %d [%d %d]',nc,count,nt+1,Nt);
                #printf('\n %d %d [%d %d] t=%f g=%e r=%e',nc,count,nt+1,Nt,tt[nc-1]*1000,jg,jr);
            printf('\n') ;  

        printf('end \n') ; 
        
        jj=np.empty((3,np.size(tt)));
        jj[0][:]=jrs;
        jj[1][:]=jrs1;
        jj[2][:]=jrs2;


        res_ec={
                't':tt,
                'jj':np.transpose(jjec),
                'ids':self.J_ids.copy()
                };


        res_rcoil={
                't':tt,
                'j':jrs,
                'j1':jrs1,
                'j2':jrs2,
                'jj':jj,
                'ec':res_ec
                };

                
        return res_rcoil; 

    def decay_ms_ec2(self,EMF=[],xinit=None):
        
        solver=self.solver;
        xx=self.solver.xx;

        if not xinit is None:
            xx[:]=xinit;        
        

        ts=self.times;
        tt=ts.times;

        nd,Nf=xx.shape;        
        NField=self.NField
        NT=tt.size;
        NJ=Nf-NField;
        
        #gcoil[NField+1]=-1;


        cEMF=self.j_parse_ec(EMF,NField,Nf);

        solver.reset_J(cEMF);
        
        jj=np.zeros((NT,nd,NJ),dtype=np.double)

        
        
        #print('start...2')
        count=np.size(tt);
        nc=0;
        for t in ts:
            solver  = t.solver;
            Nt = t.nt;
            for nt in range(Nt):
                
                solver.step();
                
                jjet=self.jcoil(xx);
                jj[nc]=jjet;
                nc+=1;
                printf('\rec: %d %d [%d %d]',nc,count,nt+1,Nt);
                #printf('\n %d %d [%d %d] t=%f g=%e r=%e',nc,count,nt+1,Nt,tt[nc-1]*1000,jg,jr);
            printf('\n') ;  

        #printf('end \n') ; 
        
        

        res={
                't':tt,
                'jj':np.transpose(jj),
                'ids':self.J_ids.copy()
                };
                
        return res; 
            
    def decay_ms_ec(self,EMF=1.0,xinit=None):
        solver=self.solver;
        xx=self.solver.xx;

        if not xinit is None:
            xx[:]=xinit;        
        print(xx)

        ts=self.times;

        Nf=self.A.shape[0];
        NField=self.NField
        gcoil=np.zeros([Nf],dtype = np.float64);
        gcoil[NField]=1;
        #gcoil[NField+1]=-1;
        cEMF=self.j_parse(EMF,gcoil);
        solver.reset_J(cEMF);
        
        tt=np.array(ts.tt,dtype=np.double);
        jrs=np.zeros_like(tt);
        jrs1=np.zeros_like(tt);
        jrs2=np.zeros_like(tt);

        jgs=np.zeros_like(tt);
        jgs1=np.zeros_like(tt);
        jgs2=np.zeros_like(tt);
        
        
        count=np.size(tt);
        nc=0;
        for t in ts:
            solver  = t.solver;
            Nt = t.nt;
            for nt in range(Nt):
                solver.step();
                a=xx[0,0:NField];
                na=norm(a);
                jjet=self.jcoil(xx);
                jg=jjet[0,0]
                jr=jjet[0,1]

                jgs[nc]=jjet[0,0]; 
                jgs1[nc]=jjet[1,0];
                jgs2[nc]=jjet[2,0];                                 


                jrs[nc]=jjet[0,1]; 
                jrs1[nc]=jjet[1,1];
                jrs2[nc]=jjet[2,1];                                 
                nc+=1;
                printf('\rec: %d %d [%d %d]',nc,count,nt+1,Nt);
                #printf('\n %d %d [%d %d] t=%f g=%e r=%e',nc,count,nt+1,Nt,tt[nc-1]*1000,jg,jr);
            printf('\n') ;  

        printf('end \n') ; 
        
        jj=np.empty((3,np.size(tt)));
        jj[0][:]=jrs;
        jj[1][:]=jrs1;
        jj[2][:]=jrs2;

        res={
                't':tt,
                'j':jrs,
                'j1':jrs1,
                'j2':jrs2,
                'jg':jgs,
                'jg1':jgs1,
                'jg2':jgs2,
                'jj':jj
                };
                
        return res; 
        
            
    def decay_ms(self,J=1.0,xinit=0.0):
        solver=self.solver;
        #C=self.C.tocsc();

        #gcoil=C*self.gcoil.astype('d');
        #rcoil=C*self.rcoil.astype('d');

        gcoil=self.mesh.form_volume(self.gcoil_name,eps=1e-4,scheme=1,with_deform=0);
        rcoil=self.mesh.form_volume(self.rcoil_name,eps=1e-4,scheme=1,with_deform=0);

        if not xinit is None:
            self.solver.xx[:]=xinit;

        ts=self.times;
        cJ=self.j_parse(J,gcoil);
        solver.reset_J(cJ);
        
        tt=np.array(ts.tt,dtype=np.double);
        jrs=np.zeros_like(tt);
        jrs1=np.zeros_like(tt);
        jrs2=np.zeros_like(tt);
        count=np.size(tt);
        nc=0;
        for t in ts:
            solver  = t.solver;
            Nt = t.nt;
            for nt in range(Nt):
                solver.step();
                jrs[nc]=solver.form(rcoil,1+0); 
                jrs1[nc]=solver.form(rcoil,1+1);
                jrs2[nc]=solver.form(rcoil,1+2);                                 
                nc+=1;
                printf('\r %d %d [%d %d]',nc,count,nt+1,Nt);
            printf('\n') ;  

        printf('end \n') ; 
        
        jj=np.empty((3,np.size(tt)));
        jj[0][:]=jrs;
        jj[1][:]=jrs1;
        jj[2][:]=jrs2;

        res={
                't':tt,
                'j':jrs,
                'j1':jrs1,
                'j2':jrs2,
                'jj':jj
                };
                
        return res; 
            

    def decay_ms2(self,J=1.0,xinit=0.0,xzero=True):
        solver=self.solver;
        #C=self.C.tocsc();

        #gcoil=C*self.gcoil.astype('d');
        #rcoil=C*self.rcoil.astype('d');

        gcoil=self.mesh.form_volume(self.gcoil_name,eps=1e-4,scheme=1,with_deform=0);
        rcoil=self.mesh.form_grad_xy(self.rcoil_name,eps=1e-4,with_deform=0)[0];

        if not xinit is None:
            self.solver.xx[:]=xinit;

        ts=self.times;
        cJ=self.j_parse(J,gcoil);
        solver.reset_J(cJ);
        
        tt=np.array(ts.tt,dtype=np.double);
        jrs=np.zeros_like(tt);
        jrs1=np.zeros_like(tt);
        jrs2=np.zeros_like(tt);
        count=np.size(tt);
        nc=0;
        for t in ts:
            solver  = t.solver;
            Nt = t.nt;
            for nt in range(Nt):
                solver.step();
                jrs[nc]=solver.form(rcoil,1+0); 
                jrs1[nc]=solver.form(rcoil,1+1);
                jrs2[nc]=solver.form(rcoil,1+2);                                 
                nc+=1;
                printf('\r %d %d [%d %d]',nc,count,nt+1,Nt);
            printf('\n') ;  

        printf('end \n') ; 
        
        jj=np.empty((3,np.size(tt)));
        jj[0][:]=jrs;
        jj[1][:]=jrs1;
        jj[2][:]=jrs2;

        res={
                't':tt,
                'j':jrs,
                'j1':jrs1,
                'j2':jrs2,
                'jj':jj
                };
                
        return res; 
                
        
        
        
        
    def init_z(self,**kws):
        self.solver=None
        print('init_z...')
        (opts,A,C,V)=self.pro_init(**kws);
        z_scale=opts['z_scale'];
        zs=to_list(opts['z']);
        zs=[z_scale*z for z in zs];
        self.zs=zs;
        self.solver=LIPA_solver_ctx_z((A,C),z=zs,parallel=opts['parallel']);

    def init_f(self,**kws):
        self.solver=None
        print('init_z...')
        (opts,A,C,V)=self.pro_init(**kws);
        jpi=np.pi*1j;
        zs=[jpi*f for f in opts['fs']];
        self.zs=zs;
        self.solver=LIPA_solver_ctx_z((A,C),z=zs,parallel=opts['parallel']);


    def spectrum(self,gcoil='lgcoil',rcoil='srcoil',ndiff=1,raw_fields=False):
        xf=self.xf;
        eps=self.rcoil_inflate;
        jgcoil=self.mesh.form_volume(gcoil,eps=eps,scheme=1,with_deform=0);
        jrcoil=self.mesh.form_volume(rcoil,eps=eps,scheme=1,with_deform=0);
        self.xf=xf=self.solver.make(jgcoil,xf);
        zs=self.zs;
        c=len(xf);

        if raw_fields:
            if ndiff:
                i=iter(zs);
                r=[ np.power(next(i),ndiff)*r for r in xf];
                return r;
            else:
                return xf;

        
        r=np.zeros(c,dtype=np.complex128);
        
        for k in range(c):
            jr=np.power(zs[k],ndiff)*jrcoil;
            r[k]=np.dot(jr,xf[k]);
        return r;



    def init(self,**kws):
        print('init...')
        tic=sw.Tic();
        
        (opts,A,C,V)=self.pro_init(**kws);

        #(A,C)=(A.tocsc(),C.tocsc())
        #A.data=-A.data;
        tic.start()
        
        self.dt=opts.get('dt',1.0);
        LIPA_solver=cast_str(opts.pop('LIPA_solver',LIPA_solver_ctx));
        lu_factory=cast_str(opts.pop('lu_factory',lu_solver_factory()));        
        print(LIPA_solver)    

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

    def decay(self,count=100,J=1.0,fnparray=False,fctx=0,ftnull=True):
        solver=self.solver;
        C=self.C.tocsc();

        #gcoil=C*self.gcoil.astype('d');
        #rcoil=C*self.rcoil.astype('d');

        gcoil=self.mesh.form_volume(self.gcoil_name,eps=1e-4,scheme=1,with_deform=0);
        rcoil=self.mesh.form_volume(self.rcoil_name,eps=1e-4,scheme=1,with_deform=0);


        solver.reset_J(J*gcoil);
        
        jrs=[0.0] if ftnull else [];
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
                printf('\r %d of %d ',n+1,count);
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
        #gcoil=self.gcoil.astype('d');
        #rcoil=self.rcoil.astype('d');
        field_step=others.get('field_step',-1);
        field_step=int(field_step);


        #jgcoil=J*(C*gcoil);
        #jrcoil=C*rcoil;
        eps=self.rcoil_inflate

        jgcoil=self.mesh.form_volume(self.gcoil_name,eps=eps,scheme=1,with_deform=0);
        jrcoil=self.mesh.form_volume(self.rcoil_name,eps=eps,scheme=1,with_deform=0);

        solver.reset_J(jgcoil);

        jfactor=1.0;

        def jsol():
            jf=-(A*solver.xn[0]) +(jfactor)*jgcoil;
            return (np.dot(jrcoil,jf),jf);


        
        
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


def repair_regions(mesh):
    def repair_rect(r):
        if r[0]>r[2]:
            (r[0],r[2])=(r[2],r[0]);
        if r[1]>r[3]:
            (r[1],r[3])=(r[3],r[1]);

    for rg in mesh['regions']:
        for r in rg['rects']:
            repair_rect(r);
    pass


class MID3barrier_decay(MIDdecay0):
    def __init__(self,mesh,gcoil='gcoil',rcoil='rcoil',rcoil_inflate=1e-3,f_repair_regions=False):
        mesh=reparse_mesh(mesh)

        if f_repair_regions:
            repair_regions(mesh);

        #mesh=mesh_region_disabled_off(mesh);
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
        regions_new=copy.copy(self.mesh.regions.regions);

        if kws.get('clsdis',False):
            rs,regions_new=regions_new,[];
            for r in rs:
                if not r.get('disabled',0):
                    regions_new+=[r];


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
    
    def init_ms(self,**kws):
        #barriers=kws.get('barriers',{});
        #update_nominal_barriers(self.mesh.regions.regions,barriers)
        
        if kws.get('barriers',{}):
            self.mesh_reconstruct(**kws);
        super(MID3barrier_decay,self).init_ms(**kws);
        return self;


    def init_z(self,**kws):                
        if kws.get('barriers',{}):
            self.mesh_reconstruct(**kws);
        super(MID3barrier_decay,self).init_z(**kws);
        return self;

    def init_f(self,**kws):
                
        if kws.get('barriers',{}):
            self.mesh_reconstruct(**kws);
        super(MID3barrier_decay,self).init_f(**kws);
        return self;

    def reset_field(self,x=0):
        self.solver.x=x;


