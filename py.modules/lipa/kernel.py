
import numpy as np
#import ctypes
import os;
from scipy import sparse as sp
from scipy.sparse import csc_matrix
from scipy.sparse import  linalg as sla
from asyn.ipc_matrix import ipc_array,ipc_matrix,schema,ipc_unpack
import numbers
import copy
#import lipa.qpl as qpl
import asyn.SharedWorker as sw
from utils.printf import printf

from .qpl import *
from p23 import *

from . import pade_exp_poles_res
#import lipa.pade_exp_poles_res as pade_exp_poles_res

norm=np.linalg.norm;
Tic=sw.Tic;



def zero_like(C,dtype=None):
    t= dtype if dtype else C.dtype;
    return csc_matrix(C.shape,dtype = t);

def eye_like(C,dtype=None):
    t= dtype if dtype else C.dtype;
    return sp.identity(C.shape[0],dtype =t);

def spmul(a,N,z=1):
    if isinstance(a,numbers.Number):
        return sp.identity(N,dtype =np.complex128).multiply(a*z);
    else:
        return a.multiply(z);

def calc_Az(AC,zt,Bt):
    
    bbt=-np.complex128(1)/Bt;
    zn=zt*bbt;    
    Az=AC[0].multiply(bbt)
    #.astype(np.complex128);
    N=Az.shape[0];
    C=AC[1:]
    if C:
        for c in C:
            Az+=spmul(c,N,zn)#zn*c;
            zn*=zt;
    else:
        Az+=zn*eye_like(Az)
    return (Az,C,N);


class lus:
    def __init__(self,solve,A=None):
        self.solve=solve
        self.A=A
    def step(self,xin,xout):
        xout[:]=self.solve(xin);
        return xout;

class lus_r:
    def __init__(self,solve):
        self.solve=solve
    def step(self,xin,xout):
        xout[:]=self.solve(xin.astype('d')).astype('complex128');
        return xout;


class lu_solver_factory:
    def __init__(self,options={}):
        #self.options=extend(options,dict(Equil=False, IterRefine='SINGLE'))
        #self.options=extend(options,dict(Equil=False, IterRefine='NO'))
        options=copy.copy(options);
        #options=dict(Equil=False, IterRefine='NO');
        #options=dict(Equil=False, IterRefine='SINGLE');
        #options={'permc_spec':'MMD_AT_PLUS_A'};#'MMD_ATA'
        #options={'permc_spec':'COLAMD'};        
        #options=dict(SymmetricMode=True);
        self.options=options;

    def __call__(self,Deltaz):
        #lu=sla.splu(Deltaz,permc_spec='COLAMD',options=self.options);
        #Deltaz=Deltaz._real().tocsc();
        #Deltaz=Deltaz.astype('d');
        lu=sla.splu(Deltaz,permc_spec=None,options=self.options);
        #return lus_r(lu.solve).step;
        return lus(lu.solve,Deltaz).step;

class bicgstab_s:
    def __init__(self,so):
        self.sargs=so;
    def step(self,xin,xout):
        
        (y,err)=sla.bicgstab(b=xin,**self.sargs);
        self.err=err;
        if err:
            raise Exception('bicgstab error='+str(err));
        xout[:]=y;
        return xout;


class bicgstab_solver_factory:
    def __init__(self,options={}):
        print('step....')
        options=copy.copy(options);
        self.options=options;
        so={
            'A':None,             
            'x0':None,
            'tol':1e-10,
            'maxiter':None,
            'xtype':None,
            'M':None,
            'callback':None,
            }
        self.sargs=so
        self.err=0;
        print(so);
        
    def __call__(self,Deltaz):
        
        #lu=sla.splu(Deltaz,permc_spec='COLAMD',options=self.options);
        #Deltaz=Deltaz._real().tocsc();
        #Deltaz=Deltaz.astype('d');
        so=copy.copy(self.sargs);
        so['A']=Deltaz;
        self.A=Deltaz;
        if so['M']=='d':
            d=Deltaz.diagonal();
            d=1.0/np.abs(d);
            M=sp.diags([d],[0]);
            so['M']=M;
        return bicgstab_s(so).step;


class LIPA_polus_base(object):
    def __init__(self,polus_res,dt,AC,**opts):    
    #def __init__(self,polus_res,dt,AC,SolverFactory=lu_solver_factory(),**opts):
        self.dt=dt=np.double(dt);
        bdt=np.double(1)/dt;
        zt=polus_res[0]*bdt
        self.zt=zt
        self.Bt=Bt=bdt*polus_res[1]
        self.Bt_zt=polus_res[1]/polus_res[0];
        (Az,self.C,N)=calc_Az(AC,zt,Bt);

        #az=Az.todense();
       # print('resz',Bt,'Az/Bt',az,'Az',-az*Bt)

        self.N=N
        #self.solver_step=SolverFactory(Az);
        self.qpj=qpJ(dt=dt);
        self.x_base=np.empty(N,np.complex128);
        self.xbz=np.empty(N,np.complex128);
        self.tmp_x=np.empty(N,np.complex128);
        self.Az=Az;
        self.ncorr=0;    
        fCp=opts.get('fCp',True);
        self.calc_Cz=calc_Cx1 if fCp else calc_Cx;
        
    def correction(self,x,b):

        def renorm(c):
            err=np.linalg.norm(c,1)/np.size(c);
            if err>0:
                s=np.log2(err);
                n=int(s);
                if n<0:
                    n2=2<<(-n)
                    n2=np.double(n2);
                    c/=n2;
                    return (n2,err);         
            return (1.0,err);
             
        r=b-self.Az*x;
        (sc,err)=renorm(r);        
        c=self.solver_step(r,self.tmp_x);
        if sc!=1:
            c*=sc;
        x+=c;
        r=b-self.Az*x;
        (sc,err)=renorm(r);        

        return (sc,err)

    def reset_corr(self,n):
        (self.ncorr,n)=(n,self.ncorr)
        return n;

    def step(self,xbz,nd,xx):
        z=zt=self.zt;
        x0=self.solver_step(xbz,xx[0])

        for j in range(self.ncorr):
            (sc,err)=self.correction(x0,xbz);


        x_base=self.x_base#*0;

        for k in range(1,nd):
            g=zt*x0+x_base[k-1];
            xx[k][:]=g;
            x0=g;
        '''
        sg=1.;
        #x_base[:]=0;
        xbg=np.zeros(self.N,np.complex128);
                
        for k in range(1,nd):
            
            #xx[k][:]=z*(x0+sg*x_base[k-1]);
            g=x_base[k-1];            
            xx[k][:]=z*x0+g+xbg;
            xbg+=sg*z*g;                                
            #xx[k][:]=z*(x0+x_base);
            z*=zt;        
        '''
        return xx;

            
    def reset_J(self,J):
        self.qpj=qpJ(J,self.dt);
    
    def set_x_base(self,v): 
        x_base=self.x_base;              
        if v.size>x_base.size:
            x_base=self.x_base=np.empty(shape=v.shape,dtype=np.complex128);            
        x_base[:]=self.Bt*v;
        #x_base[-1,:]=0;
        

    def setJ(self,xbz):
        self.qpj.set_J(self.zt,xbz);

    def set_xbz(self,xx,x0=0):#,x0_base=0):
        zt=self.zt
        xbz=self.xbz;
        #self.x_base[:]=self.Bt*x_base;
        #self.x_base[:]=self.Bt_zt*xx[0];
        #self.set_x_base(self.Bt_zt*xx);
        self.set_x_base(xx);
        xbz[:]=x0;
        #self.qpj.set_J(zt,xbz)
        self.setJ(xbz);
        self.calc_Cz(xbz,zt,self.C,xx,self.tmp_x)
        #calc_Cx1(xc,z,CC,xx):
        return xbz;

class LIPA_polus(LIPA_polus_base):
    def __init__(self,polus_res,dt,AC,SolverFactory=lu_solver_factory(),**opts):
        super(LIPA_polus,self).__init__(polus_res,dt,AC,**opts);
        self.solver_step=SolverFactory(self.Az);

class LIPA_polus_delta(LIPA_polus_base):
    def __init__(self,polus_res,dt,AC,dAC,solver_step,**opts):
        super(LIPA_polus_delta,self).__init__(polus_res,dt,AC,**opts);
        self.dAC=dAC;
        self.solver_step=solver_step;
        (zt,Bt)=(self.zt,self.Bt)
        (self.dAz,self.dC,N)=calc_Az(dAC,zt,Bt);

    def setJ(self,xbz):
        xxz,xx0=self.xxz_xx0;
        (zt,xbz,dAz,dC)=(self.zt,self.xbz,self.dAz,self.dC);        
        self.xbz[:]=-dAz*xx0;
        self.calc_Cz(xbz,zt,dC,xxz,self.tmp_x);


    def set_xxz_xx0(self,xxz_xx0):
        self.xxz_xx0=xxz_xx0;
        #self.calc_Cz(self.xbz,self.zt,self.C,xxz,self.tmp_x)


class LIPA_poles_thread(object):
    def __init__(self,polus_res_list,dt,AC,xx0,xc0,nd,xx,SolverFactory=lu_solver_factory(),**opts):

        AC=[convert_tocsc(a) for a in AC];
        if 0:
            for a in AC:
                a.eliminate_zeros();
                    
        self.solvers=[LIPA_polus(p_r,dt,AC,SolverFactory,**opts) for p_r in polus_res_list]
        #print('solvers ok pid:',os.getpid())
        self.xc0=xc0;
        self.xx0=xx0;
        self.xx=xx;
        self.nd=nd;
        #print('init ok pid:',os.getpid(),'type(xx)', type(xx))
        self.freal=freal=(xx.dtype==np.float64) or (xx.dtype==np.float32);
        #print('next init ok pid:',os.getpid(),'type(xx)', type(xx))
        #self.xc0_c=np.empty_like(xc0,dtype=np.complex128)
        self.xx_c=np.empty_like(xx,dtype=np.complex128)
        
        if freal:
            prj =lambda x: x.real
        else:
            prj =lambda x: x
        self.prj=prj
    
    def reset_J(self,J):
        for solver in self.solvers:
            solver.reset_J(J)

    def reset_corr(self,n=0):
        for solver in self.solvers:
            solver.reset_corr(n);

    def step(self):
        #xc0_c=self.xc0_c;
        #xc0_c[:]=self.xc0;
        xc0=self.xc0;
        xx=self.xx;
        xx0=self.xx0;
        xx_c=self.xx_c;
        prj=self.prj;
        f=False;
        nd=self.nd
        #testp=np.zeros(nd);
        
        for solver in self.solvers:
            solver.set_xbz(xx0,xc0)#xx[0])
        
        for solver in self.solvers:            
            #solver.set_xbz(xx0,xc0);
            xz=solver.xbz;
            solver.step(xz,nd,xx_c);
            #testp+=prj(solver.testbz);
            if f:
                xx[:]+= prj(xx_c)
            else:
                xx[:]= prj(xx_c)
                f=True;
        #self.testp=testp;

class LIPA_forms(object):
    def __init__(self,vxn):
        self.vxn=vxn;
    @property
    def xn(self):
        return self.vxn;
    @xn.setter
    def xn(self,v):
        self.vxn[:]=v;
    @property
    def x(self):
        return self.vxn[0];
    @x.setter
    def x(self,v):
        self.vxn[:]=0;
        self.vxn[0]=v;
    def form(self,lf,dn=0):
        v=self.vxn[dn];
        return np.dot(v,lf);
    def zeros(self):
        return np.zeros_like(self.vxn[0]);

def convert_tocsc(a):
    if a is None:
        return a;
    if type(a) in [list,tuple]:
        return [convert_tocsc(i) for i in a]
    elif isinstance(a,numbers.Number):
        return csc_matrix([a])
    elif type(a)==np.ndarray:
        return csc_matrix(a)
    else:
        return a.tocsc();
    



class LIPA_solver_st(LIPA_poles_thread,LIPA_forms):
    #def __init__(self,AC,dt=1,nd=1,fcomplex=False,SolverFactory=lu_solver_factory(),**opts):
    def __init__(self,AC,SolverFactory=lu_solver_factory(),**opts):
        
        tic=Tic(); 
        AC=tolist(AC)
        #AC=[convert_tocsc(i) for i in AC];
        AC=convert_tocsc(AC);
        fcomplex=opts.get('fcomplex',False)
        pade_nm=opts.get('pade_nm',(2,4))
        dAC=opts.get('dAC',[]);
        #dAC=[convert_tocsc(i) for i in dAC];

        dAC=convert_tocsc(dAC);

        fCp=(not bool(dAC)) and opts.get('fCp',True);

        
        #print('fCp',fCp);

        fp=pade_nm[0]<pade_nm[1];

        self.flag_pade=fp;



        dt=opts.get('dt',np.double(1))
        nd=opts.get('nd',1)
        nd=max(len(AC)-1,1,nd);
        
        polus_res_list=pade_exp_poles_res.get(pade_nm,not fcomplex)
        N=AC[0].shape[0];        
        #float_type=np.float64 if not fcomplex else np.complex128;
        float_type= np.complex128 if fcomplex else np.float64;
        
        saturn=opts.get('ipc_saturn',None)
        ffinished=saturn is None;
        if not fp:
            self.xx00=np.zeros((nd,N),dtype=float_type);
            self.step=self._step_fp
        else:
            self.step=self._step

        if ffinished:
            xc0=np.zeros(N,dtype=float_type)
            xx=np.zeros((nd,N),dtype=float_type)
            
        else:
            class fake_group:
                pass
            
            self.group=fake_group()
            self.group.result=True;
            (xc0,xx)=(saturn[0].value,saturn[1].value);
        xx0=xx;
        super(LIPA_solver_st,self).__init__(polus_res_list,dt,AC,xx0,xc0,nd,xx,SolverFactory,**{'dAC':dAC,'fCp':fCp})
        LIPA_forms.__init__(self,xx0);
        self.C=AC[1:]
        self.tinit=tic.sec();
        self.tloading=0;  
        self.dAC=dAC;
        self.calc_Cx0=calc_Cx0 if fCp else calc_Cx0_fake;
        
    def step_t(self,ncount=1):
        self.step(ncount,True)
        t=self.tStep
        printf('\r tic=%f sec                ',t)
        return t;
        
    def _step_fp(self,ncount=1,ftic=True):
        if ftic:
            tic=Tic();        
        xx0=self.xx0        
        xx00=self.xx00

        cCx0=self.calc_Cx0;
                
        for c in range(ncount):
            xx00[:]=xx0;
            
            cCx0(self.xc0,self.C,xx0);
            super(LIPA_solver_st,self).step();            
            #
            xx0+=xx00;
            #xx0[0:2,:]+=xx00[0:2,:];
            #xx0[2:,:]=xx00[2:,:];
            #xx0[2:,:]=0;
            
        if ftic:
            self.tStep=tic.sec()/ncount;
        return self.xx0;

    def _step(self,ncount=1,ftic=True):
        if ftic:
            tic=Tic();     
                    
        cCx0=self.calc_Cx0;

        for c in range(ncount):
            cCx0(self.xc0,self.C,self.xx0);
            super(LIPA_solver_st,self).step();  
                 
        if ftic:
            self.tStep=tic.sec()/ncount;
        return self.xx0;
        """
    @property
    def xn(self):
        return self.xx0;
    @xn.setter
    def xn(self,v):
        self.xx0[:]=v;
    @property
    def x(self):
        return self.xx0[0];
    @x.setter
    def x(self,v):
        self.xx0[:]=0;
        self.xx0[0]=v;
        """
    def reset_J(self,J):
        if type(J)==np.ndarray:
            J=({'qp':({'c':(1,)},),'x':J.flatten()},);
        super(LIPA_solver_st,self).reset_J(J);













"""
class lu_poles_node:
    @staticmethod
    @staticmethod
    def calc_poly(C,z):
        zn=np.complex128(1);
        cz=self.zero_like(C[0],np.complex128);
        for c in C:
            cz+=c*zn;
            zn*=z;
        return cz;

    def reset_J(self,J):
        self.qpj=qpJ(J,self.dt);

#,AC,polus,res,xz,c0,cz,freal=True,lu_solver_class=lu_polus):
    def __init__(self,opts):
        opts=extend(opts,{'freal':True,'dt':1,'la_solver':{'class':lu_polus,'options':{}}});
        self.dt=opts['dt'];
        AC=opts['AC'];
        z=self.polus=polus=opts['polus'];
        res=opts['res'];
        
        A=AC[0];
        C=AC[1:];
        
        if len(C)==0:
            Cz=self.eye_like(C,np.complex128);
        else:
            Cz=self.calc_poly(C,z);
        
        br=np.complex128(1)/res;
        Az =br*(A+z*Cz);
        lu_solver_class=opts['la_solver']['class']
        lu_solver_opts=opts['la_solver'].get('options',{})
        self.lu=lu_solver_class(Az,lu_solver_opts);
        self.Cz = Cz;
        self.polus=polus;
        data=opts['data']
        xz=data['xz'];
        cz=data.get('cz',xz);
        c0=data['c0'];
        xb=self.zero_like(cz,np.complex128);
        (self.xz,self.c0,self.cz,self.xb)=(xz,c0,cz,xb);

    def step(self):
        
        self.lu.step(self,xb,xz);
        if not self.xz is self.cz:
            self.cz[:]=self.Cz*xz;
        return True;


class lipa_solver_real_st:
    def __reset_AC_x0(self,x0,AC,l):
        AC=tolist(AC);
        N=AC[0].shape[0];
        czs=[];xzs=[];
        self.flc=flc=len(AC)>1;
        if k in range(l):
            z=np.zeros(N,dtype=np.complex128)
            czs.append(z)
            if flc:
                z=np.zeros(N,dtype=np.complex128)
            xzs.append(z)
        c0=np.zeros(N,dtype=self.xtype);
        if x0==None:
            x0=np.zeros(N,dtype=self.xtype);
        else:
            lx0=tolist(x0)
            M=min(len(x0s),len(Ac)-1);
            for m in xrange(0,M):
                c0+=AC[m+1]*lx0[m];
        return (Ac,x0,c0,czs,xzs);

    def __init__(self,n,pade,AC,dt=1,x0=None,la_solver={'class':lu_polus,'options':{}}):
        l=pade.count;
        self.freal=freal=pade.fhalf;
        self.xtype= np.float64 if freal else np.complex128
        seltf.dt=dt;
        zs=pade.poles;
        rs=pade.res;
        opts=extend({'AC':AC,'freal':freal,'dt':dt,'la_solver':la_solver});
        nodes=[];
        for k in xrange(l):
            opts['polus']=zs[k];
            opts['res']=rs[k];
            nodes.append(lu_poles_node(opts));
        self.nodes=nodes;
    def reset_J(self,J):
        for nd in self.nodes:
            nd.reset_J(J)
        
"""
