#
import os
import numpy as np
from asyn.ipc_matrix import ipc_array,ipc_matrix,schema,ipc_unpack
#from asyn.SharedWorker import SharedWorker,AsynGroup,Tic
import asyn.SharedWorker as sw
from . import pade_exp_poles_res
from . import kernel
from .kernel import LIPA_solver_st,convert_tocsc,lu_solver_factory

from .qpl import *

#from . import lipa
#l=lipa.LIPA_poles_thread
#kernel.LIPA_poles_thread

Tic=sw.Tic;


class LIPA_poles_thread_pp(kernel.LIPA_poles_thread):
    def __init__(self,polus_res_list,dt,AC,xx0,xc0,nd,xxp,np,Jsbuf,SolverFactory):
        pass
        AC=ipc_unpack(AC);
        xc0=ipc_unpack(xc0);
        xxp=ipc_unpack(xxp);
        xx0=ipc_unpack(xx0);
        xx=xxp[np];
        #print('create: LIPA_poles_thread_pp pid:',os.getpid())
        tic=Tic()
        super(LIPA_poles_thread_pp,self).__init__(polus_res_list,dt,AC,xx0,xc0,nd,xx,SolverFactory);
        #print('Ok:',os.getpid())
        self.Jsbuf=ipc_unpack(Jsbuf);
        self.t_factor=tic.sec();


    def reset_J(self,J):
        if J:
            Jsbuf=self.Jsbuf
            J=[{'qp':j['qp'],'x':Jsbuf[j['x']]} for j in J];
        super(LIPA_poles_thread_pp,self).reset_J(J);


class params(object):
    pass

class LIPA_poles_thread_pp2(object):
    def __init__(self,polus_res_list,dt,AC,xx0,xc0,nd,xxp,np,Jsbuf,SolverFactory):
        pass
        p=params();
        (p.polus_res_list,p.dt,p.AC,p.xx0,p.xc0,p.nd,p.xxp,p.np,p.Jsbuf,p.SolverFactory)=(polus_res_list,dt,AC,xx0,xc0,nd,xxp,np,Jsbuf,SolverFactory);
        self.params=p;

    def init(self):
        
        p=self.params;
        AC=ipc_unpack(p.AC);
        xc0=ipc_unpack(p.xc0);
        xxp=ipc_unpack(p.xxp);
        xx0=ipc_unpack(p.xx0);
        xx=xxp[p.np];
        #print('create: LIPA_poles_thread_pp pid:',os.getpid())
        tic=Tic()
        self.lipa_poles_thread=kernel.LIPA_poles_thread(p.polus_res_list,p.dt,AC,xx0,xc0,p.nd,xx,p.SolverFactory);
        #print('Ok:',os.getpid())
        self.Jsbuf=ipc_unpack(p.Jsbuf);
        self.t_factor=tic.sec();


    def reset_J(self,J):
        if J:
            Jsbuf=self.Jsbuf
            J=[{'qp':j['qp'],'x':Jsbuf[j['x']]} for j in J];
        self.lipa_poles_thread.reset_J(J);
    def reset_corr(self,n=0):
        self.lipa_poles_thread.reset_corr(n);


    def step(self):
        self.lipa_poles_thread.step();




class LIPA_solver_pp(kernel.LIPA_forms):
    def __init__(self,AC,SolverFactory=lu_solver_factory(),**opts):
        tic=sw.Tic();
        AC=tolist(AC)
        #AC=[convert_tocsc(i) for i in AC];
        pade_nm=opts.get('pade_nm',(2,4))
        fcomplex=opts.get('fcomplex',False)
        fsmd=opts.get('fsmd',False)
        cpus=opts.get('cpus',sw.affinity2cpus())
        dt=opts.get('dt',np.double(1))
        nd=opts.get('nd',1)
        fda=opts.get('affinity_disable',False);
        self.nJsmax=nJsmax=opts.get('nJs',1);
        nd=max(len(AC)-1,1,nd);
        N=AC[0].shape[0];
        fp=pade_nm[0]<pade_nm[1];
        self.flag_pade=fp;
        self.tStep=-1;
        

#    (self,pade_nm,AC,dt=1,nd=1,affinity=0,fcomplex=False,SolverFactory=kernel.lu_solver_factory(),**opts):

        float_type= np.complex128 if fcomplex else np.float64;

        if not fp:
            self.xx00=np.zeros((nd,N),dtype=float_type);
        
        saturn=opts.get('ipc_saturn',None)
        self.ffinished=ffinished=saturn is None;
        
        if ffinished:
            ipcxc0=ipc_array(schema((N,),dtype=float_type))
            #ipcxx=ipc_array(schema((nd,N),dtype=float_type))
            ipcxx0=ipc_array(schema((nd,N),dtype=float_type))
            ipcxJs=ipc_array(schema((nJsmax,N),dtype=float_type))
        else:
            (ipcxc0,ipcxx0,ipcxJs)=saturn;
        
        
        self.xc0=ipcxc0.value;
        self.xx0=xx0=ipcxx0.value;
        self.xJs=ipcxJs.value;
        
        ipcAC=[ipc_matrix(a,fsmd) for a in AC];
        ipcC0=ipcAC[1:2];
        #cc=ipcC0[0].value;
        self.C0=[convert_tocsc(c.value) for c in ipcC0];
        #self.C0=[convert_tocsc(c) for c in ipcC0];
        poles_res=pade_exp_poles_res.get(pade_nm,not fcomplex);
        packets=sw.thread_packets(poles_res,cpus);
        self.packets=packets;
        self.Np=Np=len(packets);
        ipcxxp=ipc_array(schema((Np,nd,N),dtype=float_type));
        self.xxp=ipcxxp.value;
        super(LIPA_solver_pp,self).__init__(xx0);
        #kernel.LIPA_forms.__init__(self,self.xx0);
        """
        self.gg=[]
        for n in range(Np):
             pt = packets[n]
             g=LIPA_poles_thread_pp(pt[0],dt,ipcAC,ipcxx0,ipcxc0,nd,ipcxxp,n,ipcxJs,SolverFactory)
             gg.append(g)
        """
        self.job=job = sw.affinity.create_job() if sw.f_affinity else None
        tic.start()
        """
        pt = packets[0]
        am= 0 if fda else pt[1];
        ll=LIPA_poles_thread_pp2(pt[0],dt,ipcAC,ipcxx0,ipcxc0,nd,ipcxxp,0,ipcxJs,SolverFactory);
        ll.init(0);
        """

        self.group=group=sw.AsynGroup();
        for n in range(Np):
             pt = packets[n]
             am= 0 if fda else pt[1];
             group<<sw.SharedWorker(LIPA_poles_thread_pp2,affinity_mask=am,job=job)\
             (pt[0],dt,ipcAC,ipcxx0,ipcxc0,nd,ipcxxp,n,ipcxJs,SolverFactory)#\
             #.register_methods(('step','reset_J'));
        
        self.tloading=tic.sec();

        if ffinished:
            group.results

        
        tic.start()
        if ffinished:
            for g in group:
                group<<g.call('init');
            group.results;
        
        self.tinit=tic.sec();

    def reset_corr(self,n=0):
        group=self.group;

        for g in group:
            group<<g.call('reset_corr',n);
        group.results;

    def step(self,ncount=1,ftic=True):
        if ftic:
            tic=Tic();

        fp=self.flag_pade;

        xx0=self.xx0;

        for nnn in range(ncount):
            group=self.group;        
            calc_Cx0(self.xc0,self.C0,xx0);
            for g in group:
                group<<g.call('step');
            group.results;
            if fp:
                np.sum(self.xxp,axis=0,out=xx0);
            else:
                xx0+=np.sum(self.xxp,axis=0,out=self.xx00);

        if ftic and ncount:
            self.tStep=tic.sec()/ncount;

        return xx0

    def reset_J(self,J):
        if type(J)==np.ndarray:
            J=({'qp':({'c':(1,)},),'x':J.flatten()},);
        
        if J:
            xj=self.xJs;
            Jold,J=J,[]
            if type(Jold)==dict:
                Jold=(Jold,);
            for k in range(len(Jold)):
                x=Jold[k]['x'];
                xj[k][:]=x;
                q={'qp':Jold[k].get('qp',{}),'x':k};
                J.append(q);
        
        group=self.group;
        for g in group:
            group<<g.call('reset_J',J);
        group.results;
        return self
"""
        if self.ffinished:
            group.results;
        else:
            return group
"""

class solver_test_pp():
    def __init__(self,pade_nm=(2,4),**opts):
        pass

