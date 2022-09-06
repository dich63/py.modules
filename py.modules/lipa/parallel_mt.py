import os
import numpy as np
from asyn.ipc_matrix import ipc_array,ipc_matrix,schema,ipc_unpack
#from asyn.SharedWorker import SharedWorker,AsynGroup,Tic
import asyn.SharedWorker as sw
from . import pade_exp_poles_res
from . import kernel
from .kernel import LIPA_solver_st,convert_tocsc,lu_solver_factory

from .qpl import *

from .parallel import LIPA_poles_thread_pp,LIPA_solver_pp

class LIPA_solver_pp_mt(kernel.LIPA_forms):
    def __init__(self,AC,SolverFactory=lu_solver_factory(),**opts):
        tic=sw.Tic();
        AC=tolist(AC)
        #AC=[convert_tocsc(i) for i in AC];
        pade_nm=opts.get('pade_nm',(2,4))
        fcomplex=opts.get('fcomplex',False)
        fsmd=opts.get('fsmd',False)
        cpus=opts.get('cpus',sw.affinity2cpus())
        dts=opts.get('dts',np.double(1))
        nts=len(dts);
        nd=opts.get('nd',1)
        self.nJsmax=nJsmax=opts.get('nJs',1);
        nd=max(len(AC)-1,1,nd);
        N=AC[0].shape[0];

#    (self,pade_nm,AC,dt=1,nd=1,affinity=0,fcomplex=False,SolverFactory=kernel.lu_solver_factory(),**opts):

        float_type= np.complex128 if fcomplex else np.float64;
        
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
        self.group=group=sw.AsynGroup();
        for n in range(Np):
             pt = packets[n]
             group<<sw.SharedWorker(LIPA_poles_thread_pp,affinity_mask=pt[1],job=job)\
             (pt[0],dt,ipcAC,ipcxx0,ipcxc0,nd,ipcxxp,n,ipcxJs,SolverFactory)#\
             #.register_methods(('step','reset_J'));
        if ffinished:
            group.results
        self.tinit=tic.sec();

    def reset_corr(self,n=0):
        group=self.group;

        for g in group:
            group<<g.call('reset_corr',n);
        group.results;


    def step(self):
        group=self.group;
        xx0=self.xx0;
        calc_Cx0(self.xc0,self.C0,xx0);
        for g in group:
            group<<g.call('step');
        group.results;
        return np.sum(self.xxp,axis=0,out=xx0);

    def reset_J(self,J):
        if type(J)==np.ndarray:
            J=({'qp':({'c':(1,)},),'x':J},);
        
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