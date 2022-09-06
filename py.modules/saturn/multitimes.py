#
import os
import numpy as np
from asyn.ipc_matrix import ipc_array,ipc_matrix,schema,ipc_unpack
#from asyn.SharedWorker import SharedWorker,AsynGroup,Tic
import asyn.SharedWorker as sw
from lipa import pade_exp_poles_res
from lipa import kernel
from lipa.kernel import LIPA_solver_st,convert_tocsc,lu_solver_factory
from lipa.parallel import LIPA_solver_pp
from lipa.qpl import *


class LIPA_solver_saturn(kernel.LIPA_forms):
    def __init__(self,AC,SolverFactory=lu_solver_factory(),**opts):
        tic=sw.Tic();
        AC=tolist(AC)
        AC=[convert_tocsc(i) for i in AC];
        N=AC[0].shape[0];
        fcomplex=opts.get('fcomplex',False)
        self.dts=dts=opts['dt']
        self.lipas=[]
        lipas_group=[];
        
        ipcxc0=ipc_array(schema((N,),dtype=float_type))
        ipcxx0=ipc_array(schema((nd,N),dtype=float_type))
        ipcxJs=ipc_array(schema((nJsmax,N),dtype=float_type))
        
        saturn=(ipcxc0,ipcxx0,ipcxJs)
        
        LIPA_solver=opts.get('LIPA_solver',LIPA_solver_pp)

        for dt in dts:
            op=extend({"dt":dt,'ipc_saturn':saturn},opts);
            l=LIPA_solver(AC,SolverFactory,**op);
            g=l.group;
            self.lipas.append(l)
            lipas_group.append(g)
        
        (steps,rJs)=(()())
        for l in self.lipas:
            steps+=(l.step,);
            rJs+=(l.reset_J,)
        
        (self.step,self.reset_J)=(steps,rJs);

        self.xc0=ipcxc0.value;
        self.xx0=xx0=ipcxx0.value;
        self.xJs=ipcxJs.value;

        super(LIPA_solver_saturn,self).__init__(xx0);

        for g in lipas_group:
            g.result;
        self.tinit=tic.sec();

