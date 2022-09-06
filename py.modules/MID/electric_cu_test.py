import sys,os;sys.path.append('v:/ipc/py.modules')

import MID.MIDdecay_fs_ec as mid

import jsonrpc.jsonclass as jsc
import numpy as np
from utils.derr2m import derr2m

import MID.e_circuit as ec

import asyn.SharedWorker as sw

#os.environ['klu.common.scale']='1'


ec.test();

tic=sw.Tic()




Ms=1e-3;
times=[(0.02*Ms,10*Ms ),(0.5*Ms,80*Ms)];
times=[(0.02*Ms,5*Ms ),(1*Ms,200*Ms)];
#
times=[(0.02*Ms,1*Ms)];
EMF=2000;

mesh=mid.reparse_mesh(":file:V:/ipc/MID/mesh_nominal_LS_SSz.json")
exdata=jsc.decode_from_file('V:/ipc/MID/tool_decays_stainless.json')
barriers=exdata['barriers'];
d={"sigma":0.0,"mu":1.0};

mid.region_by_name_update(mesh,{"tube0":d,"tube1":d,"tube2":d})

mesh['regions'][6]['disabled']=1;
mesh['regions'][7]['disabled']=1;
mesh['regions'][8]['disabled']=1;
mesh['regions'][9]['disabled']=1;
sol=mid.MID3barrier_decay(mesh=mesh,rcoil="m_rcoil",gcoil="m_gcoil")

# old scheme...

sol.init_ms(barriers=barriers,times=times,pade_nm=[2,4]);


#cn=sol.solver.solvers[0].polus_solvers[0].sps.cond
#print('estimate condnum={:e}'.format(cn));
#sol=mid.MID3barrier_decay(mesh=mesh,rcoil="m_rcoil",gcoil="m_gcoil")
#common={"numerical_rank":-4,"scale":1}
#sol.init_ms(barriers=barriers,times=times,pade_nm=[2,4],common=common);

cn=sol.solver.solvers[0].polus_solvers[0].sps.cond

print('estimate condnum={:e}'.format(cn));

tic.start()
decs_old=sol.decay_ms(EMF,0);
print('tic=',tic.sec(),'s')

'''
decs_oldH=sol.decay_ms2(EMF,0);


h=decs_oldH['j']
j=decs_old['j']

mj=np.mean(j)
mh=np.mean(h)
J=j/mj;
H=j/mh;
drr=derr2m(J,H);
print('derr',drr*100,'%')
'''



e_circuit=[
    {'id':"s_gcoil",
     'R':1108.1,
     'coils':[{'Nw':770}]     
    },
    {'id':"m_gcoil",
     'R':301.9,
     'coils':[{'Nw':594}]     
    },
    {'id':"s_rcoil",
     'R':1e6,
     'coils':[{'Nw':1000}]     
    },
    {'id':"m_rcoil",
     'R':1e6,
     'coils':[{'Nw':1000}]     
    }
    ];

'''
e_circuit=[    
    {'id':"m_gcoil",
     'R':301.9,
     'coils':[{'Nw':594}]     
    }
    ];

'''

print('e_circuit:',e_circuit)

Rg=e_circuit[0]['R'];
U=EMF*Rg;


sol=mid.MID3barrier_decay(mesh=mesh,rcoil="s_rcoil",gcoil="s_gcoil")

# new scheme...

jsc.compress_mode(1)

common={"numerical_rank":-4,"scale":1}
common={"numerical_rank":4,"scale":1}
common={}

sol.init_ms(barriers=barriers,e_circuit=e_circuit,times=times,pade_nm=[2,4],common=common);


cn=sol.solver.solvers[0].polus_solvers[0].sps.cond
print('estimate condnum={:e}'.format(cn));


tic.start()
decs_new=sol.decay_ms_ec2({'U':U,'id':'m_gcoil'},0);
print('tic=',tic.sec(),'s')


_r={'A':sol.A.tocsc(),'C':sol.C.tocsc()}
jsc.encode_to_file(_r,'c:/temp/fro_test_ec3.json');

#sys.exit(0)


#r={'old':decs_old,'new':decs_new,'experiment':exdata['SHORT']};

#old={'t',decs_ols['']

to=decs_old['t'];
jo=decs_old['j'];



ix=decs_new['ids']['m_rcoil'];
jj=decs_new['jj'];
jn=jj[ix,0,:];
tn=decs_new['t'];

te=exdata['SHORT']['time']
je=exdata['SHORT']['decay']

_r={'old':{'t':to,'j':jo},'new':{'t':tn,'j':jn},'experiment':{'t':te,'j':je},'econd':cn};


s=jsc.encode_to_file(_r,'v:/ipc/mid/ec_testr1.json');
print(s)