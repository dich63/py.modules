import sys;sys.path.append('v:/ipc/py.modules')

import MID.MIDdecay_fs_ec as mid

import jsonrpc.jsonclass as jsc
import numpy as np
from utils.derr2m import derr2m

import asyn.SharedWorker as sw

tic=sw.Tic()


Ms=1e-3;
times=[(0.02*Ms,10*Ms ),(0.5*Ms,80*Ms)];
times=[(0.02*Ms,1*Ms)];
EMF=2000;

mesh=mid.reparse_mesh(":file:V:/ipc/MID/mesh_nominal_LS_SSz.json")
exdata=jsc.decode_from_file('V:/ipc/MID/tool_decays_stainless.json')
barriers=exdata['barriers'];
d={"sigma":0.0,"mu":1.0};

mid.region_by_name_update(mesh,{"tube0":d,"tube1":d,"tube2":d})



sol=mid.MID3barrier_decay(mesh=mesh,rcoil="s_rcoil",gcoil="s_gcoil")

# old scheme...

sol.init_ms(barriers=barriers,times=times,pade_nm=[2,4]);
tic.start()
decs_old=sol.decay_ms(EMF,0);
print('tic=',tic.sec(),'s')


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


e_circuit_mix=[    
    {'id':"s_gcoil",
     'R':1108.1,
     'coils':[{'Nw':770}]     
    },
    {'id':"m_gcoil",
     'R':301.9,
     'coils':[{'Nw':594}]     
    }
    ];


print('e_circuit:',e_circuit)

Rg=e_circuit[0]['R'];
U=EMF*Rg;

sol=mid.MID3barrier_decay(mesh=mesh,rcoil="s_rcoil",gcoil="s_gcoil")

# new scheme...

sol.init_ms(barriers=barriers,e_circuit=e_circuit,times=times,pade_nm=[2,4]);

tic.start()
decs_new=sol.decay_ms_ec2({'U':U,'id':'s_gcoil'},0);
print('tic=',tic.sec(),'s')
sol=mid.MID3barrier_decay(mesh=mesh,rcoil="s_rcoil",gcoil="s_gcoil")

# mix scheme...

sol.init_ms(barriers=barriers,e_circuit=e_circuit_mix,times=times,pade_nm=[2,4]);

tic.start()
decs_mix=sol.decay_ms_ec_mix({'U':U,'id':'s_gcoil'},0);
print('tic=',tic.sec(),'s')
jsc.compress_mode(1)

#r={'old':decs_old,'new':decs_new,'experiment':exdata['SHORT']};

#old={'t',decs_ols['']

to=decs_old['t'];
jo=decs_old['j'];

tmix=decs_mix['t'];
jmix=decs_mix['j'];

ix=decs_new['ids']['s_gcoil'];
jj=decs_new['jj'];
jn=jj[ix,0,:];
tn=decs_new['t'];

te=exdata['SHORT']['time']
je=exdata['SHORT']['decay']

_r={'old':{'t':to,'j':jo},'mix':{'tmix':to,'jmix':jo},'new':{'t':tn,'j':jn},'experiment':{'t':te,'j':je}};


s=jsc.encode_to_file(_r,'v:/ipc/mid/ec_test.json');
print(s)
