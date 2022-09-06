import sys,os;sys.path.append('v:/ipc/py.modules')

import MID.MIDdecay_fs_ec as mid

import jsonrpc.sparse_marshal
import jsonrpc.jsonclass as jsc
import numpy as np
from utils.derr2m import derr2m

import MID.e_circuit as ec

import asyn.SharedWorker as sw
import MID.dc as dcc

#os.environ['klu.common.scale']='1'



ec.test();

tic=sw.Tic()


qq=dcc.dc_profile(.4,17.0,1)

Ms=1e-3;
times=[(0.02*Ms,10*Ms ),(0.5*Ms,80*Ms)];
times=[(0.02*Ms,5*Ms ),(1*Ms,200*Ms)];
#times=[(0.02*Ms,1*Ms)];
EMF=2000;

barriers={
"pipe1":{
                   "d":50.0,
                   "th":11.9891,
                   "th0":0.0,
                   "sigma":5e6,
                   "mu":80.0,
                   "z2r":1.5,
                   "dcr":1,
                   "dc":0#2                  

                }
}


e_circuit=[
    {'id':"gcoil",
     'R':1108.1,
     'coils':[{'Nw':770}]     
    },
    {'id':"rcoil",
      'R':1e6,
     'coils':[{'Nw':1000}]     
    }
    ];


fn=":file:V:/work/MID3D/MEDIUM_1_decentralized.json"
fn=":file:V:/work/MID3D/1111.json"
#fn=":file:V:/work/MID3D/mesh_nominal_SSz.json"
fn=":file:V:/work/MID3D/MEDIUM_2.json"
#fn=":file:V:/work/MID3D/MEDIUM_2_2_d.json"

fn=":file:V:/work/MID3D/emp_8S_1.json"
fn=":file:V:/work/MID3D/mesh_nominal_SSz_morph.json"
#
fn=":file:V:/work/MID3D/MEDIUM_1_0_morph.json"
fn=":file:V:/work/MID3D/MEDIUM_2.json"
mesh=mid.reparse_mesh(fn)
e_circuit=mesh['e_circuit'];

rpipe1=mid.region_by_name(mesh,'pipe1');

jsc.compress_mode(1)

#mid.region_by_name_update(mesh,{"pipe1":{"sigma":0.0,"mu":1000.0,"mu_r":1000.0}})

mesh['sigma_mu_bg']=(1.11,1.0);
sol=mid.MID3barrier_decay(mesh=mesh,rcoil="rcoil",gcoil="gcoil")

#mid.region_by_name_update(mesh,{"pipe1":{"z2r":2}})
#



barriers={
"pipe1":{
                   "d":73.025,
                   "th":5.512,
                   "th0":5.512,
                   "sigma":3789413.815,
                   "mu":109.616             

                },
"pipe2":{
                   "d":244.47,
                   "th":10.033,
                   "th0":10.033,
                   "sigma":9966615.039,
                   "mu":78.376             

                }


}

barriers={
"pipe1":{
                   "dc":0.1

                },
"pipe2":{
                   "d":244.47,
                   "th":10.033,
                   "th0":10.033,
                   "sigma":9966615.039,
                   "mu":78.376             

                }


}


sol.init_ms(fmbc=1,barriers=barriers,e_circuit=e_circuit,times=times);
newmesh=sol.mesh_reconstruct(barriers=barriers,clsdis=0);
#jsc.encode_to_file(newmesh,"V:/work/MID3D/MEDIUM_2_0_morph.json");
#jsc.encode_to_file(newmesh,"V:/work/MID3D/mesh_nominal_SSz_morph.json");

bzr=sol.getBzBy()



#jsc.encode_to_file({'A':sol.A.tocoo(),'C':sol.C.tocoo()},"V:/work/MID3D/AC.json");
#jsc.encode_to_file(newmesh,"V:/work/MID3D/MEDIUM_1_decentralized_morph.json");


#exit(1)
mesh=mid.reparse_mesh(":file:V:/ipc/MID/mesh_nominal_LS_SSz.json")
exdata=jsc.decode_from_file('V:/ipc/MID/tool_decays_stainless.json')
barriers=exdata['barriers'];
d={"sigma":0.0,"mu":1000.0};
da=d;
da={"sigma":0.0,"mu":1000.0,"mu_r":1000.0};
mid.region_by_name_update(mesh,{"tube0":d,"tube1":d,"tube2":d})

'''
mesh['regions'][6]['disabled']=1;
mesh['regions'][7]['disabled']=1;
mesh['regions'][8]['disabled']=1;
mesh['regions'][9]['disabled']=1;
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

sol.init_ms(fmbc=1,barriers=barriers,e_circuit=e_circuit,times=times,pade_nm=[2,4],common=common);


cn=sol.solver.solvers[0].polus_solvers[0].sps.cond
print('estimate condnum={:e}'.format(cn));


tic.start()
decs_new=sol.decay_ms_ec2({'U':U,'id':'m_gcoil'},0);
print('tic=',tic.sec(),'s')
jsc.encode_to_file(decs_new,'c:/temp/jj.json')

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