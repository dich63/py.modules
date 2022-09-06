

from ctypes import *
import numpy as np

def clipj(j):
    j=np.clip(-np.array(j),0,np.Inf)
    j[j==0]=1e-5;
    return j


#import hack_mkl


import os,sys;sys.path.append(os.getenv('ltx_js_root')+'py.modules')







from klu_lipa import *

import matplotlib.pyplot as plt
from matplotlib.pyplot import cm,bar,legend
from lipa.solvers import *
from scipy.sparse import csc_matrix,coo_matrix
from parallel.sparse import  sparse_solver_invoker
import MID.MIDdecay_fs_ec as mmd

def ctx_test():
    dt=1
    A=[[0,0,0],
       [0,1,1],
       [0,0,1]]

    C=[[1,0,0],
       [0,0,0],
       [0,0,0]]

    b=11*np.array([0,0,1.0])

    sA=coo_matrix(A,dtype=np.float64);
    sC=coo_matrix(C,dtype=np.float64);
    sAC=sC-sA;

    si=sparse_solver_invoker(sAC);
    err=si.factorize();
    err=si.factorize.status;
    si.y=b;
    y=si.y;
    err=si.solve()
    err=si.solve()
    err=si.solve()

    sol=LIPA_solver_ctx((sA,sC),dt=dt,fcomplex=False,pade_nm=(0,1),nd=2);
    err=sol.status;
    sol.reset_J(b);
    xx=sol.xx[0];
    sol.step();
    sol.step();
    return 0

def ec_test0(mesh):
    
    barriers={}
    jcc=[
        {'R':1e0,'coils':[
            {'Nw':1000,'name':"gcoil"}]},
            {'R':55,'coils':[
                {'Nw':100,'name':"rcoil1"},
                {'Nw':-100,'name':"rcoil2"}
                ]
             }
        ]

    count=100;
    Ms=1e-3;
    tb=0.314*Ms;
    T=15.5*Ms;
    dt=T/count;
    td=0.3*Ms;

    mdec=mmd.MID3barrier_decay(mesh,rcoil='srcoil',gcoil='gcoil');
    
    #    jcc={};
    #mdec.init_ms(barriers=barriers,times=[[dt,T]],e_circuit=jcc,pade_nm=(0,1));
    mdec.init_ms(barriers=barriers,times=[[dt,T]],e_circuit=jcc);
    s=mdec.solver.solvers[0]
    cn=s.polus_solvers[0].sps.cond
    print np.log10(cn)
    #mdec.init(barriers=barriers,dt=dt,e_circuit=jcc,pade_nm=(0,1));
    st=hex(mdec.solver.status)
    isp=mdec.solver.solvers[0].istep;
    err=isp()
    #ss01=mdec.decay_ms()
    A,C=(mdec.A,mdec.C)
    sA=A.tocoo()
    sC=C.tocoo()
    sAC=sC/dt-sA;

    si=sparse_solver_invoker(sAC);
    err=si.factorize();
    err=si.factorize.status;
    cn=si.cond;

    sol=LIPA_solver_ctx((sA,sC),dt=dt,fcomplex=False,pade_nm=(0,1),nd=2);
    err=sol.status;


    ss01=mdec.decay_ms_ec(737.111);
    ss00=mdec.decay_ms_ec(0,None);
    return 0;

def ec_test(mesh):
    
    barriers={}
    jcc=[
        {'R':10,'coils':[
            {'Nw':1000,'name':"gcoil"}]},
            {'R':1000,'coils':[
                {'Nw':-100,'name':"rcoil1"},
                {'Nw':100,'name':"rcoil2"}
                ]
             }
        ]

    count=500;
    Ms=1e-3;
    ms=10*1e-6;
    T=10*Ms/1;
    #T=10*ms/1;
    dt=T/count;
    
    r=mesh['regions'][6]
    r['sigma']=r['sigma']*1e-1;
    r['mu']=1;
    mdec=mmd.MID3barrier_decay(mesh,rcoil='srcoil',gcoil='gcoil');
    
    #    jcc={};
    #mdec.init_ms(barriers=barriers,times=[[dt,T]],e_circuit=jcc,pade_nm=(0,1));
    mdec.init_ms(barriers=barriers,times=[[dt,T]],e_circuit=jcc,pade_nm=(2,4));
    s=mdec.solver.solvers[0]
    cn=s.polus_solvers[0].sps.cond
    print np.log10(cn)    
    st=hex(mdec.solver.status) 
    


    ss01=mdec.decay_ms_ec(737.01);
    ss00=mdec.decay_ms_ec(0,None);
    return 0;

########
nominal_mesh=':file::y:/mid/3ESHORT_3_old.json'
#nominal_mesh=':file::y:/mid/3ESHORT_3_new.json'
#nominal_mesh=':file::t:/3ESHORT_3_old_z.json'
#nominal_mesh=':file::V:\ipc\MID\mesh_nominal_2z.json'
nominal_mesh=':file::V:\ipc\MID\8s.json'

nominal_mesh=mmd.reparse_mesh(nominal_mesh)
#ctx_test();
ec_test(nominal_mesh);

barriers={
            "pipe1":{
                   "d":46.0,
                   "th":1.0,
                   "th0":1.0,
                   "sigma":0,
                   "mu":1.0
                },
            "pipe2":{
                   "d":56.,
                   "th":1.0,
                   "th0":1.0,
                   "sigma":0,
                   "mu":1
                },
            "pipe3":{
                   "d":300.,
                   "th":18.,
                   "th0":18.0,
                   "sigma":5e6,
                   "mu":40.0
                },
            "sensor":'LS',
            'zmax':2500
            }

#barriers=None;

count=300;
Ms=1e-3;
tb=0.314*Ms;
T=15.5*Ms;
dt=T/count;
td=0.3*Ms;


mdec=mmd.MID3barrier_decay(nominal_mesh,gcoil='sgcoil',rcoil='srcoil');

mdec.init_ms(barriers=barriers,times=[[dt,T]]);

import jsonrpc.jsonclass as jsncls
#m=mdec.mesh.mesh
m=mdec.new_mesh

jsncls.compress_mode(1)
jsncls.encode_to_file(m,'t:/001.json')

b=np.ones(nominal_mesh['vxs'].shape[0],dtype='d');

J=({'Q':[1],'b':b},{'Q':[-1],'z':-1/td,'b':b})

J=({'Q':[1]},{'Q':[-1],'z':-1/td})
J=({'Q':[1]},)
#mdec.solver.reset_J(J);
#mdec.init_ms(barriers=barriers,times=[[dt1,T1],[dt0,T0]]);

ss0=mdec.decay_ms();

#ss=mdec.decay_ms(J);
'''
t,jn=ss0['t']/Ms,ss0['j']
plt.plot(t,clipj(jn), color='k',label=r'$\tau= 0$ Ms')


plt.plot(t,jn, color='k',label=r'$\tau= 0$ Ms')
t,jn=ss['t']/Ms,ss['j']
plt.plot(t,jn, color='r',label=r'$\tau=J$ Ms')


plt.legend()
plt.title(r' J(t)~ $1-e^{-\frac{t}{\tau}}$')
print('save...')
plt.savefig('zzz.svg')
plt.savefig('zzz.png')
print('ok')
os.system('start zzz.svg')
'''

'''
td=0.01*Ms;
J=({'Q':[1]},{'Q':[-1],'z':-1/td})
ss001=mdec.decay_ms(J);
'''

td=0.05*Ms;
J=({'Q':[1]},{'Q':[-1],'z':-1/td})
ss01=mdec.decay_ms(J);





td=0.1*Ms;
J=({'Q':[1]},{'Q':[-1],'z':-1/td})
ss3=mdec.decay_ms(J);


td=0.5*Ms;
J=({'Q':[1]},{'Q':[-1],'z':-1/td})
ss05=mdec.decay_ms(J);

td=1*Ms;
J=({'Q':[1]},{'Q':[-1],'z':-1/td})
ss10=mdec.decay_ms(J);

td=15*Ms;
J=({'Q':[1]},{'Q':[-1],'z':-1/td})
ss50=mdec.decay_ms(J);



#from matplotlib.pyplot import cm,bar,legend
t,jn=ss0['t']/Ms,ss0['j']
f=t>(tb/Ms)
f=t>0
t,jn=t[f],jn[f];
plt.plot(t,clipj(jn), color='k',label=r'$\tau= 0$ Ms')

t,jn=ss01['t']/Ms,ss01['j']
t,jn=t[f],jn[f];
plt.plot(t,clipj(jn), color='r',label=r'$\tau=0.05$Ms')

'''
t,jn=ss001['t']/Ms,ss01['j']
t,jn=t[f],jn[f];
plt.plot(t,clipj(jn), color='m',label=r'$\tau=0.01$Ms')
'''




t,jn=ss3['t']/Ms,ss3['j']
t,jn=t[f],jn[f];
plt.plot(t,clipj(jn), color='b',label=r'$\tau=0.1$Ms')


t,jn=ss05['t']/Ms,ss05['j']
t,jn=t[f],jn[f];
plt.plot(t,clipj(jn), color='g',label=r'$\tau=0.5$Ms')

t,jn=ss10['t']/Ms,ss10['j']
t,jn=t[f],jn[f];
plt.plot(t,clipj(jn), color='y',label=r'$\tau=1.0$Ms')

t,jn=ss50['t']/Ms,ss50['j']
t,jn=t[f],jn[f];
plt.plot(t,clipj(jn), color='c',label=r'$\tau=15.0$Ms')


plt.xscale('log')
plt.yscale('log')
plt.legend()
plt.title(r' J(t)~ $1-e^{-\frac{t}{\tau}}$')
print('save...')
plt.savefig('zzz.svg')
plt.savefig('zzz.png')
print('ok')
os.system('start zzz.svg')




'''  
    plt.plot(t,clipj(j), color='b',label='genmesh')
    
    plt.plot(t,clipj(jd), color='g',label='deform',marker='.',ls='',markersize=4)
'''    
   

'''
windll.kernel32.SetDefaultDllDirectories(0x00000000);





import os,sys;sys.path.append(os.getenv('ltx_js_root')+'py.modules')

import jsonrpc.jsonclass as jc
import jsonrpc.sparse_marshal





from klu_lipa import *

o=jc.decode_from_file('t:/test.json')
print(o)
solver=LIPA_solver_ctx(**o)


j=o['Jo']
solver.reset_J(j)
print(solver)
'''