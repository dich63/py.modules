import os
import sys;sys.path.append('v:/ipc/py.modules')

"""
qin=mp.JoinableQueue()
"""
print(os.getpid(),':',__name__)


def testLU(A,nlu=2,ntic=10):
    
    import asyn.SharedWorker as sw
    import numpy as np

    A=A.tocsc();
    N=A.shape[0];
    b=np.random.rand(N,1);
    tic=sw.Tic();
    for k in range(nlu):
        fac=sla.splu(A);

    tlu=tic.sec()/nlu;
    tic.start();
    for k in range(ntic):
        x=fac.solve(b);

    ti=tic.sec()/ntic;
    return (tlu,ti);

if __name__=='__main__':

    import jsonrpc.jsonclass as jc
    import jsonrpc.sparse_marshal
    import MID.mesh_regions as mr
    import asyn.SharedWorker as sw
    import numpy as np
    import scipy.sparse as sp
    from scipy.sparse import  linalg as sla
    import MID.MIDdecay0 as mmd
    from lipa.parallel import LIPA_solver_st,LIPA_solver_pp
     

    import MID.clib as cl
    import MID.projective2d as prq
    from collections import namedtuple

    tic=sw.Tic();

    o=jc.decode_from_file('P:/KKmeshes/BBBBB/testAC.json',1);

    #Az=o.AC[0].tocsc()+1j*o.AC[1].tocsc();
    Az=o.AC[0]+1j*o.AC[1];
    Ac=Az.tocoo();
    pf=testLU(Az);
    print('Az: perf',pf);
    bA=sp.block_diag((Ac,Ac,Ac,Ac));
    pf=testLU(bA);
    print('bA: perf',pf);
    os.system('pause')
    exit()


    oo=jc.decode_from_file('o:/GT/nominal_centroid_rz22.json')
    jc.compress_mode(1,1);
    oo['vxs']=np.array(oo['vxs'],dtype='float32')
    jc.encode_to_file(oo,'o:/GT/o1f.json')
    oo2=jc.decode_from_file('o:/GT/znominal.json')

    ss='''klklklklklk !
    kkkkkkk!
    '''
    params={
        "tbg":
        {
            "sigma":0
            ,
            "mu":1
            }
        }
 


    print(ss,params)

    r2q=prq.rect2quad
    r=[21,-1500,52,1500]
    r2=[21,-1500,57,1500]
    q,q2=r2q(r),r2q(r2);
    t0=prq.quad2quad_projective(q,q2)
    #os.system('pause')
    qq=cl.vxs2D_projective_deformation([0,0,0,0],[t0],np.array(q))


    tic=sw.Tic()

    norm=np.linalg.norm;
    spnorm=sla.norm;
    fn="V:/MID.js/mid_data_u.json"
    #fn="V:/MID.js/mesh.json"
    #fn="V:/MID.js/three_barrier_far.json"
    #fn="V:/MID.js/outmesh.txt"
    #fn="V:/smesh.txt"
    count=100
    dt=0.001
    tic.start()
    ms=mmd.MIDdecay0(':file:'+fn);
    t=tic.sec()
    print(t,'sec MIDdecay0')

    delta=2
    ms['tubing']['rects'][0][2]+=delta
    ms['tubing2casing']['rects'][0][0]+=delta

    tic.start()
    ms.init(dt=dt,LIPA_solver=LIPA_solver_pp);
    t=tic.sec()
    print(t,'sec :ms.init(LIPA_solver=LIPA_solver_st)')
    tic.start()
    jst=ms.decay(count)
    t=tic.sec()
    print(t,'sec :count=',count)
    tic.start()

    ms.init(dt=dt,LIPA_solver=LIPA_solver_pp);
    t=tic.sec()
    print(t,'sec :ms.init(LIPA_solver=LIPA_solver_pp)')
    tic.start()
    jpp=ms.decay(count)
    t=tic.sec()
    print(t,'sec :count=',count)
    print('end:')
    njst=norm(np.array(jst))
    njpp=norm(np.array(jpp))
    nerr=norm(np.array(jpp)-np.array(jst))
    print('norms:',nerr,njst,njpp)
    print('...')
    print('st:',jst)
    print('pp:',jpp)
    #os.system('pause')
    """
    s=open(fn,'r').read();
    
    mesh=jc.decode(s)
    s=0
    
    tic.start()
    mm=mr.MID_mesh(mesh)
    t=tic.sec()
    print(t,'sec')
    (ss,bms)=mm.sigma_mu_pair();
    (fck,fak)=(mesh['fck'],mesh['fak'])
    print(ss)
    print(fck)
    print(bms)
    print(fak)
    #os.system('pause')
    tic.start()
    trsmask=mm.trs_mask;
    fmask=mesh['fmask']
    nmm=norm(trsmask-fmask)
    fC=ss[trsmask];
    fA=bms[trsmask];
    fc=mesh['fc']
    fa=mesh['fa']
    ncc=norm(fC-fc)
    naa=norm(fA-fa)
    """
    """
    tic.start()
    mm=mr.MID_mesh(mesh)
    t=tic.sec()
    print(t,'MID_mesh sec')
    tic.start()
    sAC=mm.makeAC();
    t=tic.sec()
    (A,C)=[s.tocsc() for s in sAC]
    print(A.shape ,' t=',t,'sec')
    os.system('pause')
    jg=mm.region_mask('sgcoil')
    jd=mm.region_mask('srcoil')
    mA=mesh['matrix']['A'];
    mC=mesh['matrix']['C'];
    nA=spnorm(mA-A)
    nC=spnorm(mC-C)
    print('norms AC ',nA,nC,'...')

    
    print(t,'sec')
    a=11
    """