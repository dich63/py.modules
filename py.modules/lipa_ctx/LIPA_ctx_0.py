#
from parallel.sparse import *
from lipa.kernel import LIPA_forms,pade_exp_poles_res,printf
import parallel.sparse  as psp

from parallel.qp_context import * 
import asyn.SharedWorker as sw

Tic=sw.Tic;
norm=np.linalg.norm
    
    
class polus_invoker_t(object):
    def __init__(self,g,hjet,polus_res,xx0,sources, fcomplex,opts):
         #hinvoke_batch(

        self.opts=opts=copy.copy(opts)

        fzero=opts.fzero;
        fzero=1;   
        polus,res=polus,res=polus_res;
        
        dt=opts.dt;
        parallel=opts.parallel;

        self.zt,self.Bzt=zt,Bzt=polus/dt,res/dt;
        
        self.nd=opts.nd;
        #self.bzt=1.0/zt;
        self.hjet=hjet;
        #(iAz,Az)=hjet.zm(z=zt,scale=Bzt,asyn=0);        
        
        #
        (iAz,Az)=hjet.zm(z=zt,scale=1.0,asyn=parallel);
        
        '''
        Az=hjet.zm(z=zt,scale=1.0,asyn=0);
        
        #Az=hjet.zm(z=zt,scale=1.0,asyn=0);
        n=xx0.shape[1]
        '''
        self.Az=Az;
        
        n=hjet.n;
        opts.n=n;
        
        self.refAz=refAz=Az.ref_context();        
        
        #iAz=hinvoke_context();        
        opts.y=self.y=y=create_buf(n,fzero=fzero);
        
        
        self.sps=sps=sparse_solver_invoker(refAz,to_dict(opts));
        #self.sps=sps=sparse_solver_invoker(Az,to_dict(opts));
        
        ifactor=hinvoke_batch((iAz,sps.factorize),asyn=False,links=(refAz,sps))
        
        #ifactor=hinvoke_context();
        
        if opts.parallel:
            g(ifactor)
        else:
            ifactor()

        self.g=g;         
        self.xx0=xx0;
        self.xxz=create_buf_like(xx0,fzero=fzero);        
        self.jz=create_buf(n,fzero=fzero);
        self.buf0=buf0=create_buf(n,fzero=fzero);
        self.jetbuf=jetbuf=create_buf_like(xx0,fzero=fzero);        
        self.sources=sources;

        #self.y0=y0=create_buf(n,fzero=fzero);

        #self.Cz_invoker=hjet.gaxpy_jet(zt,fnd=xx0,y=y,buf=buf0,freal= not fcomplex).invoker;        
        #self.Cz_invoker()
        
        '''
        ff=np.ones(shape=(5,n),dtype=np.complex)
        #create_buf((5*n),fzero=1,dtype=np.complex);
        self.xx0=xx0=ff.reshape((5,n));
        '''
        self.Cz_invoker=hjet.gaxpy_jet(zt,fnd=xx0,y=y,buf=buf0,freal= not fcomplex).invoker;
               
        
        self.istep=self.__istep_init();
        self.istep=hinvoke_context();
        '''
        
        '''
        self.ifactorize=ifactor;
        
    def __istep_init(self):        
        #self.y[:]=self.jz;

        hjet,y,xx0,xxz,zt,jz=self.hjet,self.y,self.xx0,self.xxz,self.zt,self.jz;
        
        parallel=self.opts.parallel;

        jetbuf=self.jetbuf       
        sources=self.sources

        nd=self.nd;

        #izjz=LA_op_invoker.zero(jz);                
        #isrc_app=source.apply_source(zt,jz);
        
        icpi=LA_op_invoker.copy(y,jz);       
        isrc_app=sources.apply_sources(zt,y);    
        #isrc_app()
        #        isrc_app=hinvoke_context();
        iczi=self.Cz_invoker;        
        #        iczi=hinvoke_context();
        
        
        #ijet=hjet.make_jet(zt,xx0,jetbuf)

        #
        iszi=self.sps.solve; 
        #iszi=hinvoke_context();
        icpib=LA_op_invoker.copy(xxz[0],y);
        
        #iop=[izjz,isrc_app,icpi,iczi,iszi,icpib];
        iop=[icpi,isrc_app,iczi,iszi,icpib];
        #iop=[icpi,isrc_app,iszi,icpib];
        
        if nd>1:
            ijet=hjet.make_jet(zt,ff=xx0,jet=xxz,offset=0,asyn=parallel)        
            #
            iop.append(ijet);
            #jj=hinvoke_batch(iop,asyn=False)       
            #jj()
            #ijet()

        isrc_tstep=sources.transform_step();
        iop.append(isrc_tstep);

        i=hinvoke_batch(iop,asyn=False)        
        #i=hinvoke_batch((icpi,iczi,iszi,icpib),asyn=False)        
        #i()
        return i;
        
    
     
        

def LIPA_solver_ctx_opts_default(opts):
    o=jsobject(ext_def(opts,{
        
        'fhuge':False,
        'common':{},        
        'solver_factory':klu_context_factory,
        'pade_nm':(2,4),
        'nd':1,
        'dt':1.0,
        'fcomplex':False ,
        'asyn':False ,              
        'fzero':True,
        'xx_init':[], 
        'parallel':True,        
        }))    
    return o;
    
        
class LIPA_solver_ctx(LIPA_forms):
    
    def __init__(self,jet_or_AC,**opts):
        
        self.opts=opts=LIPA_solver_ctx_opts_default(opts);
                
        tic=Tic();
        tic.start();

        self.opts=opts;
        self.zz=opts;
        dt=opts.dt
        parallel=opts.parallel

        self.fcomplex=fcomplex=opts.fcomplex       
        self.fhuge=fhuge=opts.fhuge        
        self.pade_nm=pade_nm=opts.pade_nm   
        
         
        self.hjet=hjet=hjet_context_create(jet_or_AC,fhuge,parallel=parallel);
        
        self.fhuge=fhuge=hjet.fhuge;
        
        opts.float_type=float_type= np.complex128 if fcomplex else np.float64;
        
        self.count=count=hjet.count;
        
        nd=opts.nd=max(opts.nd,hjet.count-1);
        
        
        
        
        #self.polus_res_list=polus_res_list=pade_exp_poles_res.get(pade_nm,not fcomplex)
        
        opts.N=N=hjet.n;
        
        xx=opts.xx
        if xx is None:
            self.xx=xx=create_buf((nd,N),fzero=opts.fzero,dtype=float_type);
        

        for k in range(min(len(opts.xx_init),count)):
            xx[k][:]=opts.xx_init[k];
            
        self._j=create_buf(N);

        self.sources=sources=qp_context(N,dt);
        
        LIPA_forms.__init__(self,xx);
        self.g=g=pp_group(sync= not parallel);
        
        self.polus_res_s=prs=pade_exp_poles_res.get(pade_nm,not fcomplex)
        
        
             
        self.polus_solvers=pss=[  polus_invoker_t(g,hjet,pr,xx,sources,fcomplex,opts)    for pr in prs]
        
        
        def istep_init(parallel=1,grain=0):           
            
            i=[ps.istep for ps in pss];            
            ipoles=hinvoke_batch(i,asyn=parallel,links=i);         
                
            
            #xxzs=[ptr_array(ps.xxz) for ps in pss]
            xxzs=[ps.xxz for ps in pss]
            Bzts=[-ps.Bzt for ps in pss]
            
            

            
            
            #
            ipolessum=LA_op_invoker.linspace(xxzs,Bzts,xx,grain=grain,parallel=parallel);
            #            ipolessum=LA_op_invoker.sum(xx=xxzs,y=xx,grain=grain,parallel=1);
            #ipolessum=LA_op_invoker.copy(xx[0],xxzs[0][0]);
            #ipoles()
            #ipolessum();

            nm=opts.pade_nm;
            
            if  nm[0]<nm[1]:
                iz=LA_op_invoker.zero(xx,grain=grain,parallel=parallel);
                iop=(ipoles,iz,ipolessum)
            else:
                iop=(ipoles,ipolessum)
            
            istep=hinvoke_batch(iop,asyn=False,links=iop);                     
            
            return istep;
        
        
        
        #self.ii=istep_init
        self.istep=istep_init(parallel=opts.parallel)#,asyn=opts.asyn)
        

        #
        
        self.iwait=iwait=g.iwait;    
        self.init=iwait

        #self.istep=istep=istep_init();    
        #self.step=istep
        '''
        if opts.parallel:
            [g(p.ifactorize) for p in pss]
        else:
            [p.ifactorize() for p in pss]
        '''
        #g.join()
                

        if not opts.asyn:
            iwait();
            self.tinit=tic.sec();
            self.tloading=0;  
        

                        
    def step(self,nrep=1):   
        tic=Tic();    
        for n in range(nrep):
            self.istep()
        self.tStep=tic.sec()/nrep;
        return self;
        
    def step_t(self,nrep=1):
         self.step()
         t=self.tStep
         printf('\r tic=%f sec',t)
         return t;


    def reset_corr(self,n=1):
        pass

    def iset_sources(self,j):
        n,dt,sources=self.hjet.n,self.opts.dt,self.sources
        return sources.set_sources(j,dt,n);                        
        
    def reset_J(self,j):
        iss=self.iset_sources(j);
        self.__current=iss; 
        '''
        ia=self.sources.apply_sources(1.+1.j,y)
        iss()
        ia()
        '''
        err=iss()
        '''        
        y=self.sources.get_current();
        print('jcxt=',norm(y));
        '''
        return err;
            
        
class LIPA_solver_group_ctx(LIPA_forms):
    
    def __init__(self,jet_or_AC,**opts):

        tic=Tic();
        tic.start();
        self.tic=tic;
        self.opts=opts=LIPA_solver_ctx_opts_default(opts);
        asyn_group=opts.asyn

        dts=tolist(opts.dt);
               
        self.dts=dts;
        opts.asyn=asyn=opts.parralel;
        self.g=g=pp_group(sync=not asyn);
        self.init=init=g.iwait;


        self.hjet=hjet=hjet_context_create(jet_or_AC,fhuge=opts.fhuge,parallel=0);

        fcomplex=opts.fcomplex;
        self.fhuge=fhuge=hjet.fhuge;        
        opts.float_type=float_type= np.complex128 if fcomplex else np.float64;        
        self.count=count=hjet.count;
        

        nd=opts.nd=max(opts.nd,hjet.count-1);

        opts.N=N=hjet.n;
        
        xx=opts.xx
        if xx is None:
            self.xx=xx=create_buf((nd,N),fzero=opts.fzero,dtype=float_type);
        

        for k in range(min(len(opts.xx_init),count)):
            xx[k][:]=opts.xx_init[k];

        LIPA_forms.__init__(self,xx);

        opts.xx_init=[];
        opts.xx=xx;

        optsdt=copy.copy(opts)
        solvers=[];         
        for dt in dts:
            optsdt.dt=dt;
            s=LIPA_solver_ctx(hjet,**to_dict(optsdt));
            g(s.init); 
            solvers.append(s);                       
                         
        self.solvers=solvers;
        self.steps=[s.step for s in solvers];
        if not asyn_group:
            init();

    def stepn(self,ndt=0,nrep=1):
        solver=self.solvers[ndt];    
        tic=self.tic
        tic.start();
        for n in range(nrep):
            solver.istep()
        self.tStep=tic.sec()/nrep;
        return self;

    def step_t(self,nrep=1):
         self.step()
         t=self.tStep
         printf('\r tic=%f sec                ',t)
         return t;
         
    def step(self,nrep=1):
        return self.stepn(ndt=0,nrep=nrep);


    def reset_J(self,j):
        solvers=self.solvers;
        err=0
        for s in solvers:
            err=s.reset_J(j);            
            if err:
                return err;            
        
        return 0;


        
        
        #solver_factory=cast_str(opts.get('solver_factory',klu_context_factory))
        
#l=LIPA_solver_ctx([[1],[2]])            
    
def test_batch(AC=([1],[2],[3])):
    
    hjet=hjet_context_create(AC,parallel=True);
    g=pp_group();
    
    (invAz,Az)=hjet.zm(2,1j,asyn=True);
    (invAz2,Az2)=hjet.zm(2,10j,asyn=True);
    o=jsobject();
    o.n=hjet.n;
    pp=Az.ref_context();    
    d=o.__dict__
    sps=sparse_solver_invoker(pp,o.__dict__);
    
    hib=hinvoke_batch((invAz,invAz2,sps.factorize),asyn=False,links=(pp,sps))
    return (hib,sps,Az,Az2)

#test_batch()  
#ss=sparse_solver_invoker([2])      