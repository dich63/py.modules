#
from parallel.sparse import *
from lipa.kernel import LIPA_forms,pade_exp_poles_res,printf
import parallel.sparse  as psp


from parallel.qp_context import * 
import asyn.SharedWorker as sw

Tic=sw.Tic;
norm=np.linalg.norm
    
    
class z_invoker_t(object):
    def __init__(self,g,hjet,z,opts):
         #hinvoke_batch(

        opts=copy.copy(opts)

        fzero=opts.fzero;
        fzero=1;                 
        
        self.hjet=hjet;        
        
        (iAz,Az)=hjet.zm(z=z,scale=1.0,asyn=True);
        #print(z,iAz,Az)
        self.z=z;
        self.Az=Az;
        
        n=hjet.n;
        opts.n=n;
        
        self.refAz=refAz=Az.ref_context();        
        
        #iAz=hinvoke_context();        
        opts.y=self.y=y=create_buf(n,fzero=fzero);
        
        
        self.sps=sps=sparse_solver_invoker(refAz,to_dict(opts));

        #iAz();
        #        ifactor=sps.factorize
        #        ifactor=hinvoke_batch((sps.factorize,),asyn=False,links=(refAz,sps))
        #iAz=hinvoke_context();        
        #
        ifactor=hinvoke_batch((iAz,sps.factorize),asyn=False,links=(refAz,sps))
        
        #ifactor=hinvoke_context();
        
        if opts.parallel:
            g(ifactor)
        else:
            ifactor()

        self.g=g;                        
        self.solve=sps.solve;
        #self.__ifactor=ifactor
        
        
        
    
     
        

def LIPA_solver_ctx_opts_default(opts):
    o=jsobject(ext_def(opts,{
        
        'fhuge':False,
        'common':{},        
        'solver_factory':klu_context_factory,                 
        'asyn':False ,              
        'fzero':True,
        'parallel':True,        
        }))    
    return o;
    
def to_list(l):
    
    if type(l) in (tuple,list,np.ndarray):        
        return l;
    else:
        return (l,);

    
        
class LIPA_solver_ctx_z(object):
    
    def __init__(self,jet_or_AC,**opts):
        
        self.opts=opts=LIPA_solver_ctx_opts_default(opts);
                
        tic=Tic();
        tic.start();

        self.opts=opts;
        
        zs=to_list(opts.z);

        self.zs=zs=np.array(zs,dtype=np.complex128).flatten();

        self.count_z=len(zs);
        
        self.fhuge=fhuge=opts.fhuge                
        
        #print('jet_or_AC:',jet_or_AC)
        self.hjet=hjet=hjet_context_create(jet_or_AC,fhuge,parallel=opts.parallel);
        #print('hjet type:',type(hjet))
        self.N=hjet.n;
        
        self.fhuge=fhuge=hjet.fhuge;        
        
        self.count=count=hjet.count;        
        

        #self.g=g=pp_group_ex(opts.parallel);        
        self.g=g=pp_group_ex(1 if opts.parallel else 0);        
              
        self.z_invokers=z_invokers=[z_invoker_t(g,hjet,z,opts) for z in zs];             
            
        #
        self.iwait=iwait=g.iwait;    
        self.init=iwait        
        if not opts.asyn:
            iwait();
        

                        
    def make(self,J,xx=None):   
        c=self.count_z;
        N=self.N;
        g=self.g;
        z_invokers=self.z_invokers
        if xx is None:
            xx=np.zeros((c,N),dtype=np.complex128);
        for i in z_invokers:
            i.y[:]=J;
            #i.solve();
            g(i.solve);
        
        g.join();
            
        for k in range(c):
            xx[k]=z_invokers[k].y;
        
        return xx;    
        
        
        
