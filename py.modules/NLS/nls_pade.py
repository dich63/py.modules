#from parallel.parallel_wrapper import *
from ctypes import *
import os
import numpy as np
import NLS.nls_env as lib_env
#import context_wrapper.context_wrapper as ctx

from p23 import *
from utils.c_structs import c_update
from utils import safe_call 


#nls_lib=ctx.load_library_with_dir(os.environ['nls_wrapper_lib'])
#nls_lib=cdll.LoadLibrary(os.environ['nls_wrapper_lib'])
nls_lib=lib_env.load_lib_with_path('nls_pade_lib')
ssf_lib=lib_env.load_lib_with_path('ssf_lib')

'''
p_base_context_t  = c_void_p
pp_base_context_t  = POINTER(p_base_context_t)
'''

class base_context_t(Structure):
    _pack_=8
    _fields_ = [    
    ("weak_ref_handle",c_void_p),
    ("uuid",c_char*16),
    ("proc",c_void_p),  
    ("context",c_void_p)  
    ];
    
class nls_data_t(Structure):
    _pack_=8
    _fields_ = [
       ("ftype",c_char*8),
       ("pxx",c_void_p),
       
       ("rep",c_longlong),
       ("pp",c_longlong),
       
       ("g",c_double),
       ("alpha",c_double),
       ("alpha_i",c_double),
       
       ("dt",c_double),
       ("nt",c_longlong),
       
       ("n",c_longlong),
       ("m",c_longlong),
       
       ("omagnus",c_longlong),
       
       
       ("N",c_longlong),
       
       ("alpha_nl",c_double),
       ("alpha_nl_i",c_double),       
       
       ("w",c_double),       
       
       ("ncorr",c_longlong),
       ("flags",c_longlong),
       
       ("step_count",c_longlong),
       ("step_time",c_double),       
       ("masks",c_longlong),
       ("opts",c_void_p)
       ]    



p_base_context_t  = POINTER(base_context_t)
pp_base_context_t  = POINTER(p_base_context_t)




_create_nls_context_ex=nls_lib.create_nls_context_ex;
_create_nls_context_ex.rstype=c_uint32
_create_nls_context_ex.argtypes = (c_char_p,c_char_p,c_uint64,pp_base_context_t)

create_nls_context_ex=_create_nls_context_ex


_create_nls_context=nls_lib.create_nls_context;
_create_nls_context.rstype=c_uint32
_create_nls_context.argtypes = (c_char_p,c_uint64,pp_base_context_t)


_create_nls_magnus_context=nls_lib.create_nls_magnus_context;
_create_nls_magnus_context.rstype=c_uint32
_create_nls_magnus_context.argtypes = (c_char_p,c_uint64,pp_base_context_t)


_create_nls_magnus_exp_context=nls_lib.create_nls_magnus_exp_context;
_create_nls_magnus_exp_context.rstype=c_uint32
_create_nls_magnus_exp_context.argtypes = (c_char_p,c_uint64,pp_base_context_t)


_create_nls_magnus_exp_dec_context=nls_lib.create_nls_magnus_exp_dec_context;
_create_nls_magnus_exp_dec_context.rstype=c_uint32
_create_nls_magnus_exp_dec_context.argtypes = (c_char_p,c_uint64,pp_base_context_t)



_create_ms_context=nls_lib.create_ms_context;
_create_ms_context.rstype=c_uint32
_create_ms_context.argtypes = (c_char_p,c_uint64,pp_base_context_t)






_create_ssf_context=ssf_lib.create_ssf_context;
_create_ssf_context.rstype=c_uint32
_create_ssf_context.argtypes = (c_char_p,c_uint64,pp_base_context_t)


try:
    _create_ssf_ms_context=ssf_lib.create_ssf_ms_context;
    _create_ssf_ms_context.rstype=c_uint32
    _create_ssf_ms_context.argtypes = (c_char_p,c_uint64,pp_base_context_t)
except Exception:
    pass



_call_nls_context=nls_lib.call_nls_context;
_call_nls_context.rstype=c_uint32
_call_nls_context.argtypes = (p_base_context_t,c_char,c_void_p)


_release_nls_context=nls_lib.release_context;
_release_nls_context.rstype=c_uint32
_release_nls_context.argtypes = (p_base_context_t,)



norm=np.linalg.norm



    
'''
def nls_context_create(N,ftype=b'double'):
    ftype=to_bytes(ftype);
    p=p_base_context_t()
    err=_create_nls_context(ftype,N,byref(p));
    if err:
        raise Exception('Bad nls context');
        
    return hcontext(p);
'''

def __nm_update__(ds,dfl):
    d={}
    d.update(dfl);
    d.update(ds);    
    
    
    
    if 'nm' in d:
        nm=d['nm']
        d['n'],d['m']=nm[0],nm[1];
        
    for n in ('n','m','nt'):
        if n in d:
            d[n]=int(d[n]);
        
    if 'dz' in ds:
        d['dt']=ds['dz']
        
    return d;
    
    
class NLS_base(object):
    
    def __init__(this,N,dt=1,nt=1,nm=[8,8],g=1,alpha=0,ftype='double',factory=_create_nls_context,mp=1):
        
        #this.hctx=nls_context_create(N,ftype);
        #this._pctx=this.hctx.ptr;
        this.ftype=ftype=to_bytes(ftype);
        N=np.uint64(N);
        mp=np.uint64(mp);
        this._pctx=pctx=p_base_context_t()
        
        err=factory(ftype,N,byref(pctx));
        if err:
            raise Exception('Bad nls context');
        
        
        
        this.factory=factory
        
        #this._istart=icmd('start');
        #this._iinit=icmd('init');
        NN=N*mp;
        
        this._x=_x=np.zeros(NN,dtype=np.complex128);        
        this.dnls=dnls=nls_data_t(pxx=_x.ctypes.data, N=N);
        this._pnls=byref(dnls);
        this.reset(nm=nm,dt=dt,nt=nt);
        
        this.g=g;
        this.alpha=alpha;
        this.N=N;
        
        
            
        
    def reset0(this,nm=[8,8],dt=1,nt=1):
        dnls=this.dnls
        dnls.n,dnls.m=int(nm[0]),int(nm[1]);
        dnls.dt,dnls.nt= dt,nt;
        return _call_nls_context(this._pctx,b'i',this._pnls);
        #return this.hctx('init',this._pnls);
        
    def reset(this,**kw):
        d=__nm_update__(kw,{'nm':[8,8],'dt':1,'nt':1,'w':1.0,'step_time':np.nan,'step_count':-1})       
        c_update(this.dnls,**d);
        return _call_nls_context(this._pctx,b'i',this._pnls);
        #return this.hctx('init',this._pnls);
        
    def __del__(this):
        #print('del...')
        safe_call(lambda :_release_nls_context(this._pctx) );
        
        
        
    def conj(this):
        this._x[:]=np.conj(this._x);
        
    @property 
    def elapsed(this):
        return (this.dnls.step_time,this.dnls.step_count)
    
    @property 
    def x(this):
        return this._x;
    @x.setter
    def x(this,v):
        this._x[:]=v.flatten();
        
    @property 
    def g(this):
        return this.dnls.g;
    @g.setter
    def g(this,v):
        this.dnls.g=v;
        
        
    @property 
    def dz(this):
        return this.dnls.dt;
    @dz.setter
    def dz(this,v):
        this.dnls.dt=v;        

    @property 
    def alpha(this):
        return this.dnls.alpha+1j*this.dnls.alpha_i;
    @alpha.setter
    def alpha(this,v):
        this.dnls.alpha,this.dnls.alpha_i=np.real(v),np.imag(v);
        
        
    
    def __call__(this,rep=1,pp=True):
        d=this.dnls;
        d.rep,d.pp=rep,pp;
        err=this.step();
        if err:
            raise Exception('NLS solver error');
        return this._x;
    
    def step(this):
        return _call_nls_context(this._pctx,b's',this._pnls);
        #return call_context(this._pctx,this._istart,this._pnls);
        


class NLS_ex(NLS_base):    
    def __init__(this,solver_name,N,dt=1,nt=1,nm=[8,8],g=1,alpha=0,ftype='double'):
        
        sn=solver_name.upper();
        if sn=='SSF':
            factory=_create_ssf_context;        
        else:                    
            bname=to_bytes(solver_name);            
            factory=lambda ftype,N,pctx : _create_nls_context_ex(bname,ftype,N,pctx);
            
        super(NLS_ex,this).__init__(N,dt,nt,nm,g,alpha,ftype,factory=factory);

        
        

class NLS(NLS_base):    
    def __init__(this,N,dt=1,nt=1,nm=[8,8],g=1,alpha=0,ftype='double'):
        #,factory=_create_nls_context
        super(NLS,this).__init__(N,dt,nt,nm,g,alpha,ftype,factory=_create_nls_context);

class NLS_PM(NLS_base):    
    def __init__(this,N,dt=1,nt=1,nm=[8,8],g=1,alpha=0,ftype='double'):
        #,factory=_create_nls_context
        super(NLS_PM,this).__init__(N,dt,nt,nm,g,alpha,ftype,factory=_create_nls_magnus_context);
        
class NLS_PME(NLS_base):    
    def __init__(this,N,dt=1,nt=1,nm=[8,8],g=1,alpha=0,ftype='double'):
        #,factory=_create_nls_context
        super(NLS_PME,this).__init__(N,dt,nt,nm,g,alpha,ftype,factory=_create_nls_magnus_exp_context);
    

class NLS_PMED(NLS_base):    
    def __init__(this,N,dt=1,nt=1,nm=[8,8],g=1,alpha=0,ftype='double'):
        #,factory=_create_nls_context
        super(NLS_PMED,this).__init__(N,dt,nt,nm,g,alpha,ftype,factory=_create_nls_magnus_exp_dec_context);


class SSF(NLS_base):    
    def __init__(this,N,dt=1,nt=1,nm=[8,8],g=1,alpha=0,ftype='double'):
        #,factory=_create_nls_context
        super(SSF,this).__init__(N,dt,nt,nm,g,alpha,ftype,factory=_create_ssf_context);
        
        

class MS(NLS_base):    
    def __init__(this,N,dt=1,nt=1,nm=[8,8],g=1,alpha=0,ftype='double',fSSF=False,**kwd):
        
        pade_factory=kwd.get('pade_factory',_create_ms_context);
        ssf_factory=kwd.get('ssf_factory',_create_ssf_ms_context);
        
        factory=ssf_factory if fSSF else pade_factory;        
        super(MS,this).__init__(N,dt,nt,nm,g,alpha,ftype,factory=factory,mp=2);
        
    @property 
    def xy(this):
        return this._x.reshape((2,this.N));
    @xy.setter
    def xy(this,v):
        this._x[:]=v.flatten();
    @property 
    def x(this):
        return this._x[0:this.N];
    @x.setter
    def x(this,v):
        this._x[0:this.N]=v.flatten();
    @property 
    def y(this):
        return this._x[this.N:];
    @y.setter
    def y(this,v):
        this._x[this.N:]=v.flatten();
    
    def __call__(this,rep=1,pp=0b11):
        super(MS,this).__call__(rep,pp);
        return this.xy;




def MS_Pade(N,dt=1,nt=1,nm=[8,8],g=1,alpha=0,ftype='double'):
    return MS(N,dt,nt,nm,g,alpha,ftype,False);


def MS_SSF(N,dt=1,nt=1,nm=[8,8],g=1,alpha=0,ftype='double'):
    return MS(N,dt,nt,nm,g,alpha,ftype,True);
    
    
# ======= tests ===============

def test0(N=16*1024):
    
    import copy 
    from utils  import tic,toc
    norm=np.linalg.norm
    
    
     
    h=nls_context_create(N,ftype='q');
    data=nls_data_t(n=8,m=8,dt=1);    
    h('init',byref(data))
    x0=np.zeros(N,dtype=np.complex128)
    x0[0]=1;
    x=copy.copy(x0)
    
    data.pxx=x.ctypes.data
    data.g=0;
    data.pp=1;
    data.rep=100;    
    
    tic('start...');
    h('start',byref(data))
    toc('end:')
    
    print(norm(x),norm(x0))
    

    pass


def test1(N,g=1,rep=100,dt=1,nt=1,pp=True,ftype='d',nm=[8,-1.]):
    
    import copy 
    from utils  import tic,toc
    norm=np.linalg.norm
    
    
     
    
    #data=nls_data_t(n=8,m=-1,dt=1);    
    #h('init',byref(data))
    x0=np.zeros(N,dtype=np.complex128)
    x0[0]=1;
    
    
    nls=NLS(N,dt=dt,nt=nt,ftype=ftype,g=g,nm=nm);
    
    nls.x=x0;
    
    
    
    
    
    
    tic('start...');
    nls(rep=rep,pp=pp)
    toc('end:')
    
    x=nls.x;
    
    print(norm(x),norm(x0))
    
    nls.reset(dt=-1)
    tic('start...');
    nls(rep=rep,pp=pp)
    toc('end:')
    
    print(norm(x),norm(x0))

    return x


def test2(N,g=1,rep=100,dt=1,nt=1,pp=True,ftype='d',nm=[8,-1.]):
    
    import copy 
    from utils  import tic,toc
    norm=np.linalg.norm
    
    
     
    
    #data=nls_data_t(n=8,m=-1,dt=1);    
    #h('init',byref(data))
    x0=np.zeros(N,dtype=np.complex128)
    x0[0]=1;
    
    
    nls=SSF(N,dt=dt,nt=nt,ftype=ftype,g=g,nm=nm);
    
    nls.x=x0;
    
    
    
    
    
    
    tic('start...');
    nls(rep=rep,pp=pp)
    toc('end:')
    
    x=nls.x;
    
    print(norm(x),norm(x0))
    
    nls.reset(dt=-1)
    tic('start...');
    nls(rep=rep,pp=pp)
    toc('end:')
    
    print(norm(x),norm(x0))

    return x



if __name__=='__main__':
    
    import gc
    from p23 import *
    from utils  import *
    from utils.derr2m  import derr2m
    norm=np.linalg.norm
    
    normm=lambda x: norm(x,np.inf)
    abs=np.abs
    
    N=16*1024;
    
    test1(N);
    test2(N);
    print('pid=',os.getpid())

    
    
    
    
    