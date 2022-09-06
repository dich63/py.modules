#

from ctypes import *
from parallel.parallel_wrapper import *
import context_wrapper.context_wrapper as ctx
from jsonrpc.jsonclass import  jsobject,ext_def,to_dict,to_obj
import copy
import types

import numbers
import numpy as np
from scipy import sparse as sp


class lc_t(Structure):
    _pack_=8
    _fields_ =[
    ("grain",c_longlong),
    ("type",c_char*16),    
    ("N",c_longlong),
    ("D",c_longlong),    
    ("c",c_void_p),
    ("xx",c_void_p),
    ("pout",c_void_p)
    ]

class lc_copy_t(Structure):
    _pack_=8
    _fields_ =[
    ("grain",c_longlong),
    ("bytelen",c_longlong),
    ("dest",c_void_p),
    ("src",c_void_p)   
    ]

class lc_1_form_t(Structure):
    _pack_=8
    _fields_ =[
    ("grain",c_longlong),
    ("type",c_char*16),    
    ("N",c_longlong),
    ("offset",c_longlong),
    ("offset_inc",c_longlong),    
    ("f",c_void_p),
    ("x",c_void_p),
    ("rout",c_void_p)   
    ]


def get_splib():
    p=os.environ['sparse_context_lib'];
    #p='O:\\__ss\\suitesparse-metis-for-windows\\bin\\x64\\Debug\\tpt_klu.dll'
    #p='O:\\__ss\\suitesparse-metis-for-windows\\bin\\x64\\release\\tpt_klu.dll'
    #p="V:/Projects/tbb_wrapper/x64/Release/tbb_wrapper.dll"
    '''
    try:
        pass
        #cdll.LoadLibrary(p+'/../tbb.dll');
    except Exception:
        pass
    '''
    #
    return ctx.load_library_with_dir(p)

SMP_ISL=0x0100
SMP_ISI=0x0000
SMP_COO= 0
SMP_CSC= 1
SMP_RE = 0x010000
SMP_CX= 0x000000


class klu_common_t(Structure):
    _pack_=8
    _fields_ = [
    ("tol",c_double),
    ("memgrow",c_double),
    ("initmem_amd",c_double),
    ("initmem",c_double),
    ("maxwork",c_double),

    ("btf",c_longlong),
    ("ordering",c_longlong),  
    ("scale",c_longlong),
    #("nz",c_longlong),

    ("user_order",c_void_p),  
    ("user_data",c_void_p),  

    ("halt_if_singular",c_longlong),

    ("status",c_longlong),  
    ("nrealloc",c_longlong),

    ("structural_rank",c_longlong),  

    ("numerical_rank",c_longlong),

    ("singular_col",c_longlong),

    ("noffdiag",c_longlong),

    ("flops",c_double),
    ("rcond",c_double),
    ("condest",c_double),
    ("rgrowth",c_double),
    ("work",c_double),

    ("memusage",c_longlong),
    ("mempeak",c_longlong)

    ];


class klu_params_t(Structure):
    _pack_=8
    _fields_ = [
    ("sp_context",c_void_p),  
    ("common",c_void_p)      
     ];    

class klu_params_solve_t(Structure):
    _pack_=8
    _fields_ = [
    ("n",c_longlong),  
    ("nrs",c_longlong),  
    ("pp_in",c_void_p),      
    ("pp_out",c_void_p)
     ];    

class sparse_matrix_gaxpy_t(Structure):
    _pack_=8
    _fields_ = [('x',c_void_p),('y',c_void_p) ];

class sparse_matrix_params_t(Structure):
    _pack_=8
    _fields_ = [
    ("flags",c_longlong),
    ("n",c_longlong),
    ("m",c_longlong),  
    ("nzmax",c_longlong),
    ("nz",c_longlong),
    ("i",c_void_p),  
    ("p",c_void_p),  
    ("x",c_void_p)  
    ];

class sparse_jet_polyval_t(Structure):
    _pack_=8
    _fields_ = [
    ("z",c_double*2),
    ("scale",c_double*2),
    ("pp_out",pp_base_context_t)    
    ];
    
    
class sparse_jet_gaxpy_t(Structure):
    _pack_=8
    _fields_ = [('pp_jet',c_void_p)
    ,('p_out',c_void_p) 
    ];
    
    
class sparse_jet_jet_t(Structure):
    _pack_=8
    _fields_ = [
    ("type",c_char*8),
    ("z",c_double*2),
    ("count",c_longlong),
    ("n",c_longlong),
    ("pp_f",c_void_p),
    ("pp_jet",c_void_p)       
    ];

    

class sparse_jet_gaxpy_jet_t(Structure):
    _pack_=8
    _fields_ = [
    ("type",c_char*8),
    ("z",c_double*2),
    ("pp_f",c_void_p),
    ("py",c_void_p),    
    ("pys",c_void_p),
    ("pbuf",c_void_p)
    ];

class sparse_jet_info_t(Structure):
    _pack_=8
    _fields_ = [
    ("count",c_longlong),
    ("n",c_longlong),
    ("nz",c_longlong),
    ("fse",c_longlong),
    ("pp_sp",pp_base_context_t)    
    ];

_create_spm_context=get_splib().create_spm_context
_create_spm_context.rstype=c_int32
_create_spm_context.argtypes = (c_longlong,pp_base_context_t)

_create_spm_context_v=get_splib().create_spm_context
_create_spm_context_v.rstype=c_int32
_create_spm_context_v.argtypes = (c_longlong,c_void_p)


_create_spm_jet=get_splib().create_spm_jet
_create_spm_jet.rstype=c_int32
_create_spm_jet.argtypes = (c_longlong,c_longlong,c_void_p,pp_base_context_t)


_klu_create_context=get_splib().klu_create_context
_klu_create_context.rstype=c_int32
_klu_create_context.argtypes = (c_longlong,c_void_p,c_void_p)#pp_base_context_t)


_klu_common_default=get_splib().klu_common_default
_klu_common_default.rstype=c_int32
_klu_common_default.argtypes = (c_void_p,)#pp_base_context_t)


spm_cmp_struct=get_splib().spm_cmp_struct
spm_cmp_struct.rstype=c_int32
spm_cmp_struct.argtypes = (c_void_p,c_void_p)#pp_base_context_t)

create_qp_context=get_splib().create_qp_context
create_qp_context.rstype=c_int32
create_qp_context.argtypes = (c_void_p,)#pp_base_context_t)


def smp_check_err(hr):
    if hr:
        raise Exception('spm_context: error:'+hex(hr))
    return 0;
    
_cplxtype=np.dtype('complex128');    

'''
def _get_sparse_from_spm_context(self):
    
    s=sparse_matrix_params_t();
    smp_check_err(call_context(self.ptr,ord('i'),addressof(s)));
    n=s.n;
    nzmax=s.nzmax;
    flags=s.flags;
    c_index= c_int64  if flags&SMP_ISL  else c_int32;
    c_type=type(c_index().value);
    #print(c_type)
    
    bi=(c_index*nzmax).from_address(s.i);
    bp=(c_index*(n+1)).from_address(s.p);    
    bx=(c_double*(2*nzmax)).from_address(s.x);
    
    hh=hcontext();
    bi.hhh=hh
    i=np.frombuffer(bi,dtype=c_type);
    p=np.frombuffer(bp,dtype=c_type);
    x=np.frombuffer(bx,dtype='complex128',count=nzmax);
    
    r=sp.csc_matrix((x,i,p),shape=(n,n),copy=False);
    return r;
'''    
class hspcontext(hcontext):
    
    def __init__(self,ptr=None,ptype=None, constructor=None,addref=False):
        
        hcontext.__init__(self,ptr,ptype, constructor,addref);
        #s=sparse_matrix_params_t();
        #smp_check_err(call_context(self.ptr,ord('i'),addressof(s)));
        #self.params=s;
        
        
    def __del__(self):
        try:
            hcontext.__del__(self)
        except Exception:
            pass
    def gaxpy(self,x,y=None):        
        return spm_gaxpy(self,x,y);
    def rescale(self,scale,asyn=False):
        zz=(c_double*2)();
        zz[0],zz[1]=scale.real,scale.imag;
        if asyn:            
            return hinvoke_context(self,'rescale',zz,byref=True);
        else:
            self('rescale',byref(zz));
            return hinvoke_context();
    @property    
    def params(self):
        s=sparse_matrix_params_t();
        smp_check_err(call_context(self.ptr,ord('i'),addressof(s)));
        return s;
        
    @property    
    def sm(self):
        s=self.params
        #s=self.params               
        #s=sparse_matrix_params_t();
        #smp_check_err(call_context(self.ptr,ord('i'),addressof(s)));                                    
        n=s.n;
        nzmax=s.nzmax;
        flags=s.flags;
        c_index= c_int64  if flags&SMP_ISL  else c_int32;
        c_type=type(c_index().value);
        #print(c_type)
        
        bi=(c_index*nzmax).from_address(s.i);
        bp=(c_index*(n+1)).from_address(s.p);    
        bx=(c_double*(2*nzmax)).from_address(s.x);
        
        
        bx.__holder__=self;
        
        i=np.frombuffer(bi,dtype=c_type);
        p=np.frombuffer(bp,dtype=c_type);
        x=np.frombuffer(bx,dtype='complex128',count=nzmax);
        
        r=sp.csc_matrix((x,i,p),shape=(n,n),copy=False);
        return r;



def spm_context(sm,fhuge=0,asyn=False):
      
    type_sm=type(sm)
    
#    if (type_sm is p_base_context_t)  or issubclass(type_sm,hcontext):
#         return sm;        
    
    if issubclass(type_sm,hspcontext):
        return sm;
        
    if issubclass(type_sm,hcontext):
        return hspcontext(sm.ptr,addref=True);        
         
    if type_sm in [np.ndarray,list]:
        sm=sp.csc_matrix(sm,dtype=_cplxtype); 
        
    dtype=sm.dtype;    
    ft=SMP_CX
    if dtype==np.dtype('double'):
        ft=SMP_RE;        
    elif dtype!=_cplxtype:
        raise Exception('spm_context: only complex,double type support ')
        
    params = sparse_matrix_params_t();
    
    fmt=sm.format;
    f= fmt=='csc';
    if not f :
        (f,sm)=(False,sm.tocoo()) if fmt=='coo' else (True,sm.tocsc());
              
    
    ff= SMP_CSC if f else SMP_COO
    
    
    params.flags= ff | ft;
    nz=sm.nnz;
    (params.n,params.m)=sm.shape;
    params.nzmax=nz;
    params.nz=-1 if ff else nz;
    params.x=sm.data.ctypes.data;
    if f:
        params.i=sm.indices.ctypes.data;
        params.p=sm.indptr.ctypes.data;
    else:
        params.i=sm.row.ctypes.data;
        params.p=sm.col.ctypes.data;
    
    #pc=p_base_context_t()
    #hc=hcontext(ptype=p_base_context_t);
    #smp_check_err(_create_spm_context(fhuge,byref(pc)))    
    #smp_check_err(_create_spm_context(fhuge,hc.byref()))    
    '''
    p=p_base_context_t();
    _create_spm_context_v(fhuge,addressof(p));
    hc=hspcontext(p);
    '''
    hc=hspcontext(ptype=p_base_context_t,constructor =lambda p : _create_spm_context(fhuge,p) );
    #hc.gaxpy=types.MethodType(spm_gaxpy,hc,hcontext);
    #hc.sm=types.MethodType(_get_sparse_from_spm_context,hc,hcontext);
    if asyn:
        hi=hinvoke_context(hc,'load',params,byref=True,links=(params,sm));
        return (hi,hc);
    else:
        #smp_check_err(call_context(hc.ptr,ord('l'),addressof(params)));
        smp_check_err(hc('load',addressof(params)));
        return hc;
    



    
    
    
    

class hjet_context(hcontext):
    def __init__(self,spl,fhuge=False,parallel=False):
        
        ls=[];
        #if parallel:
        #    g=pp_group();            

        g=pp_group(sync=not parallel);           
         
        def cast_number(n):
            if isinstance(n,numbers.Number):
                return [n]
            else:
                return n
                
            
        def cast_p(p):
            if issubclass(type(p),hspcontext):
                return (None,p);        
            if p is None:
                return (None,None);
            f=(type(p) in [np.ndarray,list]);
            if f and len(p)==0:
                return (None,None);
            if f or hasattr(p,'format'):                
                
                s=spm_context(p,fhuge,asyn=1);
                #ls.append(s[1]);
                ls.append(s);
                g(s[0]);
                return s;

                '''
                if parallel:
                    ls.append(s[1]);
                    g(s[0]);
                    return s;
                else:
                    ls.append(s);                        
                    return (None,s)
                '''
            #elif issubclass(type(p),hcontext):
            #    return p.address;
                
            return (None,p)
            
        ll=[cast_number(p) for p in spl];        
        
        ll=[cast_p(p) for p in ll];
        
        #if parallel:
        #    g.join();

        g.join();
        ll=[p[1].address for p in ll];
        
        self.count=count=len(ll);
        npl=np.array(ll,dtype=np.int64);
        #return 
        p=p_base_context_t();
        smp_check_err(_create_spm_jet(fhuge,count,npl.ctypes.data,byref(p)));
        #link_object(p,(ls,spl));
        hcontext.__init__(self,p);        

        self.link((ls,spl));

        
        ci=sparse_jet_info_t();        
        smp_check_err(self('i',byref(ci)));
        
        self.n=ci.n;
        self.nz=ci.nz;
        self.struct_equival=bool(ci.fse);
        self.fhuge=fhuge 
        
        #self.zm=self.polyval         
        #self.jet=self.make_jet 
        
    def make_jet(self,z,ff,jet,n=-1,asyn=True,offset=0):
        jj=sparse_jet_jet_t();
        
        t=b'r' if ff[0].dtype==np.float64 else b'c';
        
        pff=ptr_array(ff)
        pjet=ptr_array(jet,offset)
        
        jj.z[0],jj.z[1]=z.real,z.imag
        
        jj.type=t
        jj.n=n;
        jj.count=len(pjet);
        jj.pp_f=pff.ctypes.data;
        jj.pp_jet=pjet.ctypes.data;
        
        #print('pp:'+hex(jj.pp_jet)+':[0]'+hex(pjet[0]));
        if asyn:
            i=hinvoke_context(self,'jet',params=jj,byref=True,links=(ff,jet,pff,pjet))
            
        else:
            if self('jet',byref(jj)):        
                i=hinvoke_context(raw_context=p_base_context_t())
            else:
                i=hinvoke_context();
            
        return i;
        
        
    def polyval(self,z,scale=1.0,asyn=False):
        
        
        if type(z) in [np.ndarray,list,tuple]:
            r= [self.polyval(v,asyn) for v in z];
            if asyn:
                r0,r1=[],[]
                for p in r:
                    r0.append(p[0]);
                    r1.append(p[1]);
                return (r0,r1);
            else:
                return r;
                          
            
        #fh=self.fhuge
        
        #p=p_base_context_t()        
        sjp=sparse_jet_polyval_t()        
        sjp.z[0],sjp.z[1]=z.real,z.imag                
        sjp.scale[0],sjp.scale[1]=scale.real,scale.imag
        #sjp.pp_out=pointer(p);
        
        
        #hsz=hcontext(p)#ptype=p_base_context_t,constructor =lambda p :self(fh,'p',p));
        hsz=hspcontext(p_base_context_t())
        #hsz.gaxpy=types.MethodType(spm_gaxpy,hsz,hcontext);
        #hsz.sm=types.MethodType(_get_sparse_from_spm_context,hsz,hcontext);
        #hsz.byref().contents=p;
        sjp.pp_out=hsz.byref()
        #smp_check_err(self('p',addressof(sjp)))        
        #return hsz;
        
        if asyn:
            hinvk=self.create_invoke_context('p',sjp,byref=True,links=hsz);
            #hinvk=self.create_invoke_context('p',addressof(sjp));
            #link_object(hinvk,(hsz,sjp));
            return (hinvk,hsz);
        else:
            smp_check_err(self('p',addressof(sjp)))
            return (hinvoke_context(),hsz);

    

    def gaxpy_jet(self,z,fnd=None,y=None,ys=None,buf=None,freal=True,fys=False):
        
        class gaxpy_jet_object(object):
            def __call__(self):
                return self.invoker();
            pass
        
        o=gaxpy_jet_object()
        
        sgj=sparse_jet_gaxpy_jet_t();
        
        dtype=np.complex128
        N=self.n;
        count=self.count;
        
        if  not fnd is None:
            if not  fnd[0].dtype==np.float64:
                freal=False;            
        
        if freal:
            sgj.type=b'r';
        else:
            sgj.type=b'c'; 

            dtype=np.float64
        if fnd is None:
            fnd=[np.zeros(shape=(N,),dtype=dtype) for k in range(count) ];
          
        if y is None:
            y=np.zeros(shape=(N,),dtype=np.complex128)
            
        if fys and (ys is None):
            ys=np.zeros(shape=(N,),dtype=np.complex128)
        if buf is None:
            buf=np.zeros(shape=(N,),dtype=np.complex128)
            
        o.fnd,o.y,o.ys,o.buf,o.z,o.sgj=  fnd,y,ys,buf,z,sgj
        
        lp=np.empty(shape=(count,),dtype=np.uint64);
        for k in range(count-1):                       
            lp[k]=fnd[k].ctypes.data
        
        sgj.z[0],sgj.z[1]= z.real,z.imag;
        sgj.pp_f =lp.ctypes.data;
        sgj.py=y.ctypes.data;        
        sgj.pbuf=buf.ctypes.data;
        
        if fys:
            sgj.pys=ys.ctypes.data;
        else:
            sgj.pys=0;
            
        o.invoker=hinvoke_context(context=self,cmd='x',params=sgj,byref=True,links=(sgj,lp))       
        
        return o
        
    def __del__(self):
        try:
            hcontext.__del__(self)
        except Exception:
            pass


    zm=polyval;
    jet=make_jet 

    
def hjet_context_create(spl,fhuge=False,parallel=True):
       return spl if issubclass(type(spl),hjet_context) else  hjet_context(spl,fhuge,parallel);
    
    
    
def spm_gaxpy(pc,x,y=None):
    if issubclass(type(pc),hcontext):
        pc=pc.address;
    x=np.array(x,dtype='complex',copy=False);
    if y is None:
        y=np.zeros_like(x);
    g=sparse_matrix_gaxpy_t();
    (g.x,g.y)=(x.ctypes.data,y.ctypes.data);
    smp_check_err(call_context(pc,ord('g'),addressof(g)));
    return y;
    
    

    
g_LA_op=p_base_context_t();
ctx_check_err(create_LA_op(byref(g_LA_op)));

class LA_op_invoker(object):
    @staticmethod    
    def set_lc_t(xx,cc,y,grain=0,parallel=False):
        if not parallel:
            grain=np.int64(-1);
            
        lc=lc_t()
        
        lc.grain=grain;
        lc.D=len(xx)
        lc.N=y.size;
        
            
        lc.type=b'c' if y.dtype==np.complex128 else b'r';
        
        lc.pout=y.ctypes.data;        
        
        cc=np.array(cc,dtype=np.complex128);
        lc.c=cc.ctypes.data;        
        
        #xxp=[np.array(x,dtype=np.complex128) for x in xx]
        
        npl=np.array([x.ctypes.data for x in xx],dtype=np.int64);
        
        lc.xx=npl.ctypes.data;
        
        return (lc,(xx,cc,y,npl));
    @staticmethod    
    def set_lc_copy_t(dest,src,grain=0,parallel=False):
        if not parallel:
            grain=np.int64(-1);            
        lc=lc_copy_t()
        lc.grain=grain;
        lc.bytelen=bytelen=dest.size*dest.itemsize;
        lc.dest=dest.ctypes.data;        
        
        if not src is None:
            if not src.size*src.itemsize==bytelen:
                raise Exception('lc_copy_t: error');
            lc.src=src.ctypes.data;
            
        return (lc,(dest,src))   
        
        
    @staticmethod    
    def linspace(xx,cc,y,grain=0,parallel=False):        
        lc,links=LA_op_invoker.set_lc_t(xx,cc,y,grain,parallel)
        return hinvoke_context(g_LA_op,'l',lc,byref=True,links=links);
    @staticmethod    
    def sum(xx,y,grain=0,parallel=False):        
        lc,links=LA_op_invoker.set_lc_t(xx,None,y,grain,parallel)
        return hinvoke_context(g_LA_op,'s',lc,byref=True,links=links);
        
    @staticmethod        
    def copy(dest,src,grain=0,parallel=False):        
        lc,links=LA_op_invoker.set_lc_copy_t(dest,src,grain,parallel)
        return hinvoke_context(g_LA_op,'c',lc,byref=True,links=links);
    @staticmethod        
    def swap(dest,src,grain=0,parallel=False):        
        lc,links=LA_op_invoker.set_lc_copy_t(dest,src,grain,parallel)
        return hinvoke_context(g_LA_op,'x',lc,byref=True,links=links);
    @staticmethod        
    def zero(dest,grain=0,parallel=False):        
        lc,links=LA_op_invoker.set_lc_copy_t(dest,None,grain,parallel)
        return hinvoke_context(g_LA_op,'z',lc,byref=True,links=links);

        

def relink_icontext(p,pc):    
    h=pc.contents.weak_ref_handle
    r=ctx.link_context(p.contents.weak_ref_handle,h);
    ctx.release_context(h);
    return r;


def link_object(p,o):
    p=cast_context(p);
    hp=p.contents.weak_ref_handle
    if hp:
        h=ctx.wrap_context(o);
        ctx.link_context(hp,h);
        ctx.release_context(h);
      
def spm_gaxpy_context(pc,x,y):
    pi=p_invoke_context_t();
    create_invoke_context(pc,byref(pi));
    ctx.link_context(pi.contents.weak_ref_handle,pc.contents.weak_ref_handle);
    g=sparse_matrix_gaxpy_t()
    (g.x,g.y)=(x.ctypes.data,y.ctypes.data);
    link_object(pi,(g,x,y));
    pi.contents.icmd=ord('g');
    pi.contents.params=addressof(g);        
    return pi;


def lin_comb_context(c,xx,y_out,grain=0):
    
    lc=lc_t();
    D=lc.D=len(c);
    N=lc.N=len(y_out);
    
    if not (len(xx),len(xx[0]))==(D,N):
        raise Exception('lin_comb_context invalid args');
        
    c=np.array(c,dtype='complex128');
    
    t=[];
    for x in xx:
        t.append(x.ctypes.data);
        
    xxp=np.array(t,dtype=np.longlong);
    
    t= 'c' if y_out.dtype.name=='complex128' else 'f';
    
    lc.grain=grain;
    lc.c=c.ctypes.data
    lc.xx=xxp.ctypes.data;
    lc.pout=y_out.ctypes.data;    
    lc.type=t;
    
    oh=(lc,c,xx,xxp,y_out);
    h=ctx.wrap_context(oh);
    ret=p_invoke_context_t();
    create_invoke_context(g_LA_op,byref(ret));
    
    ret.contents.icmd=ord('l');
    ret.contents.params=addressof(lc);
    
    ctx.link_context(ret.contents.weak_ref_handle,h);
    ctx.release_context(h);
    
    return ret;
    
def cast_str(o):
    if type(o) in (str,):
        o=eval(o)
    return o;

def klu_context_factory(**opts):
    common=klu_common_t()
    smp_check_err(_klu_common_default(byref(common)))    
    #struct_def(common,opts.get('common',{}))
    struct_def(common,opts)
    fhuge=opts.get('fhuge',False);
    return lambda p :_klu_create_context(fhuge,byref(common),p)


    
class sparse_solver_invoker(object):

    def __init__(self,sA,opts={}):
        
        self.opts=opts=jsobject(ext_def(opts,{
        'fhuge':False,
        'common':{},        
        'solver_factory':klu_context_factory,
        'n':None,
        'parallel':True
        }))
        
        
        
        fhuge=opts.fhuge
        common=opts.common
        sA=spm_context(sA,fhuge)
        self.sA=sA;
        y=opts.y;
        
        if y is None:
            n=opts.n
            if n is None:
                n=sA.params.n;
            y=create_buf(n,fzero=True);
        else:
            n=y.size;
        self.n=n;
        self._y=y;
        self._cond=c_double(np.nan);

        self.g=None;
        #pp_group(not self.opts.parallel);
        
        opts.n=n;

        context=opts.context  
        if context is None:
            solver_factory=cast_str(opts.solver_factory)        
            context=hcontext(ptype=p_base_context_t,constructor =solver_factory(**common));
        else:
            context=hcontext(ptype=p_base_context_t,constructor =lambda p :context('clone',p) );    
            
        self.context=context;
        
        self.analyze=hinvoke_context(context,'analyze',sA);
        self.factorize=hinvoke_context(context,'factor',sA);      
        self.solve=hinvoke_context(context,'solve',y);
        self.econd=hinvoke_context(context,'econd',addressof(self._cond));


    
    def clone(self):
        return sparse_solver_invoker(self.sA,context=self.context,n=self.n);

    def solve_group(self,bb,xx=None,grain=1,parallel=None):

        context=self.context;

        nrs=len(bb);
        n=self.n;

        grain=np.long(grain);        
        if parallel is None:
            parallel=self.opts.parallel;

        if xx is None:
            xx=bb;

        '''
        if not  ((n==self.n) and (nrs,n)==xx.shape):
            raise Exception('sparse_solver_invoker: error: invalid n');
        '''

        count=nrs/grain;
        rem =nrs%grain;

        def get_invoker(off,pnrs):
            ps=klu_params_solve_t();
            ps.n=n;
            ps.nrs=pnrs;

            px=xx[off:(off+pnrs)];
            pb=bb[off:(off+pnrs)];

            lin=ptr_array(pb);
            lout=ptr_array(px);
            ps.pp_in=lin.ctypes.data;
            ps.pp_out=lout.ctypes.data;
            hinv=hinvoke_context(context,'xsolve',params=ps,byref=True,links=(lin,lout));
            return hinv;

        hh=();

        for k in range(count):
            off=k*grain;
            h=get_invoker(off,grain);
            hh+=(h,);

        if rem:
            off=count*grain;
            h=get_invoker(off,rem);
            hh+=(h,);

        hbatch=hinvoke_batch(hh,asyn=parallel,links=(bb,xx));
        return hbatch;


        


    @property
    def phase(self):
        p=c_int64();
        self.context('phase',byref(p))
        return p.value
        
    @property
    def cond(self):
        if np.isnan(self._cond.value):
            if self.econd():
                self._cond.value=np.inf;
        return self._cond.value
    @property
    def y(self):
        return self._y;
    @y.setter
    def y(self,v):
        self._y[:]=v

'''
class sps_invoker(sparse_solver_invoker):

    def __init__(self,sA,**opts):
        super(sps_invoker,self).__init__(sA,opts);
        self._ifs=hinvoke_batch((self.factorize,self.solve));
    def __call__(self,y):
        self.y=y;
        self._ifs();
        return self.y.copy();
    @property
    def handle(self):
        return self._ifs;
        
'''        
    