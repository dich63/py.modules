# -*- coding: utf-8 -*-


from NLS.nls_pade import *

class shuttle_data_t(Structure):
    _pack_=8
    _fields_ = [       
       
       
       ("az_r",c_double),
       ("az_i",c_double),
       
       ("pa",c_void_p),
       ("pb",c_void_p),
       ("pc",c_void_p),
       ("px",c_void_p),
       ("pxo",c_void_p),          
              
       ("N",c_longlong),
       ("M",c_longlong)
       ]    


_create_shuttle_context=nls_lib.create_shuttle_context;
_create_shuttle_context.rstype=c_uint32
_create_shuttle_context.argtypes = (c_char_p,c_uint64,pp_base_context_t)

_create_cycle_shuttle_context=nls_lib.create_cycle_shuttle_context;
_create_cycle_shuttle_context.rstype=c_uint32
_create_cycle_shuttle_context.argtypes = (c_char_p,c_uint64,pp_base_context_t)





class shuttle_base(object):
    
    def __init__(this,N,ftype='double',factory=_create_shuttle_context):
        
        #this.hctx=nls_context_create(N,ftype);
        #this._pctx=this.hctx.ptr;
        this.ftype=ftype=to_bytes(ftype);
        this.N=N=np.uint64(N);                
        this._pctx=pctx=p_base_context_t()
        
        err=factory(ftype,N,byref(pctx));
        if err:
            raise Exception('Bad nls context');           
          
        this.sd=sd=shuttle_data_t();
        this._psd=byref(sd);
        this.xo=xo=np.zeros(N,dtype=np.complex128);
        sd.pxo=xo.ctypes.data;
        
        
    def __call__(this,x,tri=None,az=0.0):
        sd=this.sd;
        x=np.array(x,copy=False,dtype=np.complex128);
        sd.px=x.ctypes.data;
        if tri is None:
            return _call_nls_context(this._pctx,b's',this._psd);
        else:
           a=np.array(tri[0],copy=False,dtype=np.complex128); 
           b=np.array(tri[1],copy=False,dtype=np.complex128);
           c=np.array(tri[2],copy=False,dtype=np.complex128);
           
           sd.az_r,sd.az_i=az.real,az.imag;
           
           sd.pa=a.ctypes.data;
           sd.pb=b.ctypes.data;
           sd.pc=c.ctypes.data;
           err= _call_nls_context(this._pctx,b'f',this._psd);
           if err:
               raise Exception('shuttle solver error');
           return this.xo
            
        
        

class cycle_shuttle(shuttle_base):
    def __init__(this,N,ftype='double'):
        super(cycle_shuttle,this).__init__(N,ftype,factory=_create_cycle_shuttle_context);

class shuttle(shuttle_base):
    def __init__(this,N,ftype='double'):
        super(shuttle,this).__init__(N,ftype,factory=_create_shuttle_context);
        

