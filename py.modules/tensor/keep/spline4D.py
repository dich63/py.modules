#
import numpy as np



from numpy.linalg  import matrix_rank as rank
from numpy.linalg  import det 
from numpy.linalg  import inv as inv_mat 
from numpy.linalg  import solve as lin_solve
from numpy.linalg  import lstsq,norm

from numpy  import einsum as tensor_op

def ML_meshgrid(*kw):
    
    xx= np.meshgrid(*kw,indexing='ij');
    


def cell4D(pt=(0.,0.,0.,0.),flat=True,full=False):
    
    def make_jet(x):
        x2=x*x;        
        return np.array([[1.,x,x2,x2*x],[0.,1,2.*x,3*x2],[0.,0,2.,6*x],[0.,0,0,6]])

    
    x0,x1,x2,x3=make_jet(pt[0]);
    y0,y1,y2,y3=make_jet(pt[1]);
    z0,z1,z2,z3=make_jet(pt[2]);
    t0,t1,t2,t3=make_jet(pt[3]);
    
    rule=b'i,j,k,t->ijkt';    
    
    c0=tensor_op(rule,t0,z0,y0,x0)     
     
    cx=tensor_op(rule,t0,z0,y0,x1)
    cy=tensor_op(rule,t0,z0,y1,x0)
    cz=tensor_op(rule,t0,z1,y0,x0)
    ct=tensor_op(rule,t1,z0,y0,x0)
    
    
    
    cxx=tensor_op(rule,t0,z0,y0,x2) 
    cxy=tensor_op(rule,t0,z0,y1,x1)
    cxz=tensor_op(rule,t0,z1,y0,x1)
    cxt=tensor_op(rule,t1,z0,y0,x1)
    
    cyy=tensor_op(rule,t0,z0,y2,x0) 
    #cyx=tensor_op(rule,t0,z0,y1,x1) 
    cyz=tensor_op(rule,t0,z1,y1,x0) 
    cyt=tensor_op(rule,t1,z0,y1,x0) 
    
    
    czz=tensor_op(rule,t0,z2,y0,x0)
    #czx=tensor_op(rule,t0,z1,y0,x1)
    #czy=tensor_op(rule,t0,z1,y1,x0)
    czt=tensor_op(rule,t1,z1,y0,x0)
    
    ctt=tensor_op(rule,t2,z0,y0,x0)
    #ctz=tensor_op(rule,t1,z1,y0,x0)
    #ctx=tensor_op(rule,t1,z0,y0,x1)
    #cty=tensor_op(rule,t1,z0,y1,x0)
    
    
    cxyzt=tensor_op(rule,t1,z1,y1,x1) 
    
    
    
        #,cxx,cxy,cxz,cxt,cxyz];    
    
    #cdd=cxx+cyy+czz;
    
    ctxy=tensor_op(rule,t1,z0,y1,x1) 
    ctxz=tensor_op(rule,t1,z1,y0,x1)
    ctyz=tensor_op(rule,t1,z1,y1,x0)
    cxyz=tensor_op(rule,t0,z1,y1,x1)                 
    
    
    if full:
        
        tc=[c0,cx,cy,cz,ct,cxy,cxz,cxt,cyz,cyt,czt,cxx,cyy,czz,ctt,ctxy,ctxz,ctyz,cxyz,cxyzt];        
    else:
        #tc=[c0,cx,cy,cz,ct,cxy,cxz,cxt,cyz,cyt,czt,cxx,cyy,czz,ctt,cxyzt];
        tc=[c0,cx,cy,cz,ct,cxy,cxz,cxt,cyz,cyt,czt,ctxy,ctxz,ctyz,cxyz,cxyzt];        
    
    if flat:
        tc=[c.flatten() for c in tc]; 
    

    return tc;
       

def find_index4D(r,X,Y,Z,T):
    x,y,z,t=r[0],r[1],r[2],r[3];
    
    f=(X[0]<=x) and (x<X[-1]) and (Y[0]<=y) and (y<Y[-1]) and (Z[0]<=z) and (z<Z[-1]) and (T[0]<=t) and (t<T[-1]); 
    if f:
        ix=np.searchsorted(X,x,side='right'); 
        iy=np.searchsorted(Y,y,side='right');
        iz=np.searchsorted(Z,z,side='right');
        it=np.searchsorted(T,t,side='right');
        return (True,(ix-1,iy-1,iz-1,it-1));
    else:
        return (False,(-1,-1,-1,-1));
    
def to_a(a,dtype=np.float64):
    return np.array(a,dtype=dtype);
    
    
class bulk_4D_16(object):
    
    def __init__(this,jets,XYZT,idxs,full=True):
        
        inx=lambda ix,iy,iz,it :tuple(to_a(idxs,dtype=int)+np.array((ix,iy,iz,it),dtype=int));
             
        t_=this; 
        [X,Y,Z,T]=XYZT;
        
        #t_.ptc=ptc=XYZ[it(0,0,0)];
        
        def pxyzt(ix,iy,iz,it):
            (ix,iy,iz,it)=inx(ix,iy,iz,it);
            return to_a([X[ix],Y[iy],Z[iz],T[it]]);
                             
        t_.ptc=ptc=pxyzt(0,0,0,0);
                             
        
                             
        
                             
        def _pair(ix,iy,iz,it):            
            i=inx(ix,iy,iz,it);
            return (cell4D(pxyzt(ix,iy,iz,it)-ptc,True,full),jets[i]);
        
        (p0,j0)=_pair(0,0,0,0);
        (p1,j1)=_pair(0,0,0,1);
        (p2,j2)=_pair(0,0,1,0);
        (p3,j3)=_pair(0,0,1,1);
        
        (p4,j4)=_pair(0,1,0,0);
        (p5,j5)=_pair(0,1,0,1);
        (p6,j6)=_pair(0,1,1,0);
        (p7,j7)=_pair(0,1,1,1);
        
        
        (p8,j8)=_pair(1,0,0,0);
        (p9,j9)=_pair(1,0,0,1);
        (pA,jA)=_pair(1,0,1,0);
        (pB,jB)=_pair(1,0,1,1);
        
        (pC,jC)=_pair(1,1,0,0);
        (pD,jD)=_pair(1,1,0,1);
        (pE,jE)=_pair(1,1,1,0);
        (pF,jF)=_pair(1,1,1,1);
        
        
        
        F=np.vstack((p0,p1,p2,p3,p4,p5,p6,p7,p8,p9,pA,pB,pC,pD,pE,pF));        
        jj=np.hstack((j0,j1,j2,j3,j4,j5,j6,j7,j8,j9,jA,jB,jC,jD,jE,jF));
        
        
        #ccc=lin_solve(F,jj);
        if full:
            [ccc,err,rn,s]=lstsq(F,jj,rcond=None)  
            t_.err=err;
            t_.rank=rn;
        else:
            ccc=lin_solve(F,jj);
            
        t_.ccc=ccc.reshape((4,4,4,4));
        t_.ppp=np.ones((4,4),dtype=np.float64);
        #jet0=jets[idxs];
        t_.make=np.vectorize(lambda x,y,z,t: t_.make16((x,y,z,t)));
        
    def __call__(this,p):
        
        ccc,ppp=this.ccc,this.ppp;
        
        p=to_a(p)-this.ptc;
        #ppp[:,0]=1.0;
        ppp[:,1]=p;
        ppp[:,2]=p2=p*p;
        ppp[:,3]=p2*p;
        
        val=tensor_op(b'ijkn,i,j,k,n->',ccc,ppp[3],ppp[2],ppp[1],ppp[0]);
        
        return val;
        
        
def grid_setup_4D(f,X,Y,Z,T,axis=(0,1,2,3)):
    
    f=np.transpose(f,axis);  
  
    (X,Y,Z,T)=[ to_a(i).flatten()  for i in (X,Y,Z,T) ]
    
    px,py,pz,pt=[ i.argsort().tolist() for i in (X,Y,Z,T) ]
    
    X,Y,Z,T=X[px],Y[py],Z[pz],T[pt]        
    
    #px,py,pz=[pp[k] for k in axis ]    
    ##px,py,pz=to_a([px,py,pz])[list(axis)]
    
    f=f[:,:,:,px][:,:,py,:][:,pz,:,:][pt,:,:,:];
    
    return (f,X,Y,Z,T);
        
    
class tensor_spline_4D(object):
    
    
    def __init__(this,f,X,Y,Z,T,axis=(0,1,2,3),outval=np.nan,full=True):
        
        t_=this
        t_.outval=outval;              
        t_.axis=axis;
        t_.full=full;
        
        (f,X,Y,Z,T) = grid_setup_4D(f,X,Y,Z,T,axis)
        
        
        t_.XYZT=(X,Y,Z,T)
        
        #f=np.array(f);
        
        t_.cache=np.empty_like(f,dtype=object)
        t_.cache_count=0;
        
        faxis=(0,1,2,3)
        
        fx,fy,fz,ft=np.gradient(f,X,Y,Z,T,axis=faxis);
        
        
        fxx,fxy,fxz,fxt=np.gradient(fx,X,Y,Z,T,axis=faxis);
        fyx,fyy,fyz,fyt=np.gradient(fy,X,Y,Z,T,axis=faxis);
        fzx,fzy,fzz,fzt=np.gradient(fz,X,Y,Z,T,axis=faxis);
        ftx,fty,ftz,ftt=np.gradient(fz,X,Y,Z,T,axis=faxis);                              
        
        fxy=(fxy+fyx)/2.
        fxz=(fxz+fzx)/2.
        fxt=(fxt+fyt)/2.
        
        
        fyz=(fyz+fzy)/2.
        fyt=(fyt+fty)/2.
        fzt=(fzt+ftz)/2.
        
        
        
        del fyx,fzx,fzy,fty,ftz;
        
        '''
        ctxy=tensor_op(rule,t1,z0,y1,x1) 
        ctxz=tensor_op(rule,t1,z1,y0,x1)
        ctyz=tensor_op(rule,t1,z1,y1,x0)
        cxyz=tensor_op(rule,t0,z1,y1,x1) 
        '''
        
        fxyz,fxyt=np.gradient(fxy,Z,T,axis=(faxis[2],faxis[3]));
        fxzt=np.gradient(fxz,T,axis=faxis[3]);
        fyzt=np.gradient(fyz,T,axis=faxis[3]);  
        
        fxyzt=np.gradient(fxyz,T,axis=faxis[3]);   
        
        '''
        #full
        tc=[c0,cx,cy,cz,ct,cxy,cxz,cxt,cyz,cyt,czt,cxx,cyy,czz,ctt,ctxy,ctxz,ctyz,cxyz,cxyzt];        
        #256
        tc=[c0,cx,cy,cz,ct,cxy,cxz,cxt,cyz,cyt,czt,ctxy,ctxz,ctyz,cxyz,cxyzt];
        '''
        #
        if full:
            ts=(f ,fx,fy,fz,ft,fxy,fxz,fxt,fyz,fyt,fzt,fxx,fyy,fzz,ftt,fxyt,fxzt,fyzt,fxyz,fxyzt);        
        else:
            ts=(f ,fx,fy,fz,ft,fxy,fxz,fxt,fyz,fyt,fzt,fxyt,fxzt,fyzt,fxyz,fxyzt);        
        #ts=(f ,fx,fy,fz,ft,fxx,fxy,fxz,fxt,fyt,fyy,fyz,fyt,fzz,ftt,fxyt,fxzt,fyzt,fxyz,fxyzt);
        
        jet=np.transpose(ts,(1,2,3,4,0));
        t_.jet=jet;
        #t_.jet8=jet[:,:,:,0:8];     
        
        t_.pXYZT=np.transpose([X,Y,Z,T],(1,0));
        
        t_.make=np.vectorize(lambda x,y,z,t: t_.make16((x,y,z,t)));
        t_.pt=np.array([0.0,0,0,0]);
        
    @property    
    def pfill(this):        
        return (100.*this.cache_count)/this.cache.size;
        
    def make_projective(this,U4x4=None):
        if U4x4 is None:
            return this.make;
        else:                
            U4x4=to_a(U4x4);
            def mt(x,y,z,t):
                pt=U4x4@to_a([x,y,z,1.]);
                this.pt[0:3]=pt[0:3]/pt[3];
                this.pt[3]=t;
                return this.make16(pt);
            
            return np.vectorize(mt);
        
    
    def __call__(this,x,y,z,t):
        return this.make(x,y,z,t)
    
    def reset_cache(this):
        this.cache[:,:,:,:]=None;
        this.cache_count=0;
        
    def make16(this,p):
        t_=this;
        
        cache=t_.cache;
        
        #p=(x,y,z);
        (f,ixs)=find_index4D(p,*t_.XYZT)
        if f:
            bulk =cache[ixs];
            if bulk is None:
                cache[ixs]=bulk= bulk_4D_16(t_.jet,t_.XYZT,ixs,full=t_.full);
                t_.cache_count+=1;
                
            val=bulk(p);         
        else:
            val=t_.outval;
        
        return val;
    
## ----------------------- end --------------
def test0(num=64):
    t=np.linspace(-2,2, num=num)
    Y,Z,X=t,t,-t;
    [YY,ZZ,XX]=np.meshgrid(Y,Z,X)
    F=XX*XX*XX + YY*YY*YY + ZZ*ZZ*ZZ
    F=F*F
    ts=tesor_spline_3D(F,X,Y,Z)
    return ts

def test4Dcube(f=False):      
      
    cc=[]
    def _cube(x,y,z,t,pp):        
        c=cell4D(to_a((x,y,z,t)),1,f)
        pp+=c;
        return c;
    
    cube=lambda x,y,z,t: _cube(x,y,z,t,cc);
    
    r=cube(0,0,0,0);    
    cube(0,0,0,1);
    cube(0,0,1,0);
    cube(0,0,1,1);
    
    cube(0,1,0,0);
    cube(0,1,0,1);
    cube(0,1,1,0);
    cube(0,1,1,1);
    
    cube(1,0,0,0);    
    cube(1,0,0,1);
    cube(1,0,1,0);
    cube(1,0,1,1);
    
    cube(1,1,0,0);
    cube(1,1,0,1);
    cube(1,1,1,0);
    cube(1,1,1,1);
    
    
    
    
    return cc
    
    
    
if __name__=='__main__':
    
    from numpy import *
    from utils import *
    
    c=cell4D((0,0,0,0))
    
    c=to_a(test4Dcube(0));rc=rank(c);rc
    print('c',c.shape ,' rank=',rc)
    
    
    Np=45
    t=np.arange(Np,dtype=np.double)
    #Y,Z,X=0.1*t,10*t,t;
    Y,Z,X,T=t,t,t,t;
    [YY,ZZ,XX,TT]=np.meshgrid(Y,Z,X,T)
    
    #f=XX*YY*ZZ*XX*ZZ
    f=XX
    f=XX+YY+ZZ
    tic()
    ts=tensor_spline_4D(f,X,Y,Z,T,full=0)
    toc('tesor_spline_4D  sec:')
    
    tic()
    ff=ts(0.5,20.5,0.5,0.5)
    toc('call first  sec:')
    tic()
    ff=ts(0.5,20.7,0.1,0.5)
    toc('call second  sec:')
    print(ff)
    
    '''
    tt=np.linspace(2.1,45, num=80);
    tz=np.linspace(2.3,23.4, num=3);
    [zz,yy,xx]=np.meshgrid(tz,tt,-tt);
    
    
    
    tic()
    ff=ts(xx,yy,zz);
    toc('call first  sec:')
    del ff
    tic()
    ff=ts(xx,yy,zz);
    toc('call second sec:')
    del ff
    
    
    
    pt=ts.make_projective();    
    tic()
    ff=pt(xx,yy,zz);
    toc('projective None sec:')
    
    del ff
    U=eye(4);
    pt=ts.make_projective(U);    
    tic()
    ff=pt(xx,yy,zz);
    toc('projective sec:')
    '''

        
         

         
        
        
            
        
        