#
import numpy as np



from numpy.linalg  import matrix_rank as rank
from numpy.linalg  import det 
from numpy.linalg  import inv as inv_mat 
from numpy.linalg  import solve as lin_solve

from numpy  import einsum as tensor_op


def cell4D(pt=(0.,0.,0.,0.),flat=True,full=False):
    
    def make_jet(x):
              
        return np.array([1.,x])

    
    x0=make_jet(pt[0]);
    y0=make_jet(pt[1]);
    z0=make_jet(pt[2]);
    t0=make_jet(pt[3]);
    
    rule=b'i,j,k,t->ijkt';    
    
    c0=tensor_op(rule,t0,z0,y0,x0)     
     
 
    
    
        #,cxx,cxy,cxz,cxt,cxyz];    
    
    #cdd=cxx+cyy+czz;
    
    if flat:
        tc=[c0.flatten()]; 
    else:
        tc=[c0];

    return tc;
       

def find_index4D(r,X,Y,Z,T):
    
    x,y,z,t=r[0],r[1],r[2],r[3];
    
    f=(X[0]<=x) and (x<X[-1]) and (Y[0]<=y) and (y<Y[-1]) and (Z[0]<=z) and (z<Z[-1])  and (T[0]<=t) and (t<T[-1]); 
    if f:
        ix=np.searchsorted(X,x,side='right'); 
        iy=np.searchsorted(Y,y,side='right');
        iz=np.searchsorted(Z,z,side='right');
        it=np.searchsorted(T,t,side='right')
        return (True,(ix-1,iy-1,iz-1,it-1));
    else:
        return (False,(-1,-1,-1,-1));
                
    
def to_a(a,dtype=np.float64):
    return np.array(a,dtype=dtype);
    
    
class bulk_4D_16(object):
    
    def __init__(this,jets,XYZ,idxs):
        
        inx=lambda ix,iy,iz,it :tuple(to_a(idxs,dtype=int)+np.array((ix,iy,iz,it),dtype=int));
             
        t_=this; 
        [X,Y,Z]=XYZ;
        
        #t_.ptc=ptc=XYZ[it(0,0,0)];
        
        def pxyzt(ix,iy,iz,it):
            (ix,iy,iz,it)=itf(ix,iy,iz,it);
            return to_a([X[ix],Y[iy],Z[iz],T[it]]);
                             
        t_.ptc=ptc=pxyz(0,0,0);
                             
        
                             
        
                             
        def _pair(ix,iy,iz,it):            
            i=inx(ix,iy,iz,it);
            return (cell4D(pxyzt(ix,iy,iz,it)-ptc),jets[i]);
        
        (p0,j0)=_pair(0,0,0);
        (p1,j1)=_pair(0,0,1);
        (p2,j2)=_pair(0,1,0);
        (p3,j3)=_pair(0,1,1);
        
        (p4,j4)=_pair(1,0,0);
        (p5,j5)=_pair(1,0,1);
        (p6,j6)=_pair(1,1,0);
        (p7,j7)=_pair(1,1,1);
        
        
        
        F=np.vstack((p0,p1,p2,p3,p4,p5,p6,p7));        
        jj=np.hstack((j0,j1,j2,j3,j4,j5,j6,j7));
        
        
        ccc=lin_solve(F,jj);
        
        t_.ccc=ccc.reshape((4,4,4));
        t_.ppp=np.ones((3,4),dtype=np.float64);
        #jet0=jets[idxs];
        
    def __call__(this,p):
        
        ccc,ppp=this.ccc,this.ppp;
        
        p=to_a(p)-this.ptc;
        #ppp[:,0]=1.0;
        ppp[:,1]=p;
        ppp[:,2]=p2=p*p;
        ppp[:,3]=p2*p;
        
        val=tensor_op(b'ijk,i,j,k->',ccc,ppp[2],ppp[1],ppp[0]);
        
        return val;
        
        
def grid_setup(f,X,Y,Z,axis=(0,1,2)):
    
    f=np.transpose(f,axis);  
  
    (X,Y,Z)=[ to_a(i).flatten()  for i in (X,Y,Z) ]
    
    px,py,pz=[ i.argsort().tolist() for i in (X,Y,Z) ]
    
    X,Y,Z=X[px],Y[py],Z[pz]        
    
    #px,py,pz=[pp[k] for k in axis ]    
    ##px,py,pz=to_a([px,py,pz])[list(axis)]
    
    f=f[:,:,px][:,py,:][pz,:,:];
    
    return (f,X,Y,Z);
        
    
class tensor_line_4D(object):
    
    
    def __init__(this,f,X,Y,Z,axis=(0,1,2),outval=np.nan):
        
        t_=this
        t_.outval=outval;              
        t_.axis=axis;
        
        (f,X,Y,Z) = grid_setup(f,X,Y,Z,axis)
        
        
        t_.XYZ=(X,Y,Z)
        
        f=np.array(f);
        
        t_.cache=np.empty_like(f,dtype=object)
        t_.cache_count=0;
        
        faxis=(0,1,2)
        
        fx,fy,fz=np.gradient(f,X,Y,Z,axis=faxis);
        
        fxx,fxy,fxz=np.gradient(fx,X,Y,Z,axis=faxis);
        fyx,fyy,fyz=np.gradient(fy,X,Y,Z,axis=faxis);
        fzx,fzy,fzz=np.gradient(fz,X,Y,Z,axis=faxis);
        
        fxy=(fxy+fyx)/2.
        fxz=(fxz+fzx)/2.
        fyz=(fyz+fzy)/2.
        
        del fyx,fzx,fzy
        
        fxyx,fxyy,fxyz=np.gradient(fxy,X,Y,Z,axis=faxis);
        fxzx,fxzy,fxzz=np.gradient(fxz,X,Y,Z,axis=faxis);
        fyzx,fyzy,fyzz=np.gradient(fyz,X,Y,Z,axis=faxis);
        
        fxyz=(fxyz+fxzy+fyzx)/3.;
        
        del fxyx,fxyy,fxzx,fxzy,fxzz,fyzx,fyzy,fyzz;
        #ddf=fxx+fyy+fzz;
        
        
        
        
        
        jet=np.transpose((f,fx,fy,fz,fxy,fxz,fyz,fxyz,fxx,fyy,fzz),(1,2,3,0));
        t_.jet=jet;
        t_.jet8=jet[:,:,:,0:8];     
        
        t_.pXYZ=np.transpose([X,Y,Z],(1,0));
        t_.make=np.vectorize(lambda x,y,z: t_.make8((x,y,z)));
        
    @property    
    def pfill(this):        
        return (100.*this.cache_count)/this.cache.size;
        
    def make_projective(this,U4x4=None):
        if U4x4 is None:
            return this.make;
        else:                
            U4x4=to_a(U4x4);
            def mt(x,y,z):
                pt=U4x4@to_a([x,y,z,1.]);
                pt=pt[0:3]/pt[3];
                return this.make8(pt);
            
            return np.vectorize(mt);
        
    
    def __call__(this,x,y,z):
        return this.make(x,y,z)
        
    def make8(this,p):
        t_=this;
        
        cache=t_.cache;
        
        #p=(x,y,z);
        (f,ixs)=find_index3D(p,*t_.XYZ)
        if f:
            bulk =cache[ixs];
            if bulk is None:
                cache[ixs]=bulk= bulk_3D_8(t_.jet8,t_.XYZ,ixs);
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
    Np=256
    t=np.arange(Np,dtype=np.double)
    Y,Z,X=0.1*t,10*t,t;
    Y,Z,X=t,t,-t;
    [YY,ZZ,XX]=np.meshgrid(Y,Z,X)
    
    c=to_a(test4Dcube(1))
    rc=rank(c)
    print('c',c.shape ,' rank=',rc)
    
    
    '''
    #f=XX*YY*ZZ*XX*ZZ
    f=XX
    f=XX+YY+ZZ
    tic()
    ts=tesor_spline_3D(f,X,Y,Z)
    toc('tesor_spline_3D  sec:')
    
    ff=ts(-0.5,0.5,0.5)
    print(ff)
    
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
        
         

         
        
        
            
        
        