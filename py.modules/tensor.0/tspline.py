#
import numpy as np



from numpy.linalg  import matrix_rank as rank
from numpy.linalg  import det 
from numpy.linalg  import inv as inv_mat 
from numpy.linalg  import solve as lin_solve



def tsp_init(this,D=2,deg=3):            
    return np.ones([deg+1 for k in range(D)])*np.inf;


def idx3p (ii,kk=(0,0,0)):
    return tuple(np.array(ii,dtype=int)+np.array(kk,dtype=int))


def   tesor_spline_3D_skeleton():
    n0=np.ones(4,dtype=np.float64);
    #n1=np.array([1.,2.,3.,0.]);
    #n2=np.array([2.,6.,0.,0.]);
    
    n1=np.array([0.,1.,2.,3.]);
    n2=np.array([0.,0.,2.,6.]);
     
    c0=np.einsum("i,j,k->ijk",n0,n0,n0)
     
     
    cx=np.einsum("i,j,k->ijk",n0,n0,n1)
    cy=np.einsum("i,j,k->ijk",n0,n1,n0)
    cz=np.einsum("i,j,k->ijk",n1,n0,n0)
     
    cxy=np.einsum("i,j,k->ijk",n0,n1,n1)
    cxz=np.einsum("i,j,k->ijk",n1,n0,n1)
    czy=np.einsum("i,j,k->ijk",n1,n1,n0)
     
    cxx=np.einsum("i,j,k->ijk",n0,n0,n2)
    cyy=np.einsum("i,j,k->ijk",n0,n2,n0)
    czz=np.einsum("i,j,k->ijk",n2,n0,n0)
     
    cdd=cxx+cyy+czz;
    
    
    tc=(c0,cx,cy,cz,cdd,cxy,cxz,czy,cxx,cyy,czz);    
    
    #np.einsum("ijk,i,j,k->ijk",cx,x,x,x)    

    return tc;

__tesor_spline_3D_skeleton_64x64__=tesor_spline_3D_skeleton()[0:8];


def jet3D(x):
    x2=x*x;
    x3=x2*x;
    return np.array([[1.,x,x2,x3],[0.,1,2.*x,3*x2],[0.,0,2.,6*x],[0.,0,0,6*x]])

def get3Dcc11(pt,flat=True):
    
    x0,x1,x2,x3=jet3D(pt[0]);
    y0,y1,y2,y3=jet3D(pt[1]);
    z0,z1,z2,z3=jet3D(pt[2]);
    
    c0=np.einsum("i,j,k->ijk",z0,y0,x0)
     
     
    cx=np.einsum("i,j,k->ijk",z0,y0,x1)
    cy=np.einsum("i,j,k->ijk",z0,y1,x0)
    cz=np.einsum("i,j,k->ijk",z1,y0,x0)
     
    cxy=np.einsum("i,j,k->ijk",z0,y1,x1)
    cxz=np.einsum("i,j,k->ijk",z1,y0,x1)
    czy=np.einsum("i,j,k->ijk",z1,y1,x0)
     
    cxx=np.einsum("i,j,k->ijk",z0,y0,x2)
    cyy=np.einsum("i,j,k->ijk",z0,y2,x0)
    czz=np.einsum("i,j,k->ijk",z2,y0,x0)
    
    cxyz=np.einsum("i,j,k->ijk",z1,y1,x1) 
    #cdd=cxx+cyy+czz;
    
    
    tc=[c0,cx,cy,cz,cxy,cxz,czy,cxx,cyy,czz,cxyz];    
    
    if flat:
        tc=[c.flatten() for c in tc];
    
    #np.einsum("ijk,i,j,k->ijk",cx,x,x,x)    

    return tc;

def get3Dcc8(pt,flat=True):
    tc=get3Dcc11(pt,flat);
    #tc[7]+= tc[8] +  tc[9]+  tc[10]
    tc[7]= tc[10];

    return tc[0:8];
    

    
    

def get_cc8(x,y,z,flat=True):
    
    global __tesor_spline_3D_skeleton_64x64__
    
    def pd3(x):
        x2=x*x;
        return np.array([1.,x,x2,x*x2]);
    
    xx,yy,zz=pd3(x),pd3(y),pd3(z);
    if flat:
        cc8=[ np.einsum("ijk,i,j,k->ijk",c,zz,yy,xx).flatten() for c in __tesor_spline_3D_skeleton_64x64__ ];
    else:
        cc8=[ np.einsum("ijk,i,j,k->ijk",c,zz,yy,xx) for c in __tesor_spline_3D_skeleton_64x64__ ];
    
    return cc8;



def getSQ_cc8(sh=(0.,0,0)):
    ad=np.add;
    c0=get3Dcc8(ad((0.,0.,0.),sh));
    c1=get3Dcc8(ad((0.,0.,1.),sh));
    c2=get3Dcc8(ad((0.,1.,0.),sh));
    c3=get3Dcc8(ad((0.,1.,1.),sh));
    
    c4=get3Dcc8(ad((1.,0.,0.),sh));
    c5=get3Dcc8(ad((1.,0.,1.),sh));
    c6=get3Dcc8(ad((1.,1.,0.),sh));
    c7=get3Dcc8(ad((1.,1.,1.),sh));
    
    
    return np.vstack((c0,c1,c2,c3,c4,c5,c6,c7));

getSQ_cc=getSQ_cc8

def getSQ_cc11():
    c0=get3Dcc11((0.,0.,0.));
    c1=get3Dcc11((0.,0.,1.));
    c2=get3Dcc11((0.,1.,0.));
    c3=get3Dcc11((0.,1.,1.));
    
    c4=get3Dcc11((1.,0.,0.));
    c5=get3Dcc11((1.,0.,1.));
    c6=get3Dcc11((1.,1.,0.));
    c7=get3Dcc11((1.,1.,1.));
    
    
    return np.vstack((c0,c1,c2,c3,c4,c5,c6,c7));

    
def find_index(x,xx):    
    
    if (xx[0]<=x) and (x<xx[-1]):
        return np.searchsorted(xx, x,side='right');
    else:
        return -1;
    
def find_index3D(r,X,Y,Z):
    x,y,z=r[0],r[1],r[2];
    
    f=(X[0]<=x) and (x<X[-1]) and (Y[0]<=y) and (y<Y[-1]) and (Z[0]<=z) and (z<Z[-1]); 
    if f:
        ix=np.searchsorted(X,x,side='right'); 
        iy=np.searchsorted(Y,y,side='right');
        iz=np.searchsorted(Z,z,side='right');
        return (True,(ix-1,iy-1,iz-1));
    else:
        return (False,(-1,-1,-1));
    
    
    
class bulk_3D_8(object):
    
    def __init__(this,jets,XYZ,idxs):
        
        it=lambda ix,iy,iz :tuple(np.array(idxs,dtype=int)+np.array((ix,iy,iz),dtype=int));             
        t_=this; 
        [X,Y,Z]=XYZ;
        
        #t_.ptc=ptc=XYZ[it(0,0,0)];
        
        def pxyz(ix,iy,iz):
            (ix,iy,iz)=it(ix,iy,iz);
            return np.array([X[ix],Y[iy],Z[iz]]);
                             
        t_.ptc=ptc=pxyz(0,0,0);
                             
        
                             
        
                             
        def _pair(ix,iy,iz):
            i=it(ix,iy,iz);            
            return (get3Dcc8(pxyz(ix,iy,iz)-ptc),jets[i]);
        
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
        
        p=np.array(p,dtype=np.float64)-this.ptc;
        #ppp[:,0]=1.0;
        ppp[:,1]=p;
        ppp[:,2]=p2=p*p;
        ppp[:,3]=p2*p;
        
        val=np.einsum("ijk,i,j,k->",ccc,ppp[2],ppp[1],ppp[0]);
        
        return val;
        
        
        
        
    
class tesor_spline_3D(object):
    def __init__(this,f,X,Y,Z,faxis=(2,1,0),outval=np.inf):
        
        t_=this
        t_.outval=outval;

        
        t_.XYZ=(X,Y,Z)=[ (lambda a : np.array(a,dtype=np.float64)  )(t)  for t in (X,Y,Z) ]
        
        
        f=np.array(f);
        
        t_.cache=np.empty_like(f,dtype=object)
        
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
        
        
    def __call__(this,x,y,z):
        t_=this;
        
        cache=t_.cache;
        
        p=(x,y,z);
        (f,ixs)=find_index3D(p,*t_.XYZ)
        if f:
            bulk =cache[ixs];
            if bulk is None:
                cache[ixs]=bulk= bulk_3D_8(t_.jet8,t_.XYZ,ixs);
                
            val=bulk(p);         
        else:
            val=t_.outval;
        
        return val;
    
    
def test0():
    t=np.arange(128,dtype=np.double)
    Y,Z,X=0.1*t,10*t,t;
    [YY,ZZ,XX]=np.meshgrid(Y,Z,X)
    
    
if __name__=='__main__':
    
    t=np.arange(64,dtype=np.double)
    Y,Z,X=0.1*t,10*t,t;
    Y,Z,X=t,t,t;
    [YY,ZZ,XX]=np.meshgrid(Y,Z,X)
    
    f=XX*YY*ZZ*XX*ZZ
    
    ts=tesor_spline_3D(f,X,Y,Z)
    ff=ts(0.5,0.5,0.5)
    print(ff)
    

        
        
         

         
        
        
            
        
        