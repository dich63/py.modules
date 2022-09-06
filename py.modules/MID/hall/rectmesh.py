#

import numpy as np
from jsobj import *
from MID.hall.hsolver import *

def gen_def_mesh0(Nrx,Nry):
    Nx,Ny=Nrx+1,Nry+1
    nx=1j*np.ones((1,Nx));
    ny=np.ones((Ny,1));
    
    #np.reshape(range(NxNy),(Nx,Ny))
    vxs=ny*np.reshape(range(Nx),(1,Nx))+np.reshape(range(Ny),(Ny,1))*nx;
    ivxs=np.reshape(range(Nx*Ny),(Ny,Nx))
    rects=np.empty((Nrx*Nry,4),dtype=int);
    irects=np.reshape(range(Nrx*Nry),(Nry,Nrx))
    
    for i in range(Nrx):
        for j in range(Nry):
            rects[irects[j,i]]=[ivxs[j,i],ivxs[j,i+1],ivxs[j+1,i],ivxs[j+1,i+1]];
            
    evxs_v=np.empty((Nry,Nx,5),dtype=int);        
    irx=-2;
    ivx=0;
    f=0;  
    
    
    evx=np.empty((Nry+1,Nrx,2),dtype=int);
    erx=np.empty_like(evx);
    ebx=np.empty((Nry+1,Nrx),dtype=int);
    eix=np.reshape(range((Nry+1)*Nrx),(Nry+1,Nrx));   
    for iy in range(Nry+1):
        f=int((iy>0) & (iy<Nry));
        for ix in range(Nrx):
            ebx[iy,ix]=not f;
            evx[iy,ix]=[ivxs[iy,ix],ivxs[iy,ix+1]]
            if f:
                erx[iy,ix]=[irects[iy-1,ix],irects[iy,ix]]
            else:
                erx[iy,ix]=[-1,irects[iy,ix]]  if iy==0 else [irects[iy-1,ix],-1]
                    
            
            
    evy=np.empty((Nry,Nrx+1,2),dtype=int);
    ery=np.empty_like(evy);
    eby=np.empty((Nry,Nrx+1),dtype=int);  
    eiy=np.reshape(range(Nry*(Nrx+1)),(Nry,Nrx+1));      
    eiy+=(Nry+1)*Nrx;
    for ix in range(Nrx+1):
        f=int((ix>0) & (ix<Nrx));
        for iy in range(Nry):
            eby[iy,ix]=not f;
            evy[iy,ix]=[ivxs[iy,ix],ivxs[iy+1,ix]]
            if f:
                ery[iy,ix]=[irects[iy,ix-1],irects[iy,ix]] 
            else:
                ery[iy,ix]=[-1,irects[iy,ix]] if ix==0 else [irects[iy,ix-1],-1]

    
    
    
    
    return {
            'vxs':vxs,'rects':rects
            ,'irects':irects,'ivxs':ivxs
            ,'evx':evx,'erx':erx,'eix':eix
            ,'evy':evy,'ery':ery,'eiy':eiy
            }

    #vxs=np.array()
    pass

def e2d(a):
    s=a.shape;
    return a.reshape((int(np.prod(s[0:-1])),s[-1]));

def z_rescale(vx,Nrx,Nry,lt,br):
    def rsc(x,N,l,r):
        return l+((r-l)/N)*x;
    vr=rsc(vx.real,Nrx,lt.real,br.real)+1j*rsc(vx.imag,Nry,lt.imag,br.imag)
    return vr;     
   
def set_mu(mesh,mu=1,rect=[-1,1,-1,1]):
    pass
    
def reset_zc_rect(mesh):    
    zrc=mesh.zs[mesh.rects];
    zrc=np.sum(zrc,1)
    zrc*=0.25;    
    return zrc;



def gen_def_mesh1(Nrx=4,Nry=3,lt=-1-1j,br=1+1j):    
    rm=gen_def_mesh0(Nrx,Nry)
    o=jso(rm);
    r=jso()
    #r.rm=o;
    r.evx,r.erx,r.eix, r.evy,r.ery,r.eiy=[e2d(s) for s in (o.evx,o.erx,o.eix,o.evy,o.ery,o.eiy)]
    r.rects=o.rects
    r.mu=np.ones(r.rects.shape[0],dtype=np.float64);
    r.irects=o.irects;
    r.vxs=o.vxs.flatten();
    r.zs=z_rescale(r.vxs,Nrx,Nry,lt,br);
    r.ejx=np.zeros_like(r.eix,dtype=np.float64);
    r.ejy=np.zeros_like(r.eiy,dtype=np.float64);
    return (r,o)
    
def gen_def_mesh2(Nrx=4,Nry=3,domain=None):
    
    d=domain
    if d is None:
        lt,br=0+0j,Nrx+1j*Nry;
    else:        
        lt,br=d[0]+1j*d[1],d[2]+1j*d[3];
        
    (r,o)=gen_def_mesh1(Nrx,Nry,lt,br);
    
    m=jso();
    
    evx,erx,eix,evy,ery,eiy,m.zs,m.rects,m.irects=r.evx,r.erx,r.eix, r.evy,r.ery,r.eiy,r.zs,r.rects,r.irects
    
    
    ne=evx.shape[0];
    '''
    
    evy=evy+ne;
    f=ery<0;
    ery=ery+ne;
    ery[f]=-1;
    eiy=eiy+ne;
    '''
    
    m.eix,m.eiy=eix,eiy;
    evxy=np.concatenate((evx,evy),0)
    erxy=np.concatenate((erx,ery),0)
    m.evxy,m.erxy=evxy,erxy
    #evxy=100*evxy;
    es=np.concatenate((evxy,erxy),1)
    m.edges=es.reshape((es.shape[0],2,2))
    #m.edges=es
    ox=1*np.ones(evx.shape[0],dtype=np.complex128);
    oy=1j*np.ones(evy.shape[0],dtype=np.complex128);
    m.ef=np.concatenate((ox,oy))
    
    
    
    
    m.mu=r.mu

    m.zrc=reset_zc_rect(m);
    
    return (m,r,o)
    
def make_rectangulation(rl,zl):
    pass


def make_interp_xy(mesh):
    m=mesh;
    nx,ny=m.irects.shape;
    iix,iiy=np.arange(nx),np.arange(ny);
    
    xc=m.zs[m.rects[m.irects[0,:]]]
    x=xc.real;
    xmin,xmax=x[0,0],x[-1,-1];

    yc=m.zs[m.rects[m.irects[:,0]]]
    y=yc.imag;
    ymin,ymax=y[0,0],y[-1,-1];

    xx=x[:,0];
    yy=y[:,0];

    

    

    

    return  (xx,yy,xmin,xmax,ymin,ymax)


def make_interp_index(mesh):

    irects=mesh.irects;
    (xx,yy,xmin,xmax,ymin,ymax)=make_interp_xy(mesh);    

    fb=lambda x,y : (xmin<x) and (x<xmax) and (ymin<y) and (y<ymax)

    def get_index(z):
        x,y=z.real,z.imag;
        if fb(x,y):
            ix=np.nonzero(xx<x)[0];
            iy=np.nonzero(yy<y)[0];
            return irects[iy[-1],ix[-1]];
            pass
        else:
            return -1;
    return get_index;

def zc_in_rects_mask(zc,rects):
    x,y=zc.real,zc.imag;
    fs=np.zeros_like(zc,dtype=bool);
    for r in rects:
       f=(r[0]<x) & (x< r[2]) & (r[1]<y) & ( y< r[3]);
       fs=fs | f;

    return fs;


    


def MeshData2mesh(meshdata):
    X, Y = np.meshgrid(meshdata.rl, meshdata.zl);
    Ny,Nx=X.shape;
    print('Ny=',Ny,'Nx=',Nx)
    mesh=gen_def_mesh2(Nx-1,Ny-1)[0];
    z=X+1j*Y;

    mesh.zs=zs=z.flatten();
    mesh.zrc=zrc=reset_zc_rect(mesh);
    e=mesh.edges[:,0];
    fv=mesh.ef.real==0.0;
    #fv=mesh.ef.imag==0.0;
    ze=zs[e];
    ze=np.sum(ze,1)/2;
    J=np.zeros_like(ze);
    mu=mesh.mu;
    for rg in meshdata.regions:
        f=zc_in_rects_mask(zrc,rg.rects);
        mu[f]=rg.mu;
        if rg.name=='gcoil1':
            f=zc_in_rects_mask(ze,rg.rects);
            f=fv & f;
            J[f]=1;
        if rg.name=='gcoil2':
            f=zc_in_rects_mask(ze,rg.rects);
            f=fv & f;
            J[f]=-1;           


    return mesh,J;




if __name__=='__main__':
    Nrx=4
    Nry=3
    r=gen_def_mesh2(Nrx,Nry)
    m=r[0]
    m=jso(m)
    print(m)
    