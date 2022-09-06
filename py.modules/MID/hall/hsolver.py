#
import numpy as np
from jsobj import *
from scipy import sparse as sp
import scipy.sparse.linalg as la
from numpy.polynomial.polynomial import polyval
import math


def node_points(zz,Np):
    n = Np + 1
    a = zz[0]
    d = (zz[1] - a) / n    
    z = [ a + d * k for k in range(1,n)]
    return z

def cri(z,Ndeg,vh):
    N = Ndeg + 1
    zn = np.empty(N,dtype=np.complex)
    t = vh
    zn[0] = t
    for k in range(1,N):
        t*=z
        zn[k] = t
    return (zn.real,zn.imag)


def get_nrows(Np):
    return 2 * len(node_points((0,1j),Np))

def make_edge_row(mesh,nedge=0,Np=1,Ndeg=-1):
    
    Ndeg = 4 * Np - 1 if Ndeg < 0 else Ndeg

    edge = mesh.edges[nedge]
    ef = mesh.ef[nedge]
    ip = edge[0]
    
    zp = mesh.zs[ip]
    ze = node_points(zp,Np)  
    Ntb = 2 * (Ndeg + 1)    
    ir = edge[1]
    [f0,f1] = ir >= 0
    ff = f0 and f1
    irb = Ntb * ir
    mus = mesh.mu
    
    #if f[0] and f[1] :
    data = []
    irows = []
    icols = []
    nn = np.arange(Ntb)
    no2 = np.ones(2 * Ntb,dtype=np.int)
    no = no2[0:Ntb]
    nrows = 0
    bcdn = 1.e-0
    for z in ze:
        (ar,ai) = cri(z,Ndeg,ef)
        
        
        if ff:
            mu0,mu1 = mus[ir[0]],mus[ir[1]]
            rm = mu0 / mu1
            mm = (mu0 + mu1) / 2
            
            
            
            data +=[ar,-ai,-ar,ai]

            if rm < 1:
                data +=[ai,ar,-ai * rm,-ar * rm] 
            else:
                rm = 1. / rm
                data +=[ai * rm,ar * rm,-ai,-ar] 



            
            
            
            inn0,inn1 = irb[0] + nn,irb[1] + nn              
            icols+=[inn0,inn1,inn0,inn1]
            
            nnr1 = nrows * no2 
            nnr2 = (nrows + 1) * no2
            irows+=[nnr1,nnr2]           
            
        elif f0:
            ar*=bcdn
            ai*=bcdn
            mu0 = mus[ir[0]]
            data +=[ar,-ai]                            
            data +=[ai,ar]                 
            inn0 = irb[0] + nn
            icols+=[inn0,inn0]
            
            nnr1 = nrows * no 
            nnr2 = (nrows + 1) * no
            irows+=[nnr1,nnr2]           
        else:
            ar*=bcdn
            ai*=bcdn
            #if ef.imag:
            #    ar*=0;
            

            mu1 = mus[ir[1]]
            data +=[ar,-ai]                            
            data +=[ai,ar]                 
            inn1 = irb[1] + nn
            icols+=[inn1,inn1]
            
            nnr1 = nrows * no 
            nnr2 = (nrows + 1) * no
            irows+=[nnr1,nnr2]           
            
        nrows+=2
    
    irows = np.concatenate(irows,None)
    icols = np.concatenate(icols,None)
    data = np.concatenate(data,None)
    return (data,irows,icols,nrows)


def make_edge_row_c(mesh,nedge=0,Np=1,Ndeg=-1):
    
    Ndeg = 4 * Np - 1 if Ndeg < 0 else Ndeg

    edge = mesh.edges[nedge]
    ef = mesh.ef[nedge]
    ip = edge[0]
    
    zp = mesh.zs[ip]
    ze = node_points(zp,Np)  
    Ntb = 2 * (Ndeg + 1)    
    ir = edge[1]
    [f0,f1] = ir >= 0
    ff = f0 and f1
    irb = Ntb * ir
    mus = mesh.mu
    
    #if f[0] and f[1] :
    data = []
    irows = []
    icols = []
    nn = np.arange(Ntb)
    no2 = np.ones(2 * Ntb,dtype=np.int)
    no = no2[0:Ntb]
    nrows = 0
    bcdn = 1.e-0
    zrc = mesh.zrc

    #zc0,zc1=zrc[ir[0]],zrc[ir[1]]

    for z in ze:
        
        
        
        if ff:
            zc0,zc1 = zrc[ir[0]],zrc[ir[1]]

            (ar0,ai0) = cri(z - zc0,Ndeg,ef)
            (ar1,ai1) = cri(z - zc1,Ndeg,ef)

            mu0,mu1 = mus[ir[0]],mus[ir[1]]
            rm = mu0 / mu1
            mm = (mu0 + mu1) / 2
            bmm = 1 / mm

            bm0 = mm / mu0
            bm1 = mm / mu1

            #bmm=np.sqrt(bm0*bm1);

            data +=[ar0,-ai0,-ar1,ai1] 
            data +=[bm0 * ai0,bm0 * ar0,-bm1 * ai1,-bm1 * ar1] 
            
            '''
            #
            data +=[bmm*ar0,-bmm*ai0,-bmm*ar1,bmm*ai1]; 
            data +=[bm0*ai0,bm0*ar0,-bm1*ai1,-bm1*ar1]; 

            #data +=[ai0/mu0,ar0/mu0,-ai1/mu1,-ar1/mu1]; 

            #mu0,mu1=mu0/mm,mu1/mm
            #data +=[ai0*mu1,ar0*mu1,-ai1*mu0,-ar1*mu0]; 
            
            

            data +=[ar0,-ai0,-ar1,ai1];
            if rm<1:
                data +=[ai0,ar0,-ai1*rm,-ar1*rm]; 
            else:
                rm=1./rm;
                data +=[ai0*rm,ar0*rm,-ai1,-ar1]; 
            '''            

            
            
            
            inn0,inn1 = irb[0] + nn,irb[1] + nn              
            icols+=[inn0,inn1,inn0,inn1]
            
            nnr1 = nrows * no2 
            nnr2 = (nrows + 1) * no2
            irows+=[nnr1,nnr2]           
            
        elif f0:

            zc0 = zrc[ir[0]]
            (ar0,ai0) = cri(z - zc0,Ndeg,ef)

            ar0*=bcdn
            ai0*=bcdn
            mu0 = mus[ir[0]]
            data +=[ar0,-ai0]                            
            data +=[ai0,ar0]                 
            inn0 = irb[0] + nn
            icols+=[inn0,inn0]
            
            nnr1 = nrows * no 
            nnr2 = (nrows + 1) * no
            irows+=[nnr1,nnr2]           
        else:

            zc1 = zrc[ir[1]]
            (ar1,ai1) = cri(z - zc1,Ndeg,ef)

            ar1*=bcdn
            ai1*=bcdn
            #if ef.imag:
            #    ar*=0;
            

            mu1 = mus[ir[1]]
            data +=[ar1,-ai1]                            
            data +=[ai1,ar1]                 
            inn1 = irb[1] + nn
            icols+=[inn1,inn1]
            
            nnr1 = nrows * no 
            nnr2 = (nrows + 1) * no
            irows+=[nnr1,nnr2]           
            
        nrows+=2
    
    irows = np.concatenate(irows,None)
    icols = np.concatenate(icols,None)
    data = np.concatenate(data,None)
    return (data,irows,icols,nrows)

            
    
def make_edge_row_o(mesh,nedge=0,Np=1,Ndeg=-1):    

    o = jso({'nedge':nedge})
    (o.data,o.irows,o.icols,o.nrows) = make_edge_row(mesh,nedge,Np,Ndeg)
    return o
    
    
def make_J_edges(N,nrows):
    NN = N * nrows
    d = np.ones(NN,dtype=np.float64)
    r = np.arange(NN)
    c = np.arange(N).repeat(nrows)
    return sp.coo_matrix((d,(r,c)))

def make_J_edges2(N,nrows):    
    NN = N * nrows
    d = np.ones(NN,dtype=np.float64)
    d[0:-1:2] = 0
    r = np.arange(NN)
    c = np.arange(N).repeat(nrows)
    return sp.coo_matrix((d,(r,c)))


def make_J_edges1(N,nrows):    

    N,nrows,n2 = int(N),int(nrows),int(nrows / 2)        
    d = np.ones(N * n2,dtype=np.float64)    
    r = np.arange(1,N * nrows,2)
    c = np.arange(N).repeat(n2)

    return sp.coo_matrix((d,(r,c)))


def make_H_edges(mesh,Np=1,Ndeg=-1):
    mesh = to_jso(mesh)
    #sm=sp.csc_matrix(sm,dtype=_cplxtype);
    N = mesh.edges.shape[0]
    nrows = get_nrows(Np)
    
    sj = make_J_edges1(N,nrows)


    (data,irows,icols) = ([],[],[])
    
    Ne = 0
    for ne in range(N):
        (d,ir,ic,nrows) = make_edge_row_c(mesh,ne,Np,Ndeg)
        #(d,ir,ic,nrows)=make_edge_row(mesh,ne,Np,Ndeg);
        ir+=Ne
        Ne+=nrows
        data+=[d]
        irows+=[ir]
        icols+=[ic]
        pass
    irows = np.concatenate(irows,None)
    icols = np.concatenate(icols,None)
    data = np.concatenate(data,None)
    sm = sp.coo_matrix((data,(irows,icols)),dtype=np.float64) 
    #sj=sj*sm;
    return (sm,sj)
    pass

def make_RJ(mesh,Np=1,Ndeg=-1):
    (sm,sj) = make_H_edges(mesh,Np,Ndeg)
    smt = sm.T
    (sm,sj,smt) = [v.tocsc() for v in (sm,sj,smt) ]
    R = smt * sm
    J = smt * sj
    return (R,J)

def make_reg_G(Ndeg,rfu):
    w = [rfu(k) for  k in range(Ndeg) ]


def make_RJr(mesh,Np=1,Ndeg=1,rfu=math.factorial):

    (sm,sj) = make_H_edges(mesh,Np,Ndeg)
    smt = sm.T
    (sm,sj,smt) = [v.tocsc() for v in (sm,sj,smt) ]
    R = smt * sm
    J = smt * sj
    return (R,J)


def generate_HC(mesh,je,Np=1,Ndeg=-1):
    
    __magic__ = (1 + 1j) / np.sqrt(2)
    mesh = to_jso(mesh)
    (R,J) = make_RJ(mesh,Np,Ndeg)
    b = J * je
    hc = la.spsolve(R,b)
    Nr = mesh.rects.shape[0]    
    hc = hc.reshape((Nr,2,Ndeg + 1))
    c = hc[:,0,:] + 1j * hc[:,1,:]
    c = c.reshape((Nr,Ndeg + 1))
    #c*=__magic__;
    return c



def in_rect(z,rz):
    xmin,xmax = rz[0].real,rz[3].real
    ymin,ymax = rz[0].imag,rz[3].imag
    x,y = z.real,z.imag
    f = (xmin <= x) and (x < xmax) and (ymin <= y) and (y < ymax) 
    return f

def get_irect(z,rzs,c):
    k = 0
    for rs in rzs:
        if in_rect(z,rs):
            return k
        k+=1
    return -1

def get_value0(z,rzs,c):
    
    k = get_irect(z,rzs,c)
    if k >= 0:
        return polyval(z,c[k])
    return 0

def get_value(zs,rzs,c):
    t = type(np.array(0))
    if type(zs) == t:
        shape = zs.shape
        ps = zs.flatten()
        rs = np.empty_like(ps)
        k = 0
        for z in ps:
            rs[k] = get_value0(z,rzs,c)
            k+=1
        return rs.reshape(shape)
    else:
        return get_value0(zs,rzs,c)

def fun_app(zs,gf):

    t = type(np.array(0))

    if type(zs) != t:
        return gf(zs)
    else:
        shape = zs.shape
        ps = zs.flatten()
        rs = np.empty_like(ps)
        k = 0
        for z in ps:
            rs[k] = gf(z)
            k+=1
        return rs.reshape(shape)

        

def make_HCz(mesh,je,Np=1,Ndeg=-1):
    mesh = to_jso(mesh)
    c = generate_HC(mesh,je,Np,Ndeg)
    rs = mesh.rects
    rzs = mesh.zs[rs]
    return (c,rzs)

def HCz_functor(mesh,je,Np=1,Ndeg=-1):
    (c,rzs) = make_HCz(mesh,je,Np,Ndeg)
    return lambda z: get_value(z,rzs,c)



def HCz_functor2(mesh,je,Np=1,Ndeg=-1):
    from MID.hall.rectmesh import make_interp_index
    (c,rzs) = make_HCz(mesh,je,Np,Ndeg)
    zrc = mesh.zrc
    r_index = make_interp_index(mesh)

    def getv(z):
        k = r_index(z)
        if k >= 0:
            return polyval(z - zrc[k],c[k])
        else:
            return 0
    
    return  lambda z: fun_app(z,getv)



def getHxy(mesh,cr,z):
    pass


if __name__ == '__main__':
    
    pass