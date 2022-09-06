
# from __future__ import division
# try:
#     from runtime import *
# except ImportError:
#     from smop.runtime import *

from .pade_laplace_LU import pade_laplace_LU_

from .emtClass import Dat

def mid_utils_(*args,**kwargs):
    class MidU:
        pass

    midu = MidU()

    # midu.evalif=evalif_
    # midu.morph_grid_D1=morph_D1_
    # midu.c0=c0_
    # midu.gmid=gmid_
    # midu.gmidc=gmidc_
    midu.rescale_mesh=rescale_mesh_
    midu.make_scheme=make_scheme_
    # midu.mesh_gen=mesh_gen_
    # midu.bsd=bsd_
    # midu.expandc1=expandc1_
    # midu.expandc2=expandc2_
    # midu.expand1=expand1_
    # midu.expand2=expand2_
    # midu.crop_bound=crop_bound_
    midu.mesh_mask=mesh_mask_
    midu.set_region_params=set_region_params_
    # midu.initPL=initPL_
    # midu.mask_rect=mask_rect_
    # midu.mask_rects=mask_rects_
    # midu.mirrz2=mirrz2_
    # midu.geom_param=geom_param_
    # midu.renorm_f=renorm_f_
    midu.opt_def=opt_def_
    # midu.derr2m=derr2m_
    # midu.rerr=rerr_
    # midu.rerrn=rerrn_
    midu.pade_laplace=pade_laplace_
    midu.options_default=options_default_
    return midu

def derr2m_(x=None,y=None,epsa=None,*args,**kwargs):
    varargin = cellarray(args)
    nargin = 3-[x,y,epsa].count(None)+len(args)

    if not epsa:
        epsa=np.spacing(1)
    axx=abs_(y) + abs_(x)
    axx[axx < epsa]=1
    xy=y - x
    xy[abs_(xy) < epsa]=0
    err=xy / axx
    return err

def renorm_f_(f=None,t=None,r=None,*args,**kwargs):
    varargin = cellarray(args)
    nargin = 3-[f,t,r].count(None)+len(args)

    fr=f[:]
    if nargin < 3:
        fr=fr / np.mean(abs_(fr))
    else:
        ii=arange_(1,numel_(t))
        ii=ii[(t >= r[1]) and (t <= r[2])]
        fr=fr / np.mean(abs_(fr[ii]))
    return fr

def pade_laplace_(*args,**kwargs):
    pls=pade_laplace_LU_()
    return pls

def opt_def_(m=None,md=None,*args,**kwargs):
    pls=pade_laplace_LU_()
    m=pls.opt_def(m,md)
    return m

def mask_rect_(sy=None,r=None,*args,**kwargs):

    r=r.ravel()
    xc=sy.Xc
    yc=sy.Yc
    b=(r[0] < xc) & (xc < r[2]) & (r[1] < yc) & (yc < r[3])
    b=b.astype(int)
    bi=np.logical_not(b).astype(int)
    return b#,bi

def mask_rects_(sy=None,rs=None,*args,**kwargs):

    rs=evalif_(rs)
    nr=1
    bs=np.zeros_like(sy.Xc)
    for k in np.arange(nr):
        r=c_a_(rs)
        b=mask_rect_(sy,r)
        bs = np.logical_or(bs, b).astype(int)
    bsi=np.logical_not(bs).astype(int)
    return bs,bsi

def mesh_dflt_(sy=None,*args,**kwargs):
    varargin = cellarray(args)
    nargin = 1-[sy].count(None)+len(args)

    n,m=size_(sy.X,nargout=2)
    t=copy_(sy)
    mu0=np.ones([n - 1,m - 1])
    zz=np.zeros([n - 1,m - 1])
    sigma0=copy_(zz)

    class M0:
        pass

    m0 = M0()

    m0.name=char('')
    m0.mu=1
    m0.sigma=0
    m0.rects=[]
    t.mu=m0.mu
    t.sigma=m0.sigma
    t.meshs=[]
    t.disabled=0
    sy=opt_def_(sy,t)
    sy=c00_(sy)
    return sy,m0,mu0,sigma0,zz

def set_region_params_(name=None,sigma=0,mu=1,rects=None):
    #varargin = cellarray(args)
    #nargin = 4-[name,sigma,mu,rects].count(None)+len(args)

    class Region:
        pass

    reg = {}

    reg["name"] = name
    reg["sigma"] = sigma

    if mu <=0:
        raise Exception("Error: Mu must bu great zero")

    reg["mu"] = mu
    reg["rects"] = rects
    reg["disabled"] = 0
    return reg

def mesh_mask_(sy=None,meshes=None,*args,**kwargs):
    varargin = cellarray(args)
    nargin = 2-[sy,meshes].count(None)+len(args)

    #if not isfield_(sy,char('X')) or not isfield_(sy,char('Y')):
    if not hasattr(sy, 'X') or not hasattr(sy, 'Y'):
        sy.X,sy.Y=np.meshgrid(sy.rl,sy.zl)
    if nargin > 1:
        sy.meshes=meshes
    sy,m0,mu_full,sigma_full,zzm=mesh_dflt_(sy,nargout=5)
    mss=sy.meshes
    nm=numel_(mss)
    if nm:
        for h in np.arange(nm):
            ms=opt_def_(mss[h],m0)
            if not ms.disabled:
                b,bi=mask_rects_(sy,ms.rects)
                mu_full = mu_full*bi + ms.mu * b
                sigma_full=sigma_full*bi + ms.sigma * b
                zzm=zzm + b
                ms.mask=b
            mss[h]=ms
        sy.overlap_map=zzm
        sy.meshes=mss
        sy.mu=mu_full
        sy.sigma=sigma_full
    return sy

def c00_(sy=None,*args,**kwargs):

    n=len(sy.zl)
    n1= 1 #len(np.array(sy.mu))
    f=(n1 >= (n + 2))
    sy=c0_(sy,f)
    return sy

def initPL_(sy=None,t=None,pade=None,*args,**kwargs):
    varargin = cellarray(args)
    nargin = 3-[sy,t,pade].count(None)+len(args)

    if not pade:
        pade=matlabarray([2,4])
    pls=pade_laplace_()
    pd=pls.get_data(pade[1],pade[2])
    sy=c00_(sy)
    #tic
    #fprintf_(char('make scheme...'))
    #toc
    sn=midu.rescale_mesh(sy)
    scheme,sn=midu.make_scheme(sn,nargout=2)
    #fprintf_(char('make AC...'))
    #toc
    A=scheme.make_AC()
    #fprintf_(char('end\\n'))
    #toc
    s.scheme=scheme
    s.sn=sn
    s.to_mesh=scheme.imesh.to_mesh
    s.J=A.J * 1e+14
    return s

def gmid_(d=None,rc=None,*args,**kwargs):
    varargin = cellarray(args)
    nargin = 2-[d,rc].count(None)+len(args)

    d2rc=d / (2 * rc)
    g=d / np.log((1 + d2rc) / (1 - d2rc))
    return g

def gmidc_(d=None,rc=None,*args,**kwargs):
    varargin = cellarray(args)
    nargin = 2-[d,rc].count(None)+len(args)

    g=gmid_(d,rc)
    gc=g - rc
    return gc,g

def evalif_(v=None,*args,**kwargs):
    varargin = cellarray(args)
    nargin = 1-[v].count(None)+len(args)

    if ischar_(v):
        v=np.eval(v)
    return v

def split_lr_(x=None,*args,**kwargs):

    eee = len(x)

    xl=x[0:eee - 1]
    xr=x[1:eee]
    return xl,xr,x

def crop_bound_(sy=None,rect=None,*args,**kwargs):
    varargin = cellarray(args)
    nargin = 2-[sy,rect].count(None)+len(args)

    lr=rect[1,1]
    rr=rect[2,1]
    lz=rect[1,2]
    rz=rect[2,2]
    return sy
def expand1_(xy=None,*args,**kwargs):
    varargin = cellarray(args)
    nargin = 1-[xy].count(None)+len(args)

    n,m=size_(xy,nargout=2)
    xy2=zeros_(2 * n - 1,m)
    for k in arange_(1,n).reshape(-1):
        xy2[2 * k - 1,:]=xy[k,:]
    xyc=zeros_(n - 1,m)
    for k in arange_(1,m).reshape(-1):
        xyc[:,k]=make_centerD1_(xy[:,k])
    for k in arange_(1,n - 1).reshape(-1):
        xy2[2 * k,:]=xyc[k,:]
    return xy2
def expandc1_(b=None,*args,**kwargs):
    varargin = cellarray(args)
    nargin = 1-[b].count(None)+len(args)

    n,m=size_(b,nargout=2)
    b2=zeros_(2 * n,m)
    for k in arange_(1,n).reshape(-1):
        b2[2 * k - 1,:]=b[k,:]
        b2[2 * k,:]=b[k,:]
    return b2
def expandc2_(b=None,*args,**kwargs):
    varargin = cellarray(args)
    nargin = 1-[b].count(None)+len(args)

    b=(expandc1_(b.T)).T
    return b
def expand2_(xy=None,*args,**kwargs):
    varargin = cellarray(args)
    nargin = 1-[xy].count(None)+len(args)

    xy=(expand1_(xy.T)).T
    return xy
def bsd1_(sy=None,*args,**kwargs):
    varargin = cellarray(args)
    nargin = 1-[sy].count(None)+len(args)

    sy.zl=expand1_(sy.zl)
    sy.X=expand1_(sy.X)
    sy.Y=expand1_(sy.Y)
    sy.mu=expandc1_(sy.mu)
    sy.sigma=expandc1_(sy.sigma)
    sy.rj=expandc1_(sy.rj)
    sy.gj=expandc1_(sy.gj)
    return sy
def bsd2_(sy=None,*args,**kwargs):
    varargin = cellarray(args)
    nargin = 1-[sy].count(None)+len(args)

    sy.rl=expand1_(sy.rl)
    sy.X=expand2_(sy.X)
    sy.Y=expand2_(sy.Y)
    sy.mu=expandc2_(sy.mu)
    sy.sigma=expandc2_(sy.sigma)
    sy.rj=expandc2_(sy.rj)
    sy.gj=expandc2_(sy.gj)
    return sy

def bsd_(sy=None,fxy=None,*args,**kwargs):
    varargin = cellarray(args)
    nargin = 2-[sy,fxy].count(None)+len(args)

    if exist_(char('fxy')) != 1:
        fxy=3
    if bitand_(fxy,1):
        sy=bsd1_(sy)
    if bitand_(fxy,2):
        sy=bsd2_(sy)
    sy=c0_(sy)
    return sy

def mesh_gen_(sy=None,*args,**kwargs):
    varargin = cellarray(args)
    nargin = 1-[sy].count(None)+len(args)
    return sy

def make_fluxes_raw_(scheme=None,*args,**kwargs):
    varargin = cellarray(args)
    nargin = 1-[scheme].count(None)+len(args)

    sf.NN=numel_(scheme.cs)
    at=scheme.imesh.at
    bc=scheme.geom.bc
    h11r=scheme.h11r
    h12r=scheme.h12r
    nn,mm=size_(h11r,nargout=2)
    znm=zeros_([nn,mm])
    fl=copy_(znm)
    i1=copy_(znm)
    i2=copy_(znm)
    s11=copy_(znm)
    s12=copy_(znm)
    s21=copy_(znm)
    s22=copy_(znm)
    znm=matlabarray([])
    testg=2
    for n in arange_(1,nn).reshape(-1):
        for m in arange_(1,mm).reshape(-1):
            h11=h11r[n,m]
            h12=h12r[n,m]
            i1[n,m]=at[n,m]
            i2[n,m]=at[n,m + 1]
            s11[n,m]=h11
            s12[n,m]=h12
            s21[n,m]=- h11
            s22[n,m]=- h12
            fl[n,m]=1
        if (bc.r0 == 2) or (bc.ri == 2):
            m=1
            h11=h11r[n,m]
            h12=h12r[n,m]
            i1[n,m]=at[n,m]
            i2[n,m]=at[n,mm]
            s11[n,m]=h11
            s12[n,m]=h12
            s21[n,m]=- h11
            s22[n,m]=- h12
            fl[n,m]=1
        else:
            if (bc.r0):
                m=1
                h11=h11r[n,m]
                h12=h12r[n,m]
                i1[n,m]=at[n,m]
                i2[n,m]=at[n,m + 1]
                s11[n,m]=h11
                s12[n,m]=0
                s21[n,m]=- h11
                s22[n,m]=- h12
                fl[n,m]=1
            if (bc.ri):
                m=copy_(mm)
                h11=h11r[n,m]
                h12=h12r[n,m]
                i1[n,m]=at[n,m]
                i2[n,m]=at[n,m + 1]
                s11[n,m]=h11
                s12[n,m]=h12
                s21[n,m]=0
                s22[n,m]=- h12
                fl[n,m]=1
    sf.fl=fl[:]
    sf.i1=i1[:]
    sf.i2=i2[:]
    sf.s11=s11[:]
    sf.s12=s12[:]
    sf.s21=s21[:]
    sf.s22=s22[:]
    h11z=scheme.h11z
    h12z=scheme.h12z
    nn,mm=size_(h11z,nargout=2)
    znm=zeros_([nn,mm])
    fl=copy_(znm)
    i1=copy_(znm)
    i2=copy_(znm)
    s11=copy_(znm)
    s12=copy_(znm)
    s21=copy_(znm)
    s22=copy_(znm)
    znm=matlabarray([])
    for m in arange_(1,mm).reshape(-1):
        for n in arange_(1,nn).reshape(-1):
            h11=h11z[n,m]
            h12=h12z[n,m]
            i1[n,m]=at[n,m]
            i2[n,m]=at[n + 1,m]
            s11[n,m]=h11
            s12[n,m]=h12
            s21[n,m]=- h11
            s22[n,m]=- h12
            fl[n,m]=1
        if (bc.zt == 2) or (bc.zb == 2):
            n=1
            h11=h11z[n,m]
            h12=h12z[n,m]
            i1[n,m]=at[n,m]
            i2[n,m]=at[nn,m]
            s11[n,m]=h11
            s12[n,m]=h12
            s21[n,m]=- h11
            s22[n,m]=- h12
            fl[n,m]=1
        else:
            if (bc.zt):
                n=1
                h11=h11z[n,m]
                h12=h12z[n,m]
                i1[n,m]=at[n,m]
                i2[n,m]=at[n + 1,m]
                s11[n,m]=h11
                s12[n,m]=0
                s21[n,m]=- h11
                s22[n,m]=- h12
                fl[n,m]=1
            if (bc.zb):
                n=copy_(nn)
                h11=h11z[n,m]
                h12=h12z[n,m]
                i1[n,m]=at[n,m]
                i2[n,m]=at[n + 1,m]
                s11[n,m]=h11
                s12[n,m]=h12
                s21[n,m]=0
                s22[n,m]=- h12
                fl[n,m]=1
    sf.fl=[[sf.fl],[fl[:]]]
    sf.i1=[[sf.i1],[i1[:]]]
    sf.i2=[[sf.i2],[i2[:]]]
    sf.s11=[[sf.s11],[s11[:]]]
    sf.s12=[[sf.s12],[s12[:]]]
    sf.s21=[[sf.s21],[s21[:]]]
    sf.s22=[[sf.s22],[s22[:]]]
    global rsf
    rsf=copy_(sf)
    return sf


def make_fluxes_(scheme=None,*args,**kwargs):
    varargin = cellarray(args)
    nargin = 1-[scheme].count(None)+len(args)

    at=scheme.imesh.at
    bc=scheme.geom.bc
    h11r=scheme.h11r
    h12r=scheme.h12r
    nn,mm=size_(h11r,nargout=2)
    fluxesr=cell_([nn,mm])
    f=matlabarray([])
    testg=2
    for n in arange_(1,nn).reshape(-1):
        for m in arange_(1,mm).reshape(-1):
            h11=h11r[n,m]
            h12=h12r[n,m]
            f.nn=[at[n,m],at[n,m + 1]]
            f.h2x2=[[h11,h12],[- h11,- h12]]
            fluxesr[n,m]=f
        if (bc.r0 == 2) or (bc.ri == 2):
            m=1
            h11=h11r[n,m]
            h12=h12r[n,m]
            f.nn=[at[n,m],at[n,mm]]
            f.h2x2=[[h11,h12],[- h11,- h12]]
            fluxesr[n,m]=f
        else:
            if (bc.r0):
                m=1
                h11=h11r[n,m]
                h12=h12r[n,m]
                f.nn=[at[n,m],at[n,m + 1]]
                f.h2x2=[[h11,0],[- h11,- h12]]
                fluxesr[n,m]=f
            if (bc.ri):
                m=copy_(mm)
                h11=h11r[n,m]
                h12=h12r[n,m]
                f.nn=[at[n,m],at[n,m + 1]]
                f.h2x2=[[h11,h12],[0,- h12]]
                fluxesr[n,m]=f
    h11z=scheme.h11z
    h12z=scheme.h12z
    nn,mm=size_(h11z,nargout=2)
    fluxesz=cell_([nn,mm])
    f=matlabarray([])
    for m in arange_(1,mm).reshape(-1):
        for n in arange_(1,nn).reshape(-1):
            h11=h11z[n,m]
            h12=h12z[n,m]
            f.nn=[at[n,m],at[n + 1,m]]
            f.h2x2=[[h11,h12],[- h11,- h12]]
            fluxesz[n,m]=f
        if (bc.zt == 2) or (bc.zb == 2):
            n=1
            h11=h11z[n,m]
            h12=h12z[n,m]
            f.nn=[at[n,m],at[nn,m]]
            f.h2x2=[[h11,h12],[- h11,- h12]]
            fluxesz[n,m]=f
        else:
            if (bc.zt):
                n=1
                h11=h11z[n,m]
                h12=h12z[n,m]
                f.nn=[at[n,m],at[n + 1,m]]
                f.h2x2=[[h11,0],[- h11,- h12]]
                fluxesz[n,m]=f
            if (bc.zb):
                n=copy_(nn)
                h11=h11z[n,m]
                h12=h12z[n,m]
                f.nn=[at[n,m],at[n + 1,m]]
                f.h2x2=[[h11,h12],[0,- h12]]
                fluxesz[n,m]=f
    fluxes=matlabarray([[fluxesr[:]],[fluxesz[:]]])
    return fluxes,fluxesr,fluxesz
def SLD_LU_cache_(scheme=None,AC=None,dt=None,padenm=None,*args,**kwargs):
    varargin = cellarray(args)
    nargin = 4-[scheme,AC,dt,padenm].count(None)+len(args)

    pls=scheme.pd
    pd=pls.get_data(padenm)
    lu=pls.LU_cache(pd,AC,dt,1,scheme.geom.fnocycref)
    return lu
def SLD_LU_cache_hyp_(scheme=None,AC=None,dt=None,padenm=None,fhyp=None,*args,**kwargs):
    varargin = cellarray(args)
    nargin = 5-[scheme,AC,dt,padenm,fhyp].count(None)+len(args)

    pls=scheme.pd
    pd=pls.get_data(padenm)
    pd.fhyp=fhyp
    lu=pls.LU_cache(pd,AC,dt,1,scheme.geom.fnocycref)
    return lu
def make_SLD_(scheme=None,fhyp=None,fcompact=None,*args,**kwargs):
    varargin = cellarray(args)
    nargin = 3-[scheme,fhyp,fcompact].count(None)+len(args)

    if exist_(char('fhyp')) != 1:
        fhyp=0
    if exist_(char('fcompact')) != 1:
        fcompact=0
    pd=scheme.pd
    AC.C,AC.J,AC.V=make_CJ_(scheme,nargout=3)
    NN=numel_(scheme.cs)
    fraw=scheme.geom.fluxes_raw
    if fraw:
        sf=make_fluxes_raw_(scheme)
        AC.A=pd.set_2x2_fluxes_raw(sf)
    else:
        fluxes=make_fluxes_(scheme)
        AC.A=pd.set_2x2_fluxes(fluxes,NN)
    if not scheme.geom.fnocycref:
        AC.jd=lambda f: get_jd_(scheme,AC,f)
        AC.jd1=lambda f1: get_jd1_(scheme,f1)
        AC.jmesh=lambda f,J: get_jmesh_(scheme,AC,f,J)
        AC.jmesh1=lambda f1: get_jmesh1_(scheme,f1)
    if (not fcompact) and (not fraw):
        AC.fluxes=fluxes
    if (fhyp):
        AC.eps=scheme.sn.epsc * AC.V
    AC.SLD_LU_cache=SLD_LU_cache
    AC.SLD_LU_cache_hyp=SLD_LU_cache_hyp
    if not scheme.geom.fnocycref:
        AC.LU_cache=lambda dt,padenm: SLD_LU_cache_(scheme,AC,dt,padenm)
        AC.LU_cache_hyp=lambda dt,padenm,hyp: SLD_LU_cache_hyp_(scheme,AC,dt,padenm,hyp)
    return AC
def get_jmesh1_(scheme=None,f1=None,*args,**kwargs):
    varargin = cellarray(args)
    nargin = 2-[scheme,f1].count(None)+len(args)

    f1=scheme.imesh.to_mesh(f1)
    jd=f1.dot(scheme.geom.Xc).dot(scheme.geom.sigma)
    return jd
def get_jmesh_(scheme=None,AC=None,f=None,J=None,*args,**kwargs):
    varargin = cellarray(args)
    nargin = 4-[scheme,AC,f,J].count(None)+len(args)

    if exist_(char('J')) != 1:
        J=[]
    jd=AC.A * (f[:])
    if not isempty_(J):
        jd=jd + J
    jd=scheme.imesh.to_mesh(jd) / scheme.vs
    return jd
def get_jd1_(scheme=None,f1=None,*args,**kwargs):
    varargin = cellarray(args)
    nargin = 2-[scheme,f1].count(None)+len(args)

    imesh=scheme.imesh
    rj=scheme.geom.rj(arange_())
    cs=scheme.cs(arange_())
    f1=f1[rj > 0]
    cs=cs[rj > 0]
    jd=sum_(cs.dot(f1))
    return jd
def get_jd_(scheme=None,AC=None,f=None,*args,**kwargs):
    varargin = cellarray(args)
    nargin = 3-[scheme,AC,f].count(None)+len(args)

    jd=AC.A * (f[:])
    rj=scheme.geom.rj(arange_())
    jd=jd[rj > 0]
    jd=sum_(jd)
    return jd
def l_r_(AC=None,f=None,J=None,*args,**kwargs):
    varargin = cellarray(args)
    nargin = 3-[AC,f,J].count(None)+len(args)

    err=- AC.A * (f[:])
    if (exist_(char('J')) == 1) and (not isempty_(J)):
        err=err + J[:]
    return err
def l_l_(AC=None,f1=None,*args,**kwargs):
    varargin = cellarray(args)
    nargin = 2-[AC,f1].count(None)+len(args)

    err=AC.C * (f1[:])
    return err
def rerrn_(AC=None,f=None,f1=None,J=None,*args,**kwargs):
    varargin = cellarray(args)
    nargin = 4-[AC,f,f1,J].count(None)+len(args)

    d=6
    return err
def rerr_(AC=None,f=None,f1=None,J=None,*args,**kwargs):
    varargin = cellarray(args)
    nargin = 4-[AC,f,f1,J].count(None)+len(args)

    err=AC.A * (f[:]) + AC.C * (f1[:])
    if (exist_(char('J')) == 1) and not isempty_(J):
        err=err + J[:]
    err=abs_(err)
    return err
def make_CJ_(scheme=None,*args,**kwargs):
    varargin = cellarray(args)
    nargin = 1-[scheme].count(None)+len(args)

    pd=scheme.pd
    cs=scheme.cs
    vs=scheme.vs
    js=scheme.jmask
    nn,mm=size_(cs,nargout=2)
    NM=nn * mm
    imesh=pd.index2D([nn,mm])
    at=imesh.at
    cd=zeros_([NM,1])
    vd=copy_(cd)
    J=copy_(cd)
    for n in arange_(1,nn).reshape(-1):
        for m in arange_(1,mm).reshape(-1):
            nm=at[n,m]
            cd[nm]=cs[n,m]
            vd[nm]=vs[n,m]
            J[nm]=js[n,m]
    NNN=arange_(1,NM)
    C=sparse_(NNN,NNN,cd)
    V=sparse_(NNN,NNN,vd)
    return C,J,V,cs,js,vs

def geom_param_(x=None,*args,**kwargs):

    eee = len(x)

    xl=x[0:eee-1]
    xr=x[1:eee]
    xc=(xl + xr) / 2
    dx=(xr - xl)
    xc0,xc1,Null =split_lr_(xc,nargout=2)
    dx0,dx1,Null =split_lr_(dx,nargout=2)

    return dx,xc,x,dx0,xc0,dx1,xc1

def make_scheme_(sn=None,*args,**kwargs):
    varargin = cellarray(args)
    nargin = 1-[sn].count(None)+len(args)

    dr,rc,r,dr0,rc0,dr1,rc1=geom_param_(sn.rl,nargout=7)
    dz,zc,z,dz0,zc0,dz1,zc1=geom_param_(sn.zl,nargout=7)
    grc,gr=gmidc_(dr,rc,nargout=2)
    if sn.geom.fmfv == 0:
        gr=copy_(rc)
        grc=0 * rc
    gr0,gr1,x=split_lr_(gr,nargout=2)
    grc0,grc1,x=split_lr_(grc,nargout=2)
    lambda_r0=(1 / 2) * (dr0 + dr1) - (grc0 + grc1)
    lambda_z0=(1 / 2) * (dz0 + dz1)

    scheme = Dat()

    scheme.lambda_r0=lambda_r0
    scheme.lambda_z0=lambda_z0
    nz=numel_(z)
    nr=numel_(r)
    sc=sn.geom.normalizator
    zz=np.zeros([nz - 1,nr - 2])
    lambda_r=np.zeros([nz - 1,nr - 2])
    h11r=np.zeros([nz - 1,nr - 2])
    h12r=np.zeros([nz - 1,nr - 2])
    #zz=matlabarray([])
    pow2=2 ** (sn.geom.distortion / 2)
    scr=sc * pow2 / (1 + pow2)
    scz=sc * 1 / (1 + pow2)
    scr=sc * pow2
    scz=sc * 1 / (pow2)
    for k in np.arange(nz - 1):
        muk0,muk1,muk=split_lr_(sn.mu[k],nargout=3)
        muk0=muk0.T
        muk1=muk1.T
        l=(1 / 2) * (dr0.dot(muk0) + dr1.dot(muk1)) - (grc0.dot(muk0) + grc1.dot(muk1))
        lambda_r[k,:]=l
        q=scr * (4 * dz[k]) / l
        h11r[k]= -q*gr0
        h12r[k]=q*gr1
    zz=np.zeros([nz - 2,nr - 1])
    lambda_z=copy_(zz)
    h11z=copy_(zz)
    h12z=copy_(zz)
    #zz=matlabarray([])
    for j in np.arange(nr - 1):
        muj0,muj1,muj=split_lr_(sn.mu[:,j],nargout=3)
        l=(1 / 2) * (dz0.dot(muj0) + dz1.dot(muj1))
        lambda_z[:,j]=l
        q=scz * (dr[j]) / l
        h11z[:,j]=- q
        h12z[:,j]=q
    sigma=sn.sigma
    vs=np.kron(np.reshape(dr, (-1, 1)),dz).transpose()
    cs=vs*sigma * sn.beta * sc
    scheme.cs=cs
    scheme.vs=vs
    jmask=sn.geom.gj
    scheme.jmask=jmask*vs * sc
    NN=numel_(cs)
    scheme.NN=NN
    scheme.imesh=sn.pd.index2D(size_(cs))
    scheme.zero=zeros_(NN,1)
    scheme.lambda_r=lambda_r
    scheme.lambda_z=lambda_z
    scheme.h11r=h11r
    scheme.h12r=h12r
    scheme.h11z=h11z
    scheme.h12z=h12z
    scheme.pd=sn.pd
    scheme.geom=sn.geom
    scheme.sn=sn
    if not scheme.geom.fnocycref:
        scheme.make_fluxes=lambda : make_fluxes_(scheme)
        scheme.make_fluxes_raw=lambda : make_fluxes_raw_(scheme)
        scheme.make_CJ=lambda : make_CJ_(scheme)
        scheme.make_AC=lambda fhyp,f: make_SLD_(scheme,fhyp,f)
        scheme.make_SLD=lambda fhyp,f: make_SLD_(scheme,fhyp,f)
        scheme.get_crop_c=lambda r: get_crop_(sn.geom,r,1)
        scheme.get_crop=lambda r: get_crop_(sn.geom,r,0)
        scheme.get_fHzEx=lambda : get_fHz_(sn.geom)
        scheme.get_fHz=lambda : get_fHz0_(sn.geom)
        scheme.get_fHzdn=lambda d,n: get_fHzdn_(sn.geom,d,n,1)
        scheme.get_fBzdn=lambda d,n: get_fHzdn_(sn.geom,d,n,0)
    scheme.make_SLD_scheme=make_SLD
    scheme.make_fluxes_scheme=make_fluxes
    scheme.make_CJ_scheme=make_CJ
    scheme.make_AC_scheme=make_SLD
    scheme.make_SLD_scheme=make_SLD
    return scheme,sn

def get_fHzdn_(sy=None,d=None,n=None,fhm=None,*args,**kwargs):
    varargin = cellarray(args)
    nargin = 4-[sy,d,n,fhm].count(None)+len(args)

    r=sy.Xc ** n
    if fhm:
        r=r / sy.mu
    if d > 0:
        fHz=lambda x: abs_((r).dot(x)) ** d
    else:
        fHz=lambda x: log_(abs_((r).dot(x)))
    return fHz
def get_fHz_(sy=None,*args,**kwargs):
    varargin = cellarray(args)
    nargin = 1-[sy].count(None)+len(args)

    rm=sy.Xc
    fHz=lambda x,d,n: abs_((rm ** n).dot(x)) ** d
    return fHz
def get_fHz0_(sy=None,*args,**kwargs):
    varargin = cellarray(args)
    nargin = 1-[sy].count(None)+len(args)

    rm=sy.Xc
    fHz=lambda x,d,n: abs_((rm ** 1).dot(x)) ** 0.5
    return fHz
def get_crop_(sy=None,rect=None,fc=None,*args,**kwargs):
    varargin = cellarray(args)
    nargin = 3-[sy,rect,fc].count(None)+len(args)

    xl=rect[1]
    xr=rect[2]
    yl=rect[3]
    yr=rect[4]
    if fc:
        x=sy.rlc
        y=sy.zlc
    else:
        x=sy.rl
        y=sy.zl
    n=numel_(x)
    nn=arange_(1,n)
    m=numel_(y)
    mm=arange_(1,m)
    nncr=nn[(xl <= x) and (x < xr)]
    mmcr=mm[(yl <= y) and (y < yr)]
    crop=lambda v: v[mmcr,nncr]
    return crop

def options_default_(*args,**kwargs):
    varargin = cellarray(args)
    nargin = 0-[].count(None)+len(args)

    t = Dat()

    t.L0=5000
    t.To=1
    t.inflat_z=np.sqrt(1)
    t.inflat_r=t.inflat_z ** 2
    t.sigma_titan_0=1.2 * 10 ** 6
    t.sigma_titan_0=1.05 * 10 ** 6
    cl=299792458
    t.clight=cl
    cgc=8.98755 * 10 ** 9
    t.cgc_sigm=cgc
    mu00=1.25663706144
    t.beta_titan_0=1.3195
    betit=4 * np.pi * (1050000.0 * cgc) / cl ** 2
    t.beta_titan_0=betit
    t.fBC=[[0,0],[1,0]]
    t.normalizator=1
    t.fmfv=1
    t.distortion=0
    t.fpoolcache=0
    t.fpool=0
    t.rj=0
    t.gj=0

    t.bc = Dat()

    t.bc.zt = 1
    t.bc.zb = 1
    t.bc.r0 = 1
    t.bc.ri = 1

    t.scaler=1
    t.fnocycref=1
    t.fluxes_raw=1
    return t

def rescale_mesh_(sy=None,sy_dflt=None,*args,**kwargs):
    varargin = cellarray(args)
    nargin = 2-[sy,sy_dflt].count(None)+len(args)

    if not (hasattr(sy,'X') and hasattr(sy,'Y')):
        sy.X,sy.Y=np.meshgrid(sy.rl,sy.zl)
    pd=pade_laplace_()
    opt_def=pd.opt_def
    if nargin < 2:
        sy_dflt=options_default_()

    sy=opt_def(sy,sy_dflt)
    meter=1000
    L0=sy.L0
    L0meter=L0 / meter
    dlf=(L0meter ** 2) / sy.To
    beta=sy.scaler * sy.beta_titan_0 * dlf
    epsc=(L0meter ** 3) / sy.To / (sy.clight ** 2)
    epsc=dlf / (sy.clight ** 2) / (4 * np.pi)

    sn = Dat()

    sn.epsc=epsc
    sn.beta=beta
    Lr=max_(sy.rl) / sy.inflat_r
    Lz=max_(sy.zl) / sy.inflat_z
    rl=sy.rl / L0
    rl=rl ** 2
    zl=sy.zl / L0
    sn.rl=rl
    sn.zl=zl
    sn.mu=sy.mu
    sy.volm=(L0meter ** 3)
    sn.sigma=sy.sigma / sy.sigma_titan_0
    tt,sy.rlc,x,dx0,xc0,dx1,xc1=geom_param_(sy.rl,nargout=2)
    tt,sy.zlc,x,dx0,xc0,dx1,xc1=geom_param_(sy.zl,nargout=2)
    sy.Xc,sy.Yc=np.meshgrid(sy.rlc,sy.zlc)
    sn.geom=sy
    sn.pd=pd
    return sn

def mirrz2_(v=None,*args,**kwargs):
    varargin = cellarray(args)
    nargin = 1-[v].count(None)+len(args)

    n,m=size_(v,nargout=2)
    vv=zeros_(2 * n,m)
    vv[1:n,:]=flipud_(v)
    vv[n + 1:2 * n,:]=v
    return vv

def make_centerD1_(r=None,*args,**kwargs):
    n=numel_(r)
    rc=(r[0:n - 1] + r[1:n]) / 2
    return rc

def c0_(sy=None,ff=None,*args,**kwargs):

    if not 'ff':
        ff=0
    if ff:
        sy.mu=mirrz2_(sy.mu)
        sy.sigma=mirrz2_(sy.sigma)
        sy.rj=mirrz2_(sy.rj)
        sy.gj=mirrz2_(sy.gj)
    sy.rlc=make_centerD1_(sy.rl)
    sy.zlc=make_centerD1_(sy.zl)
    sy.Xc,sy.Yc=np.meshgrid(sy.rlc,sy.zlc)
    return sy

def c_a_(c=None,*args,**kwargs):
    c=evalif_(c)
    #if iscell_(c):
    #    c=cell2mat_(c)
    return c

def c_a0_(c=None,*args,**kwargs):
    c=evalif_(c)
    #if iscell_(c):
    #    c=cell2mat_(c)
    return c

def morph_D1_(segs=None,*args,**kwargs):
    varargin = cellarray(args)
    nargin = 1-[segs].count(None)+len(args)

    if (not ischar_(segs)) and (not iscell_(segs)):
        Gm=copy_(segs)
        return Gm,Gm1
    segs=evalif_(segs)
    segs=segs[:]
    ns=size_(segs,1)
    Gm=c_a0_(segs[1])
    Gm=Gm[:]
    if (ns < 2):
        Gm1=diff_(Gm)
        return Gm,Gm1
    ll=zeros_(1,ns)
    lr=zeros_(1,ns)
    hl=zeros_(1,ns)
    hr=zeros_(1,ns)
    Gm1=diff_(Gm)
    return Gm,Gm1
def morphgrid1_(L=None,R=None,hL=None,hR=None,*args,**kwargs):
    varargin = cellarray(args)
    nargin = 4-[L,R,hL,hR].count(None)+len(args)

    dX=R - L
    n,q=stepcalc_(dX,hL,hR,nargout=2)
    G=zeros_(1,n - 1)
    g=copy_(L)
    for k in arange_(1,n - 1).reshape(-1):
        h=q[1] * k ** 2 + q[2] * k + q[3]
        g=g + h
        G[k]=g
    return G
def stepcalc_(dX=None,hL=None,hR=None,*args,**kwargs):
    varargin = cellarray(args)
    nargin = 3-[dX,hL,hR].count(None)+len(args)

    n=2 * dX / (hR + hL)
    n=floor_(n)
    q=step_co_(n,dX,hL,hR)
    return n,q
def step_co_(n=None,dX=None,hL=None,hR=None,*args,**kwargs):
    varargin = cellarray(args)
    nargin = 4-[n,dX,hL,hR].count(None)+len(args)

    nn=arange_(1,n)
    sn1=sum_(nn)
    sn2=sum_(nn ** 2)
    A=matlabarray([[sn2,sn1,n],[(n + 1) ** 2,(n + 1),1],[0,0,1]])
    q=linsolve_(A,[[dX],[hR],[hL]])
    a=q[1]
    b=q[2]
    c=q[3]
    return q
