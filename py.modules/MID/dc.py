#
import copy
import numpy as np

norm=np.linalg.norm;

def heron_sf(a,b,c):
    p=0.5*(a+b+c);
    ss=p*(p-a)*(p-b)*(p-c);
    return 0.0 if ss<=0.0 else np.sqrt(ss);

def heron(a,b,c):
    p=0.5*(a+b+c);
    return np.sqrt(p*(p-a)*(p-b)*(p-c));

def dc_profile(r,h,R):
    
    if np.size(r)>1:
        return np.array([dc_profile(v,h,R) for v in r]);
    
    #if h>=R:  
    #   print('warnimg++',h,R)
    if r>(h+R):
        return 0.0;
    if R>=(h+r):
        return 1.0;
    
    s=heron_sf(r,h,R);

    if s<=0.0:
        print('s<0')
        return 0.0;
    
    if np.abs(np.imag(s)>2.5e-16):
        return 0.0;
    
    fi=np.arcsin(2*s/(r*h))
    z=np.real(fi)/np.pi
    if R*R>r*r+h*h:
        z=1-z;
    #return np.sqrt(z);
    return z;


    
def rc_trs(trs,vxs):
    return (vxs[trs[:,0],0]+vxs[trs[:,1],0]+vxs[trs[:,2],0])/3.0;

def mu_sigma_profile3(trs_mask,trs,vxs,dc,Ri,Ro,fa,fc,bmu,sigma=0,fbmu=1):    
    
    trs=trs.reshape((int(trs.size/3),3));

    n=trs_mask.size;
    ii0=np.arange(n);    
    ii=ii0[trs_mask];
    i0=ii[0];    
    print('tm=trs[trs_mask,:];',n,trs.shape)
    tm=trs[trs_mask,:];
    print('tm.astype')
    tm=tm.astype(dtype='uint32',copy=False);
    r=rc_trs(tm,vxs);
    #dc=0.0
    pfo=dc_profile(r,dc,Ro)
    pfi=dc_profile(r,dc,Ri)    
    pf=pfo-pfi;
    #pf=dc_profile(r,dc,Ri) ;   
    #pf0=dc_profile(r,.0,Ro)-dc_profile(r,.0,Ri)    
    #fmc=sigma*pf;
    print('bmu=',bmu,'mu=',1./bmu);
    
    if fbmu:
        #fma=bmu*pf;
        xi=bmu-1.0
        fma=1+xi-xi*pfi;
        #fma=1.0+(bmu-1.0)*(pfo-pfi);
    else:
        mu=1.0/bmu;
        fma=1.0+(mu-1.0)*pf;
        fma=1.0/fma;

    #fma[fma<=0]=1.0;
    print(fma)
    fc[ii]=sigma;
    fa[ii]=fma;
    #fa[ii]=bmu;
    return (fa,fc);


def mu_sigma_profile(trs_mask,trs,vxs,dc,Ri,Ro,fa,fc,fbmu=1):    
    
    trs=trs.reshape((int(trs.size/3),3));

    n=trs_mask.size;
    ii0=np.arange(n);    
    ii=ii0[trs_mask];
    i0=ii[0];
    (bmu,sigma)=(fa[i0],fc[i0]);
    print('tm=trs[trs_mask,:];',n,trs.shape)
    tm=trs[trs_mask,:];
    print('tm.astype')
    tm=tm.astype(dtype='uint32',copy=False);
    r=rc_trs(tm,vxs);
    #dc=0.0
    pf=dc_profile(r,dc,Ro)-dc_profile(r,dc,Ri)    
    pf0=dc_profile(r,.0,Ro)-dc_profile(r,.0,Ri)    
    fmc=sigma*pf;
    print('bmu=',bmu,'mu=',1./bmu);
    if fbmu:
        fma=1.0+(bmu-1.0)*pf;
    else:
        mu=1.0/bmu;
        fma=1.0+(mu-1.0)*pf;
        fma=1.0/fma;

    fc[ii]=fmc;
    fa[ii]=fma;
    return (fa,fc);

def mu_sigma_profile2(trs_mask,trs,vxs,dc,Ri,Ro,sdc,sRi,sRo,fa,fc,fbmu=1):    
    
    trs=trs.reshape((int(trs.size/3),3));

    n=trs_mask.size;
    ii0=np.arange(n);    
    ii=ii0[trs_mask];
    i0=ii[0];
    (bmu,sigma)=(fa[i0],fc[i0]);
    print('tm=trs[trs_mask,:];',n,trs.shape)
    tm=trs[trs_mask,:];
    print('tm.astype')
    tm=tm.astype(dtype='uint32',copy=False);
    r=rc_trs(tm,vxs);
    
    pf=dc_profile(r,dc,Ro)-dc_profile(r,dc,Ri)    
    pfs=dc_profile(r,sdc,sRo)-dc_profile(r,sdc,sRi)    
    fmc=sigma*pfs;
    print('bmu=',bmu,'mu=',1./bmu);
    if fbmu:
        fma=1.0+(bmu-1.0)*pf;
    else:
        mu=1.0/bmu;
        fma=1.0+(mu-1.0)*pf;
        fma=1.0/fma;

    fc[ii]=fmc;
    fa[ii]=fma;
    return (fa,fc);

def mu_sigma_profile_ex(trs_mask,trs,vxs,dcL,RiRo,fac,fbmu=1):
    (Ri,Ro)=RiRo;
    (dc,L)=dcL;
    (fa,fc)=fac;


    pass

      