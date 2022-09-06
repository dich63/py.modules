# -*- coding: utf-8 -*-
"""
Created on Fri Mar 25 03:51:04 2016

@author: dich
from  MID.rect_triangulation import make_triangulation_c
"""
import sys
import numpy as np
from scipy.spatial import Delaunay

if sys.version_info.major==3:
    xrange=range

def expand1(x):
    x=np.array(x,dtype=np.float64);
    n=x.size;
    x=x.reshape(n);
    x2=np.zeros(2*n-1);
    xc=(x[1:]+x[:-1])/2.0;
    for k in range(n-1):
        x2[2*k]=x[k];
        x2[2*k+1]=xc[k];
    x2[2*n-2]=x[n-1];
    return x2;
    
        
    

def make_triangulation_c(rl,zl,lbound=0,fbc_mask=(1,1,1,1),unit_scale=1.0,info={}):

    def centroid(f):
        f=np.array(f,dtype='d');
        return (f[0:-1]+f[1:])/2.        

    def make_indexer(X,offset=0):
        ind=np.empty_like(X,dtype=np.uint32);
        ind.ravel()[:]=range(0,ind.size)
        ind+=offset;
        return ind;

    def make_vxs(X, Y):
        t = np.zeros( (np.size(X), 2) ,dtype='d')
        t[:,0] = X.reshape( np.size(X) )
        t[:,1] = Y.reshape( np.size(Y) )
        t = t.reshape(np.size(t))
        return t

    def make_bc_c(rl,zl,fbc_mask):

        Nz=np.size(zl);
        Nr=np.size(rl);
        
        fbc_c=np.zeros((Nz-1,Nr-1),dtype=np.uint16)

        fbc=np.zeros((Nz,Nr),dtype=np.uint16)

        fbc[ :, 0]=fbc_mask[0]
        fbc[ :,-1]=fbc_mask[1]
        fbc[ 0, :]=fbc_mask[2]
        fbc[-1, :]=fbc_mask[3]

        fbc=np.concatenate((fbc.ravel(),fbc_c.ravel()));

        return fbc


    def rect_c_triangulation(rl,zl):
                
        (rl,zl)=[np.array(k,dtype='d').ravel() for k in (rl,zl)]

        [X, Y] = np.meshgrid(rl,zl);

        rc,zc=centroid(rl),centroid(zl);
        [Xc, Yc] = np.meshgrid(rc,zc);

        ind=make_indexer(X);
        ind_c=make_indexer(Xc,ind.size);


        (Nzc,Nrc)=ind_c.shape;
        trs=np.zeros((Nzc,Nrc,4,3),dtype=np.uint32)
        for nz in xrange(Nzc):
            for nr in xrange(Nrc):
                tq=trs[nz,nr,:,:]
                ic=ind_c[nz,nr]
                i00=ind[nz,nr]
                i01=ind[nz,nr+1]
                i10=ind[nz+1,nr]
                i11=ind[nz+1,nr+1]
                tq[0]=(ic,i10,i00)
                tq[1]=(ic,i00,i01)
                tq[2]=(ic,i01,i11)
                tq[3]=(ic,i11,i10)
                
        vxs=make_vxs(X, Y)
        vxs_c=make_vxs(Xc, Yc)
        vxs=np.concatenate((vxs,vxs_c));        
        vxs =vxs.reshape((vxs.size//2,2))    
        trs=trs.reshape((Nzc*Nrc*4,3));
        
        return (trs,vxs,ind,ind_c)
        
            

    fbc= make_bc_c(rl,zl,fbc_mask);
    (trs,vxs,ind,ind_c)=rect_c_triangulation(rl,zl);

    vxs/=unit_scale;
    if not lbound==0:
        trs+=np.int(lbound);

    vxs = np.ascontiguousarray(vxs);
    trs = np.ascontiguousarray(trs);
    fbc = np.ascontiguousarray(fbc);

    info['vxs']=vxs;     
    info['trs']=trs;
    info['fbc']=fbc;
    info['lbound']=lbound;

    #info['ind']=ind;
    #info['ind_c']=ind_c;


    return (trs,vxs,fbc)

def make_triangulation(rl,zl,lbound=0,fbc_mask=(1,1,1,1),unit_scale=1.0,info={}):
    
    def make_bc(nxny,fbc_mask):
        fbc=np.zeros(nxny,dtype=np.uint16)
        fbc[ :, 0]=fbc_mask[0]
        fbc[ :,-1]=fbc_mask[1]
        fbc[ 0, :]=fbc_mask[2]
        fbc[-1, :]=fbc_mask[3]
        return fbc

    
    [X, Y] = np.meshgrid(rl,zl)
    nxny = np.shape(X)
    fbc=make_bc(nxny,fbc_mask)
    t = np.zeros( (np.size(X), 2) )

    t[:,0] = X.reshape( np.size(X) )
    t[:,1] = Y.reshape( np.size(Y) )

    vxs = np.ascontiguousarray(t, dtype="d")
    vxs /= unit_scale

    tri = Delaunay(t)
    trs = np.ascontiguousarray(tri.simplices, dtype=np.uint32) #trs
    if lbound:
        trs+=lbound
        
    info['vxs']=vxs;     
    info['trs']=trs;
    info['fbc']=fbc;
    info['lbound']=lbound;
    return (trs,vxs,fbc)

def triangulation_create(rl,zl,lbound=0,fcentroid=0,fbc_mask=(1,1,1,1),unit_scale=1.0,info={}):
    if fcentroid:
        make_triangulation_c(rl,zl,lbound,fbc_mask,unit_scale,info)
    else:
        make_triangulation(rl,zl,lbound,fbc_mask,unit_scale,info)
    return info;
        
    