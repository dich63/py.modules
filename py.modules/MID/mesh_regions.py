# -*- coding: utf-8 -*-
"""
Created on Mon Jan 25 00:35:08 2016

@author: dich
"""


import copy
import numpy as np
import sys
#import jsonrpc.jsonclass as jc
from .clib import *
from .clib_ex import *
from .projective2d import quad2quad_projective,rect2quad
import  lipa.trisparse as trisparse
import jsonrpc.jsonclass as jc
import scipy.sparse as ssp
import MID.dc as dcc

norm=np.linalg.norm;

normm= lambda x : norm(x,np.inf)


epsilon=sys.float_info.epsilon

def set_region_indexes(regions):
    i=0
    for rgn in regions:
        i+=1
        rgn.index=i

def get_regions_params(regions):
    pass

def is_rect_eq(r1,r2):
    r1=np.array(r1,dtype='d')
    r2=np.array(r2,dtype='d')
    return norm(r1-r2)<=epsilon

def region_radius_mask(trs,vxs,rc):
    def rmask(v):
        z=v[:,0]+1j*v[:,1];
        f=np.abs(z)>rc;
        return f;

    trs=np.reshape(trs,[int(trs.size/3),3]);

    fvx=rmask(vxs)

    v0=vxs[trs[:,0],:];
    f0=rmask(v0);

    v1=vxs[trs[:,1],:];
    f1=rmask(v1);

    v2=vxs[trs[:,2],:];
    f2=rmask(v2);

    fin=~ (f0 | f1 | f2);
    fout= (f0 & f1 & f2);

    fbc=~(fin | fout );

    return (fin,fout,fbc,fvx);

def region_by_name(rgs,name,dfl=None):    
    if not not name:
        for r in rgs:
            if r['name']==name:
                return r;
    return dfl;

    
def trs_out_mask(trs,vxs,rc):
    (fin,fout,fbc,fvx)=region_radius_mask(trs,vxs,rc);
    trs=np.reshape(trs,[int(trs.size/3),3]);
    return trs[fout,:];
    
def trs_ring_mask(trs,vxs,rin=-1.0,rout=np.inf):
    (bin,bout,bbc,bvx)=region_radius_mask(trs,vxs,rin);
    (ein,eout,ebc,evx)=region_radius_mask(trs,vxs,rout);
    f=bout & ein;
    trs=np.reshape(trs,[int(trs.size/3),3]);
    return trs[f,:];

class MID_regions:
    def __init__(self,regions):
        self.regions=copy.deepcopy(regions);
        d={}
        k=0;
        for rgn in self.regions:
            k+=1;
            rgn['index']=k;
            d[rgn['name']]=rgn;
            rects=rgn['rects']
            for n in range(len(rects)):
                r=rects[n]
                if type(r)== np.ndarray:
                    rects[n]=r.reshape(4).tolist()
                rects[n]=[np.float64(x) for x in rects[n]]
        self.rgnmap=d;
    def __iter__(self):
        return iter(self.regions);
    def __getitem__(self,name):
        return self.rgnmap[name];

class MID_mesh:
    def __init__(self,mesh):
        self.mesh=mesh
        self.regions0=MID_regions(mesh['regions'])
        self.regions=MID_regions(mesh['regions'])
        self.fbc=mesh['fbc'];
        self.vxs0=mesh['vxs'];
        self.sigma_mu_bg=mesh.get('sigma_mu_bg',(0.00,1.0))
        #self.vxs=self.vxs0.copy();
        self.vxs=np.zeros_like(self.vxs0);
        self.trs=mesh['trs'];
        self.index_base=self.mesh.get('lbound',1)
        u=mesh.get('units',{'regions':1.0});
        self.unit_r=np.double(u['regions']);
        self.beta_0=np.double(mesh.get('beta_0',1.25663681152658e-06))
        #self.rect_index_pair0=self.rect_index_pair
        self.fuzzy_bc=fuzzy_bc_t(self.vxs0,self.trs);
        
    def reset_regions(self):
        self.regions=MID_regions(self.regions0)
    def reset_vertexes(self):
        self.vxs=self.mesh['vxs'];
    def rescale_rect(self,rect):
        unit_r=self.unit_r
        xb=unit_r*rect[0];
        yb=unit_r*rect[1];
        xe=unit_r*rect[2];
        ye=unit_r*rect[3];
        return (xb,xe,yb,ye);
    @property
    def rect_index_pair(self):
        return self.get_rect_index_pair(self.regions,clsdis=False)
    @property
    def rect_index_pair0(self):
        return self.get_rect_index_pair(self.regions0,clsdis=False)


    def get_rect_index_pair(self,regions,clsdis=False):
        unit_r=self.unit_r
        (rindexes,rects)=([0],[(-np.Inf,np.Inf,-np.Inf,np.Inf)])
        for rgn in regions:
            if not  ( clsdis and rgn.get('disabled',False)):
                indx=rgn['index'];
                for rect in rgn['rects']:
                    rindexes.append(indx)
                    rects.append(self.rescale_rect(rect));
        return (rindexes,rects);

    def sigma_mu_pair_anx(self,sigma_bg=0.0,mu_bg=1):

        beta_0=-np.double(self.beta_0)
        (sigmas,bmus,anx)=([beta_0*sigma_bg],[1./np.double(mu_bg)],[1.0])
        
        for rgn in self.regions:
            f=not rgn.get('disabled',False)
            #f=0;
            (sigma,mu)=(rgn.get('sigma',sigma_bg),rgn.get('mu',mu_bg)) if f else (sigma_bg,mu_bg)
            mu=np.double(mu);
            anz=np.double(rgn.get('z2r',1.0));
            #for rect in rgn['rects']:
            sanz=np.sqrt(anz);
            bmus.append(1./(mu*sanz))
            anx.append(sanz);

            sigmas.append(beta_0*sigma)
            
        return (np.array(sigmas,dtype='d'),np.array(bmus,dtype='d'),np.array(anx,dtype='d'));

    def sigma_mu_pair(self,sigma_bg=0.0,mu_bg=1):

        beta_0=-np.double(self.beta_0)        
        (sigmas,bmus)=([beta_0*sigma_bg],[1./np.double(mu_bg),1.0])
        
        for rgn in self.regions:
            f=not rgn.get('disabled',False)
            #f=0;
            (sigma,mu)=(rgn.get('sigma',sigma_bg),rgn.get('mu',mu_bg)) if f else (sigma_bg,mu_bg)
            #for rect in rgn['rects']:
            sigmas.append(beta_0*sigma)
            bmus.append(1./np.double(mu))
        return (np.array(sigmas,dtype='d'),np.array(bmus,dtype='d'));


    def decenter_profile3(self,trs_m,fa,fc,fbmu=1):       

        trs=self.trs;
        vxs=self.vxs;
        regions=self.regions;
        for rgn in regions:
            f=not rgn.get('disabled',False)
            if f:
                dcn=rgn.get('dc_pipe',None);
                f=not not dcn;
                if f:
                    rgc=region_by_name(regions,dcn,{});
                    dc=rgc.get('dcr',0.0);
                    f=dc>0;

                        
            if f:
                print('decenter_profile3',rgn)
                bmu=1.0/rgc.get('mu',1.0);
                indx=rgn['index'];
                ftm=trs_m==indx;
                ur=self.unit_r;

                rt=rgn['rects'][0];
                
                (Ri,Ro)=(rt[0],rt[2]);
                [Ri,Ro,dc]=[ur*v for v in [Ri,Ro,dc]];

                rtc=rgc['rects'][0];
                (sRi,sRo)=(rtc[0],rtc[2]);

                [sRi,sRo,sdc]=[ur*v for v in [sRi,sRo,0.0]];

                print('mu_sigma_profile start')
                print(ftm,trs,vxs,dc,Ri,Ro,fa,fc)
                dcc.mu_sigma_profile3(ftm,trs,vxs,dc,sRi,sRo,fa,fc,bmu,0,fbmu);                
                pass
        pass



    def decenter_profile(self,trs_m,fa,fc,fbmu=0):       
        
        for rgn in self.regions:
            f=not rgn.get('disabled',False)
            dc=rgn.get('dc',0);
                        
            trs=self.trs;
            vxs=self.vxs;
            if f and dc>0:
                indx=rgn['index'];
                ftm=trs_m==indx;
                ur=self.unit_r;

                rt=rgn['rects'][0];
                
                (Ri,Ro)=(rt[0],rt[2]);
                [Ri,Ro,dc]=[ur*v for v in [Ri,Ro,dc]];

                rtc=rgn['rects_c'][0];
                (sRi,sRo)=(rtc[0],rtc[2]);

                [sRi,sRo,sdc]=[ur*v for v in [sRi,sRo,0.0]];

                print('mu_sigma_profile start')
                print(ftm,trs,vxs,dc,Ri,Ro,fa,fc)
                dcc.mu_sigma_profile(ftm,trs,vxs,dc,sRi,sRo,fa,fc,fbmu);                
                pass
        pass
        

    @property
    def trs_mask(self):
        rect_index_pair0=self.get_rect_index_pair(self.regions0,clsdis=True);
        return mesh2D_in_rect_index2(self.trs,self.vxs0,rect_index_pair0,index_base=self.index_base);


    def tri2sparse(self,data):
        trs=self.trs;
        vxs=self.vxs;
        nvxs=int(vxs.size/2);
        return trisparse.trimesh2sparse(trs=trs,nvxs=nvxs,data=data,lbound=self.index_base)

    def rect_trs_mask(self,l=0,b=-np.inf,r=np.inf,t=np.inf):
        u=self.unit_r;
        pair=[[1],[u*l,u*r,u*b,u*t]];
        return mesh2D_in_rect_index2(self.trs,self.vxs0,pair,index_base=self.index_base);

    def rect_trs_mask2(self,r):
        return self.rect_trs_mask(r[0],r[1],r[2],r[3]);

    def make_volume(self):
        data=volume_to_matrix_trimesh(self.trs,self.vxs,self.index_base);
        mv=self.tri2sparse(data);
        return mv;

    def make_gradxy(self,mV):

        mV=mV.tocsc();
        vxs=self.vxs;
        n=np.int64(vxs.size/2);
        v=np.ones(n,dtype = np.float64);
        v=mV*v;
        mbv=ssp.diags(1./v,format='csc');

        gxy=grad_matrix(self.trs,vxs);

        gxy=[ self.tri2sparse(data) for data in gxy];
        gxy=[ mbv*Hi.tocsc() for Hi in gxy];
        #gxy=[ Hi.tocsc() for Hi in gxy];

        return gxy[0]+1j*gxy[1];
        

    def makeACVM(self,scheme=1,fuzzy=None):
        (fcm,fam,f_anx)=self.sigma_mu_pair_anx(self.sigma_mu_bg[0],self.sigma_mu_bg[1]);

        #fcm=1.0*fam;
        #fam=1.0*fam;
        #fcm[:]=1;
        #fam[:]=2;

        tm=self.trs_mask;
        fa=fam[tm];
        fc=fcm[tm];
        
            
            

        #fc[:]=0;
        (trs,vxs,index_base)=(self.trs,self.vxs,self.index_base)
        #fbc=copy.copy(self.fbc);
        #fbc[:]=0;
        #self.fbc=fbc;
        
        self.decenter_profile(tm,fa,fc);       
        print('decenter_profile end')

        self.decenter_profile3(tm,fa,fc,1);       
        print('decenter_profile3 end')

        #trs=np.uint32(np.round(trs))
        mmzf=mmz_form_mask(trs,vxs,tm,fam,index_base=index_base);
        #m0=mmz_form(trs,vxs,fa,index_base=index_base);
        
        #a0=to_matrix_trimesh(trs,fa,fc,vxs,self.fbc,index_base=index_base,scheme=scheme);
        #a0=to_matrix_trimesh_mask(trs,tm,(fam,fcm),vxs,self.fbc,index_base=index_base,scheme=scheme);
        #acv=to_matrix_trimesh_mask(trs,tm,(fam,fcm,f_anx),vxs,self.fbc,index_base=index_base,scheme=scheme);
        acv=to_matrix_trimesh(trs,fa,fc,vxs,self.fbc,index_base=index_base,scheme=scheme);
        #ACV=[self.tri2sparse(data).tocoo() for data in acv ]
        ACVM=[self.tri2sparse(data) for data in acv ]
        ACVM+=[mmzf];
        return ACVM

    def makeACV(self,scheme=1,fuzzy=None):
        (fcm,fam)=self.sigma_mu_pair(self.sigma_mu_bg[0],self.sigma_mu_bg[1]);
        tm=self.trs_mask;
        fa=fam[tm];
        fc=fcm[tm];
        
        if not fuzzy is None:
            a=self.fuzzy_bc.make(**fuzzy).morph_scale;
            expa=np.exp(a)
            #fc+= 0.1*self.beta_0*expa;
            fa/=expa;
            
            

        #fc[:]=0;
        #mmzf=mmz_form(self.trs,self.vxs,fa,index_base=self.index_base)
        acv=to_matrix_trimesh(self.trs,fa,fc,self.vxs,self.fbc,index_base=self.index_base,scheme=scheme);
        #ACV=[self.tri2sparse(data).tocoo() for data in acv ]
        ACV=[self.tri2sparse(data) for data in acv ]
        return ACV

    def makeAC(self,scheme=1):
        return self.makeACV(scheme)[0:2]

    def make_vxs_deformation(self):
        idT=((0,0,0),(0,0,0),(0,0,0))
        m3x3=[idT];
        (tmp,rs0)=self.rect_index_pair0;
        (tmp,rs)=self.rect_index_pair;
        Nr=len(rs0);
        i=0;
        
        indxs=[]
        rsd0=[]
        fdeform=False;
        for n in range(1,Nr):
            if not is_rect_eq(rs0[n],rs[n]):
                fdeform=True;
                print('deformation:',rs0[n],'->',rs[n])
                i+=1
                rsd0.append(rs0[n]);
                indxs.append(i);
                T=quad2quad_projective(rect2quad(rs0[n]),rect2quad(rs[n]))
                print(T)
                m3x3.append(T);

        if fdeform:
            masks=vxs2D_in_rect_index(self.vxs0,indxs,rsd0,mode=0);
            self.vxs=vxs2D_projective_deformation(masks,m3x3,self.vxs0);
        else:
            self.vxs[:]=self.vxs0;

        return self

    def region_mask_rect(self,name,eps=1e-4):
        def inflate_rect(r,eps):
            return [r[0]-eps,r[1]+eps,r[2]-eps,r[3]+eps];

        sr=self.unit_r;
        eps=eps*sr;
        if type(name) in jc.string_types:
            rs=self.regions0[name]['rects'];
        else:
            rs=name;
        rindx=np.ones(len(rs),dtype=uint32_t)
        rects=[inflate_rect(self.rescale_rect(r),eps) for r in rs];
        return (rindx,rects);

    def region_mask(self,name,eps=1e-4,cutoff_rect=[-np.inf,+np.inf,-np.inf,+np.inf]):
        """
        def inflate_rect(r,eps):
            return [r[0]-eps,r[1]+eps,r[2]-eps,r[3]+eps];

        sr=self.unit_r;
        eps=eps*sr;
        if type(name) in jc.string_types:
            rs=self.regions0[name]['rects'];
        else:

        rindx=np.ones(len(rs),dtype=uint32_t)
        rects=[inflate_rect(self.rescale_rect(r),eps) for r in rs];
        """
        """
        (rindx,rects)=([],[])
        for r in rs:
            rects.append(self.rescale_rect(r));
            rindx.append(1);
        """
        def ir(r1,r2):
            return [max(r1[0],r2[0]),min(r1[1],r2[1]),max(r1[2],r2[2]),min(r1[3],r2[3])];


        (rindx,rects)=self.region_mask_rect(name,eps);
        rects=[ir(r,cutoff_rect) for r in rects];
        return vxs2D_in_rect_index(self.vxs0,rindx,rects);

    def region_rect(self,name,eps=1.e-4):
        (ri,rs)=self.region_mask_rect(name,eps);        
        r=rs[0];
        return r; 

    def form_grad_xy(self,rect,xzero=True,eps=1.e-4,with_deform=False):

        vxs=self.vxs if with_deform else self.vxs0;

        if type(rect) is str:
            rect=self.region_rect(rect,eps);

            if xzero:
                rect[0]=-eps;

        else:
            rect=copy.copy(rect);
        


        (fx,fy)=grad_form_rect_xy(self.trs,vxs,rect);

        return (fx,fy);




    def form_volume(self,rgn,eps=1e-4,scheme=1,with_deform=0):
        (rindexes,rects)=([0],[(-np.Inf,np.Inf,-np.Inf,np.Inf)]);
        (ri,rs)=self.region_mask_rect(rgn,eps);        
        rects+=rs;
        rindexes+=ri.tolist();

        vxs=self.vxs if with_deform else self.vxs0;
        
        b0=vxs2D_in_rect_index(vxs,ri,rs);

        #b0=vxs2D_in_rect_index(vxs,ri,[[0.005, 0.5, -1.0, 1.0]]);

        tm=mesh2D_in_rect_index2(self.trs,vxs,(rindexes,rects),index_base=self.index_base,fcenter=0);
        fam=np.array([0.0,0.0],dtype=np.float64);
        fcm=np.array([0.0,1.0],dtype=np.float64);


        fa=fam[tm];
        fc=fcm[tm];

        tmp,c,tmp=to_matrix_trimesh(self.trs,fa,fc,self.vxs,self.fbc,index_base=self.index_base,scheme=scheme);
        C=self.tri2sparse(c).tocsc();

        n=C.shape[0];
        #b0=np.ones(n,dtype=np.float64);
        return C*b0;

'''
    def form_volume(self,rgn,eps=1e-4,scheme=1,with_deform=0):
        (rindexes,rects)=([0],[(-np.Inf,np.Inf,-np.Inf,np.Inf)]);
        (ri,rs)=self.region_mask_rect(rgn,eps);        
        rects+=rs;
        rindexes+=ri.tolist();

        vxs=self.vxs if with_deform else self.vxs0;
        
        b0=vxs2D_in_rect_index(vxs,ri,rs);

        #b0=vxs2D_in_rect_index(vxs,ri,[[0.005, 0.5, -1.0, 1.0]]);

        tm=mesh2D_in_rect_index2(self.trs,vxs,(rindexes,rects),index_base=self.index_base,fcenter=0);
        fam=np.array([0.0,0.0],dtype=np.float64);
        fcm=np.array([0.0,1.0],dtype=np.float64);
        #fcm=np.array([1.0,1.0],dtype=np.float64);

        fa=fam[tm];
        fc=fcm[tm];
        fbc=self.fbc;        
        tmp,c,v=to_matrix_trimesh(self.trs,fa,fc,self.vxs,fbc,index_base=self.index_base,scheme=scheme);
        C=self.tri2sparse(c).tocsc();
        
        n=C.shape[0];
        b0=np.ones(n,dtype=np.float64);
        V=self.tri2sparse(v).tocsc();
        rv=V*b0;
        r=C*b0;
        return r;
'''        