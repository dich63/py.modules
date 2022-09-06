# -*- coding: utf-8 -*-
"""
Created on Sun May 29 20:02:04 2022

@author: wwww
"""
from utils import *
from jsonrpc.jsonclass import *
from waves.fluid_solid_1FEM import *
import os,copy


def _jso_extend(d1,d2=None):    
    d={}
    if d2:
        d.update(as_dict(d2));
    if d1:
        d.update(as_dict(d1));
    return jso(d)


def reparse_mesh(mesh):
    #print('reparse_mesh:...',type(mesh))    
    if type(mesh)==str:
        l=6;
        mesh=mesh.strip()
        if mesh.find(':file:',0,l)==0:            
            l= l+1 if mesh[l:l+1]==':' else l;
            mesh=open(os.path.abspath(mesh[l:]),'r').read()
                
        mesh=decode(mesh,jslike=True)
        return mesh;
    
    #return copy.deepcopy(mesh)
    return mesh

def model2ext(model):
    
    model=reparse_mesh(model);
    
    default_params=_jso_extend(model.default_params,default_params_sf_1D())
    
    _ext=lambda x : _jso_extend(x,default_params)
    
    model.main_region.params=_ext(model.main_region.params);
    
    #regions_media=[_ext(r.media) for r in model.regions]
    model.regions=model._def("regions",[])
        
    for r in model.regions:
        r.params=_ext(r.params);   
    
    return model

model_expand=model2ext

def model2FEM1D(model,N=None,dx=None):
    
    
    
    def _get_KMG(NN,params,dx):
        
        v_p=params.v_p
        
        if v_p is not None:
            Ks=params.K_s;
            #params=jso(params);            
            params.rho_s=Ks/v_p**2
            
        FSM=make_sf_1D(**to_dict(params));
        mK,mM,mG=make_sf_FEM_data(NN,FSM,dx);#??!!!!!!!!!!!!!++
        '''
        mK,mM,mG,_=make_sf_1D_rgn(NN,dx,**to_dict(params))
        '''
        return mK,mM,mG;
    
    
    
    
    model=model2ext(model)
    
    pgn=model.main_region.polygon    
    zb,ze=pgn._def('zb',0.0),pgn.ze
    
        
    if dx is None:        
        if N is None:
            raise Exception('Either dx or N must be defined');                         
        L,dx=np.linspace(zb,ze,N,retstep=True,endpoint=False)
        index=np.arange(N);        
    else:
        N=int((ze-zb)/dx);
        index=np.arange(N);
        L=dx*index;        
    
        
    
    
    
    def _get_mask(polygon):
        b=polygon._def('zb',-np.inf);
        e=polygon._def('ze',np.inf);
        return index[(b<=L) & (L <= e)]
        
    
    
    mK,mM,mG=_get_KMG(N,model.main_region.params,dx);
    
    for rg in model.regions:
        im=_get_mask(rg.polygon)
        #print(rg.name,'=',im)
        mK[im],mM[im],mG[im]=_get_KMG(len(im),rg.params,dx);
        
    return mK,mM,mG,dx,N,model;
    
    
     


if __name__ == '__main__':
    
    #md=decode_from_file('model1D.json',1)
    #mK,mM,mG,dx=model2FEM1D(md,1000)
    fn=':file:'+__file__+'/../model1D.json'
    #fn=':file:'+__file__+'/../mod0.json'
    
    mK,mM,mG,dx,N,model=model2FEM1D(fn,16000);
    
    [smK,smG,smM]=make_sf_FEM([mK,mG,mM])
    #model.main_region.polygon.zb=-np.inf
    #encode_to_file(model,__file__+'/../model1D_full.json')