import numpy as np
import numpy.ctypeslib as npct
import ctypes
#from .ll import midlib,lib_path

NullPtr=ctypes.c_void_p(0)
#midlib = npct.load_library('MID2D.dll','V:/pycpp/MID2D/x64/Debug')
#midlib = npct.load_library('MID2D.dll','V:/pycpp/MID2D/x64/Release')
#
midlib = npct.load_library('MID2D.dll',__file__+'/../bin')
#print(__file__)

#midlib =ll.midlib

""""
ca=np.ascontiguousarray(a);
"""

pdouble_t=ctypes.POINTER(ctypes.c_double)
puint32_t=ctypes.POINTER(ctypes.c_uint32)
uint32_t=ctypes.c_uint32
uint16_t=ctypes.c_uint16
ptri_t=puint32_t
pfbc_t=ctypes.POINTER(uint16_t)
pvoid_t=ctypes.c_void_p

#tst=midlib.tst
_volume_2D_to_matrix_trimesh=midlib.volume_2D_to_matrix_trimesh
_volume_2D_to_matrix_trimesh.restype=ctypes.c_int
_volume_2D_to_matrix_trimesh.argtype=[uint32_t,\
    uint32_t,puint32_t,uint32_t,pdouble_t,pdouble_t]


_laplace_2D_to_matrix_trimesh=midlib.laplace_2D_to_matrix_trimesh
_laplace_2D_to_matrix_trimesh.restype=ctypes.c_int
_laplace_2D_to_matrix_trimesh.argtype=[pvoid_t,uint32_t,\
    uint32_t,puint32_t,pdouble_t,pdouble_t,\
    uint32_t,pdouble_t,pfbc_t,pdouble_t,pdouble_t,pdouble_t]

_laplace_2D_to_matrix_trimesh_mask=midlib.laplace_2D_to_matrix_trimesh_mask
_laplace_2D_to_matrix_trimesh_mask.restype=ctypes.c_int
_laplace_2D_to_matrix_trimesh_mask.argtype=[pvoid_t,uint32_t,\
    uint32_t,puint32_t,puint32_t,pdouble_t,\
    uint32_t,pdouble_t,pfbc_t,pdouble_t,pdouble_t,pdouble_t]


_laplace_2D_to_matrix_trimesh_anx_mask=midlib.laplace_2D_to_matrix_trimesh_anx_mask
_laplace_2D_to_matrix_trimesh_anx_mask.restype=ctypes.c_int
_laplace_2D_to_matrix_trimesh_anx_mask.argtype=[pvoid_t,uint32_t,\
    uint32_t,puint32_t,puint32_t,pdouble_t,\
    uint32_t,pdouble_t,pfbc_t,pdouble_t,pdouble_t,pdouble_t]



_vxs2D_prjtransform_in_index=midlib.vxs2D_prjtransform_in_index
_vxs2D_prjtransform_in_index.restype=ctypes.c_int
_vxs2D_prjtransform_in_index.argtype=[uint32_t,puint32_t,
                                     pdouble_t,pdouble_t,pdouble_t]

_mesh2D_in_rect_index=midlib.mesh2D_in_rect_index
_mesh2D_in_rect_index.restype=ctypes.c_int
_mesh2D_in_rect_index.argtype=[uint32_t,uint32_t,uint32_t,puint32_t,puint32_t,
                    pdouble_t,uint32_t,puint32_t,pdouble_t]


_vxs2D_in_rect_index=midlib.vxs2D_in_rect_index
_vxs2D_in_rect_index.restype=ctypes.c_int
_vxs2D_in_rect_index.argtype=[uint32_t,uint32_t,pdouble_t,puint32_t,
                              uint32_t,puint32_t,pdouble_t]


_matrix_form=midlib.matrix_form
_matrix_form.restype=ctypes.c_int
_matrix_form.argtype=[uint32_t,uint32_t,puint32_t,
                    pdouble_t,pdouble_t,pdouble_t]

_grad_form=midlib.grad_form
_grad_form.restype=ctypes.c_int
_grad_form.argtype=[uint32_t,uint32_t,puint32_t,
                    pdouble_t,pdouble_t,pdouble_t]


_i_grad_form=midlib.i_grad_form
_i_grad_form.restype=ctypes.c_int
_i_grad_form.argtype=[uint32_t,uint32_t,puint32_t,
                    pdouble_t,pdouble_t,pdouble_t]


_i_grad_form_rect=midlib.i_grad_form_rect
_i_grad_form_rect.restype=ctypes.c_int
_i_grad_form_rect.argtype=[uint32_t,uint32_t,puint32_t,
                    pdouble_t,pdouble_t,pdouble_t,pdouble_t]

_i_mmz=midlib.i_mmz
_i_mmz.restype=ctypes.c_int
_i_mmz.argtype=[uint32_t,uint32_t,puint32_t,
                    pdouble_t,pdouble_t,pdouble_t]

_i_mmz_mask=midlib.i_mmz_mask
_i_mmz_mask.restype=ctypes.c_int
_i_mmz_mask.argtype=[uint32_t,uint32_t,puint32_t,
                    pdouble_t,puint32_t,pdouble_t,pdouble_t]


def _flat_(x,dtype):
    return np.ascontiguousarray(x,dtype=dtype)

def grad_matrix(trs,vxs,index_base=0):
    trs=_flat_(trs,uint32_t);
    vxs=_flat_(vxs,np.double);
    nvx=int(vxs.size/2);
    ntrs=int(trs.size/3);
    ptrs=trs.ctypes.data_as(puint32_t);
    pvxs=vxs.ctypes.data_as(pdouble_t);
    ntrs33=ntrs*9;
    mgx=np.empty((ntrs33,),dtype=np.double)
    mgy=np.empty((ntrs33,),dtype=np.double)
    
    pmgx=mgx.ctypes.data_as(pdouble_t);
    pmgy=mgy.ctypes.data_as(pdouble_t);
    
    index_base=uint32_t(int(index_base))
    ntrs=uint32_t(int(ntrs));
    
    _grad_form(index_base,ntrs,ptrs,pvxs,pmgx,pmgy);
    
    return (mgx,mgy);

def grad_form_rect_xy(trs,vxs,rect,index_base=0):
    trs=_flat_(trs,uint32_t);
    vxs=_flat_(vxs,np.double);
    (nvx,t)=vxs.shape;
    (ntrs,t)=trs.shape;
    ptrs=trs.ctypes.data_as(puint32_t);
    pvxs=vxs.ctypes.data_as(pdouble_t);
    
    fx=np.zeros((nvx,),dtype=np.double)
    fy=np.zeros((nvx,),dtype=np.double)
    
    pfx=fx.ctypes.data_as(pdouble_t);
    pfy=fy.ctypes.data_as(pdouble_t);
    
    rect=np.array(rect,dtype=np.float64);
    prect=rect.ctypes.data_as(pdouble_t);
    
    index_base=uint32_t(int(index_base))

    ntrs=uint32_t(int(ntrs));
    
    _i_grad_form_rect(index_base,ntrs,ptrs,pvxs,prect,pfx,pfy);
    
    return (fx,fy);

def mmz_form(trs,vxs,mu,index_base=0):
    trs=_flat_(trs,uint32_t);
    vxs=_flat_(vxs,np.double);
    #(nvx,t)=vxs.shape;
    #(ntrs,t)=trs.shape;
    nvx=np.uint32(vxs.size/2);
    ntrs=np.uint32(trs.size/3);
    ptrs=trs.ctypes.data_as(puint32_t);
    pvxs=vxs.ctypes.data_as(pdouble_t);

    pmu=mu.ctypes.data_as(pdouble_t);
    
    mmz=np.zeros((nvx,),dtype=np.double)
    
    
    pmmz=mmz.ctypes.data_as(pdouble_t);  
    
    
    index_base=uint32_t(int(index_base))

    ntrs=uint32_t(int(ntrs));
    
    _i_mmz(index_base,ntrs,ptrs,pvxs,pmu,pmmz);
    
    return mmz;

def mmz_form_mask(trs,vxs,trs_mask,mu,index_base=0):

    trs=_flat_(trs,uint32_t);
    vxs=_flat_(vxs,np.double);
    trs_mask=_flat_(trs_mask,uint32_t);
    mu=_flat_(mu,np.double);

    (nvx,t)=vxs.shape;
    (ntrs,t)=trs.shape;

    ptrs=trs.ctypes.data_as(puint32_t);
    pvxs=vxs.ctypes.data_as(pdouble_t);

    ptrs_mask=trs_mask.ctypes.data_as(puint32_t);

    pmu=mu.ctypes.data_as(pdouble_t);

    
    mmz=np.zeros((nvx,),dtype=np.double)
    
    
    pmmz=mmz.ctypes.data_as(pdouble_t);  
    
    
    index_base=uint32_t(int(index_base))

    ntrs=uint32_t(int(ntrs));
    
    _i_mmz_mask(index_base,ntrs,ptrs,pvxs,ptrs_mask,pmu,pmmz);
    
    return mmz;
    
    
    #(tmp,)=trs.shape
    
def volume_to_matrix_trimesh(trs,vxs,index_base=1):

    trs=_flat_(trs,uint32_t)    
    vxs=_flat_(vxs,'d')
    
    ptrs=trs.ctypes.data_as(puint32_t);    
    pvxs=vxs.ctypes.data_as(pdouble_t);

    nvxs=int(vxs.size/2);
    ntrs=int(trs.size/3);
    ntrs33=ntrs*9    

    mv=np.empty((ntrs33,),dtype=np.double)    
    
    pmv=mv.ctypes.data_as(pdouble_t);
    
    _volume_2D_to_matrix_trimesh(uint32_t(int(index_base)),uint32_t(int(ntrs)),ptrs,uint32_t(int(nvxs)),pvxs,pmv);
    
    return mv


def to_matrix_trimesh(trs,fa,fc,vxs,fbc,index_base=1,scheme=1):
    trs=_flat_(trs,uint32_t)
    fa=_flat_(fa,'d')
    fc=_flat_(fc,'d')
    vxs=_flat_(vxs,'d')
    fbc=_flat_(fbc,uint16_t)
    ntrs=fa.size;
    nvxs=fbc.size
    if not ((ntrs==fc.size) and (3*ntrs==trs.size)):
        raise Exception('invalid triangle data size');
    if not (2*nvxs==vxs.size):
        raise Exception('invalid vertexes data size');
    ptrs=trs.ctypes.data_as(puint32_t);
    pa=fa.ctypes.data_as(pdouble_t);
    pc=fc.ctypes.data_as(pdouble_t);
    
    pvxs=vxs.ctypes.data_as(pdouble_t);
    pfbc=fbc.ctypes.data_as(pfbc_t);
    ntrs33=ntrs*9
    ma=np.zeros((ntrs33,),dtype=np.double)
    mc=np.zeros((ntrs33,),dtype=np.double)
    mv=np.zeros((ntrs33,),dtype=np.double)
    
    pma=ma.ctypes.data_as(pdouble_t);
    pmc=mc.ctypes.data_as(pdouble_t);
    pmv=mv.ctypes.data_as(pdouble_t);
    
    _laplace_2D_to_matrix_trimesh(pvoid_t(scheme),uint32_t(int(index_base)),uint32_t(int(ntrs)),ptrs,pa,pc,uint32_t(int(nvxs)),pvxs,pfbc,pma,pmc,pmv);
    
    return (ma,mc,mv)

def to_matrix_trimesh_mask(trs,mask,fac,vxs,fbc,index_base=1,scheme=1):
    trs=_flat_(trs,uint32_t)    

    fac=np.array(fac).transpose();    
    nanx=fac.shape[1];
    fac=_flat_(fac,'d');

    mask=_flat_(mask,uint32_t);

    vxs=_flat_(vxs,'d')
    fbc=_flat_(fbc,uint16_t)

    ntrs=mask.size;
    nvxs=fbc.size

    if not ((3*ntrs==trs.size)):
        raise Exception('invalid triangle data size');
    if not (2*nvxs==vxs.size):
        raise Exception('invalid vertexes data size');

    ptrs=trs.ctypes.data_as(puint32_t);
    pac=fac.ctypes.data_as(pdouble_t);
    pmask=mask.ctypes.data_as(puint32_t);
    
    
    pvxs=vxs.ctypes.data_as(pdouble_t);
    pfbc=fbc.ctypes.data_as(pfbc_t);
    ntrs33=ntrs*9
    ma=np.zeros((ntrs33,),dtype=np.double)
    mc=np.zeros((ntrs33,),dtype=np.double)
    mv=np.zeros((ntrs33,),dtype=np.double)
    
    pma=ma.ctypes.data_as(pdouble_t);
    pmc=mc.ctypes.data_as(pdouble_t);
    pmv=mv.ctypes.data_as(pdouble_t);
    
    if nanx==2:
        _laplace_2D_to_matrix_trimesh_mask(pvoid_t(scheme),uint32_t(int(index_base)),uint32_t(int(ntrs)),ptrs,pmask,pac,uint32_t(int(nvxs)),pvxs,pfbc,pma,pmc,pmv);
    elif nanx==3:
        _laplace_2D_to_matrix_trimesh_anx_mask(pvoid_t(scheme),uint32_t(int(index_base)),uint32_t(int(ntrs)),ptrs,pmask,pac,uint32_t(int(nvxs)),pvxs,pfbc,pma,pmc,pmv);
    else:
        raise Exception('invalid fACX size');

    
    return (ma,mc,mv)

def vxs2D_projective_deformation(masks,matrix3x3,vxs,vxs_out=None):
    masks=_flat_(masks,dtype=uint32_t)
    matrix3x3=_flat_(matrix3x3,dtype=np.double)
    vxs=_flat_(vxs,dtype=np.double)
    #nvxs=np.prod(masks.shape)
    nvxs=masks.size
    #if not (9*nvxs==np.prod(matrix3x3.shape)):
    #    raise Exception('invalid matrix3x3 data size');
    if not (2*nvxs==vxs.size):
        raise Exception('invalid vertexes data size');
    
    if vxs_out is None:
        vxs_out=np.empty_like(vxs)
    elif not (2*nvxs==vxs_out.size):
        raise Exception('invalid vertexes_out data size');
    
    pmasks=masks.ctypes.data_as(puint32_t);
    pm3x3=matrix3x3.ctypes.data_as(pdouble_t);
    pvxs_in=vxs.ctypes.data_as(pdouble_t);
    pvxs_out=vxs_out.ctypes.data_as(pdouble_t);
    
    _vxs2D_prjtransform_in_index(uint32_t(int(nvxs)),pmasks,pm3x3,pvxs_in,pvxs_out)
    return vxs_out


def mesh2D_in_rect_index(trs,vxs,rects_indx,rects,masks=None,fcenter=True,index_base=1):
    '''
    trs=np.array(trs,dtype=uint32_t,copy=False)
    vxs=np.array(vxs,dtype=np.double,copy=False)
    rects_indx=np.array(rects_indx,dtype=uint32_t,copy=False)
    rects=np.array(rects,dtype=np.double,copy=False)
    '''
    trs=_flat_(trs,dtype=uint32_t)
    vxs=_flat_(vxs,dtype=np.double)
    rects_indx=_flat_(rects_indx,dtype=uint32_t)
    rects=_flat_(rects,dtype=np.double)

    ntrs=int(trs.size/3)
    if masks is None:
        masks=np.zeros((ntrs,),dtype=uint32_t);
    
    nrgn=rects_indx.size
    
    fcenter=uint32_t(fcenter)
    if not (ntrs==masks.size):
        raise Exception('invalid triangle data size');
    if not (4*nrgn==rects.size):
        raise Exception('invalid rects data size');

    index_base=uint32_t(int(index_base))
    ntrs=uint32_t(int(ntrs))
    nrgn=uint32_t(int(nrgn))
    
    ptrs=trs.ctypes.data_as(puint32_t);
    pmasks=masks.ctypes.data_as(puint32_t);
    prects_indx=rects_indx.ctypes.data_as(puint32_t);
    
    pvxs=vxs.ctypes.data_as(pdouble_t);
    prects=rects.ctypes.data_as(pdouble_t);
    
    _mesh2D_in_rect_index(index_base,fcenter,ntrs,ptrs,pmasks,pvxs,nrgn,prects_indx,prects)
    
    return masks

def mesh2D_in_rect_index2(trs,vxs,ri_pair,masks=None,fcenter=True,index_base=1):
    return mesh2D_in_rect_index(trs,vxs,ri_pair[0],ri_pair[1],masks,fcenter,index_base)

def vxs2D_in_rect_index(vxs,rects_indx,rects,mode=1,masks=None):
    '''
    vxs=np.array(vxs,dtype=np.double,copy=False)
    rects_indx=np.array(rects_indx,dtype=uint32_t,copy=False)
    rects=np.array(rects,dtype=np.double,copy=False)
    '''

    vxs=_flat_(vxs,dtype=np.double)
    rects_indx=_flat_(rects_indx,dtype=uint32_t)
    rects=_flat_(rects,dtype=np.double)
    
    nvxs=int(vxs.size/2)
    if masks is None:
        masks=np.zeros((nvxs,),dtype=uint32_t);
    #nrgn=np.prod(rects_indx.shape)
    nrgn=rects_indx.size
    if not (4*nrgn==rects.size):
        raise Exception('invalid rects data size')
    pmasks=masks.ctypes.data_as(puint32_t);
    prects_indx=rects_indx.ctypes.data_as(puint32_t);    
    pvxs=vxs.ctypes.data_as(pdouble_t);
    prects=rects.ctypes.data_as(pdouble_t);
    mode=uint32_t(int(mode))
    nvxs=uint32_t(int(nvxs))
    nrgn=uint32_t(int(nrgn))
    _vxs2D_in_rect_index(mode,nvxs,pvxs,pmasks,nrgn,prects_indx,prects);
    return masks

def vxs2D_in_rect_index2(vxs,ri_pair,rects,masks=None):
    return vxs2D_in_rect_index(vxs,ri_pair[0],ri_pair[1],masks)
