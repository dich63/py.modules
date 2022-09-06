import cmath
import numpy as np
from skimage.transform import resize
import k3d

cmps=k3d.matplotlib_color_maps;
#def Graph3D(Farray,Xarray=np.array([np.float32(x/10.0) for x in range (-10,10)]),Yarray=np.array([np.float32(y/10.0) for y in range (-10,10)]),axes=['T', 'z', 'Amp'], color_map=k3d.matplotlib_color_maps.Plasma ):
def _k3D_surf(Farray,Xarray,Yarray,axes, color_map,Xmm,Ymm ):    
    
    plot = k3d.plot(axes=axes,grid_auto_fit=True);
    #xmin, xmax,= np.min(Xarray.flatten()),np.max(Xarray.flatten()) ;
    #ymin, ymax,= np.min(Yarray.flatten()),np.max(Yarray.flatten()) ;    
    xmin, xmax= Xmm;
    ymin, ymax= Ymm;    
    
    cm = color_map; 
    z=1.*Farray/np.max(Farray.flatten());    
    z=z.real.astype(np.float32);
    plt_surface = k3d.surface(z, color_map=cm, attribute=z, color_range=[np.min(z.flatten()),np.max(z.flatten())], xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax)
    #plt_surface = k3d.surface(z, color_map=cm, attribute=z, color_range=[np.min(z.flatten()),np.max(z.flatten())],bounds=[-.7, .7, -4, 5,11,111], xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax)
    #plt_surface = k3d.surface(z, color_map=cm, attribute=z, color_range=[np.min(z.flatten()),np.max(z.flatten())],bounds=[-.7, .7, -4, 5,11,111],scaling=[1,1,1/10.], xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax)
    #plt_surface = k3d.surface(z, color_map=cm, attribute=z, color_range=[np.min(z.flatten()),np.max(z.flatten())],bounds=[-1, 1, -1, 1])
    #plt_surface.grid=(-40, 40,ymin,ymax)
    #plot.grid=(-40, 40,ymin,ymax,np.min(z), np.max(z));
    scale_ratio=np.abs(xmax-xmin)/np.abs(ymax-ymin);
    if scale_ratio<1:
        scale_ratio=1/scale_ratio;
        
   #plt_surface.model_matrix=[[1.,  0.,  0.,  0.],
    #         [ 0., 1.0*scale_ratio,  0.,  0.],
    #         [ 0.,  0.,  0.1,  0.],
    #         [ 0.,  0.,  0.,  1]];
    plot += plt_surface;
    plot.render();
    #plot.display();
    return plot;


def surf(F=None,XY=None,axes=['X', 'Y', 'F'], cmp=cmps.jet,Xmm=None,Ymm=None ):
    
    if XY is None:        
        XY=np.meshgrid(np.linspace(-1.,1.,200),np.linspace(-2.,2.,300))
        
    X,Y=XY;
    
    if F is None:        
        F=X*X+Y**3;
        
    if callable(F):
        F=F(X,Y);
        
    if Xmm is None:
        Xmm=np.min(X.flatten()),np.max(X.flatten()) ;
        
    if Ymm is None:
        Ymm=np.min(Y.flatten()),np.max(Y.flatten()) ;
        
    
        
    return _k3D_surf(F,X,Y,axes, cmp,Xmm,Ymm);