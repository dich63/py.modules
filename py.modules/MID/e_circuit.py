# -*- coding: utf-8 -*-
"""
Created on Mon Aug 20 15:03:40 2018

@author: dich
"""

import copy
import numpy as np
import sys
from MID.clib import *
from  lipa.trisparse import coo_scheme

from scipy import sparse as sp
from scipy.sparse import coo_matrix,csc_matrix
from lipa.trisparse import trimesh2sparse

_epsilon=sys.float_info.epsilon

norm=np.linalg.norm;

def _to_coo_scheme(Nf,col,row,data):
    col=np.concatenate(col);
    row=np.concatenate(row);
    data=np.concatenate(data);
    return coo_scheme((Nf,Nf),row=row,col=col,data=data);


def tolist(l):
    if type(l) in (tuple,list):
        return l;
    else:
        return (l,);

def nzmask(a,eps=1e-15):    
    return np.abs(a)>eps;

def field_Mz(vxs,Mz=1):
    z=v[:,0]+1j*v[:,1];
    return z*z/np.abs(z)**3; 


def addMz(A,C,C0,vxs,mMz=None,eps=_epsilon):
    iA=[A.col];    
    jA=[A.row];
    dA=[A.data[0:A.nnz]];
    iC=[C.col];    
    jC=[C.row];
    dC=[C.data[0:C.nnz]];
    (Nf,tmp)=vxs.shape;
    (N,tmp)=A.shape;
    N1=N+1
    r=np.zeros(N,dtype = np.float64);
    r[0:Nf]=vxs[:,1];
    C0=C0.tocsc();

    data=r*C0;
    i0=np.arange(N,dtype=np.uint32);
    f=nzmask(data,eps);
    data=data[f];    
    Nc=data.size;
    i=i0[f];
    Nc=i.size;
    j=N*np.ones(Nc,dtype=np.uint32);    
    
    iC+=[i];    
    jC+=[j];
    dC+=[data];
    
    ii=range(Nf,N);
    cc=C0[ii,:];
    y=2*cc*r;
    jj=N*np.ones(N-Nf,dtype=np.uint32);    
    
    iA+=[ii];    
    jA+=[jj];
    dA+=[y];
    
    
    
    if not mMz is None:
        print('norm mMz=',norm(mMz))
        f=nzmask(mMz,eps);
        data=mMz[f];
        i=i0[f];
        Nc=i.size;
        j=N*np.ones(Nc,dtype=np.uint32);    
        iA+=[i];
        jA+=[j];
        dA+=[data];
        pass
    


    iA+=[[N]];    
    jA+=[[N]];
    dA+=[[1.0]];

    Ax=_to_coo_scheme(N1,col=iA,row=jA,data=dA);
    Cx=_to_coo_scheme(N1,col=iC,row=jC,data=dC);
    return (Ax,Cx);



class coil_chain_t:
    def __init__(self,ACV,units={'a':1,'e':0.5,'r':1}):
        self.units=units;
        self.ACV=ACV;
        
        pass

    def make_coil(self,Nw,mask,nj,eps=1e-15):
        data=np.array(mask,dtype = np.float64 );
        N=mask.size;
        data=self.mV*data;
        f=nzmask(data,eps);
        data=data[f];
        zm=Nw/np.sum(data);        
        data*=zm;
        Nc=data.size;
        i=np.arange(N,dtype=np.uint32);
        i=i[f];
        j=nj*np.ones(Nc,dtype=np.uint32);
        return (data,i,j);
    def getunits(self):
        us=self.units
        return (us['a'],us['e'],us['r'])

    def make_ecj(self,jc,nj,eps=1e-15):
        mV=self.mV;
        N=mV.shape[0];
        dAs,iAs,jAs=[],[],[];
        dCs,iCs,jCs=[],[],[];        

        (ua,ue,ur)=self.getunits()

        coils=jc['coils']
        R=jc['R'];
        
        
        for cl in coils:
            Nw=cl['Nw'];
            mask=cl['mask'];
            mask=np.array(mask,dtype = np.uint32);
            (eta,i,j)=self.make_coil(Nw,mask,nj,eps);

            iAs+=[i];
            jAs+=[j];
            dAs+=[ua*eta];

            iCs+=[j];
            jCs+=[i];
            dCs+=[ue*eta];
        
        aR=ur*np.array([R],dtype=np.float64);
        anj=np.array([nj],dtype=np.uint32); 

        cc=jc.get('C',0.0);
        fcc=not not cc;

        if fcc:
            
            (dCs0,iCs0,jCs0)=[ copy.copy(c)  for c in   (dCs,iCs,jCs)];

            mCR=+cc*R;

            dCCs=[mCR*dA for dA in dAs]            
            dCs+=dCCs;
            iCs+=iAs;
            jCs+=jAs;
        else:
            (dCs0,iCs0,jCs0)=(dCs,iCs,jCs);




        iAs+=[anj];
        jAs+=[anj];
        dAs+=[aR];

        return ((iAs,jAs,dAs),(iCs,jCs,dCs),(iCs0,jCs0,dCs0))


        

    def make_cicuit(self,jcc,rho=None,eps=1e-15):


        A,C,V=self.ACV;
        self.mV=V.tocsc();
        Nj=len(jcc);
        N=V.shape[0];
        Nf=N+Nj;  
        iAs=[A.col];    
        jAs=[A.row];
        dAs=[A.data[0:A.nnz]];

        iCs=[C.col];    
        jCs=[C.row];
        dCs=[C.data[0:C.nnz]];

        (dCs0,iCs0,jCs0)=[ copy.copy(c)  for c in   (dCs,iCs,jCs)];

        for k in range(Nj):
            ((iA,jA,dA),(iC,jC,dC),(iC0,jC0,dC0))=self.make_ecj(jcc[k],N+k,eps);

            iAs+=iA;
            jAs+=jA;
            dAs+=dA;

            iCs+=iC;
            jCs+=jC;
            dCs+=dC;

            iCs0+=iC0;
            jCs0+=jC0;
            dCs0+=dC0;

        

        Aex=_to_coo_scheme(Nf,row=iAs,col=jAs,data=dAs);
        Cex=_to_coo_scheme(Nf,row=iCs,col=jCs,data=dCs);
        Cex0=_to_coo_scheme(Nf,row=iCs0,col=jCs0,data=dCs0);

        
        (ua,ue,ur)=self.getunits()

        #bur=1./ur;
        #jcurr = lambda x: bur*x[:,N:];
        jcurr = lambda x: ur*x[:,N:];

        return (Aex,Cex,Cex0,jcurr,N);



class coil_t:
    def __init__(self,ACV,coil_mask,Nwind,cR_2pi):

        #cR_2pic=tolist(cR_2pic)
        #Nj=len(cR_2pic);
        A,C,V=ACV;
        
        mV=V.tocsc();

        coil_mask=np.array(coil_mask,dtype=int);

        N=coil_mask.size;
        N1=N+1;
        eta=np.array(coil_mask,dtype=np.double);
        eta=mV*eta;
        zm=Nwind/np.sum(eta);
        eta*=zm;
        f=np.abs(eta)>1.e-8;
        edata=eta[f];
        ii=np.arange(N);
        icol=ii[f];
        irow=N*np.ones(icol.size,dtype=icol.dtype);

        tcol=np.append(A.col,icol);
        trow=np.append(A.row,irow);

        
        nnz=C.nnz;

        tdata=np.append(C.data[0:nnz],edata);

        Cex=to_coo_scheme((N1,N1),row=trow,col=tcol,data=tdata);

        tcol,trow=np.append(trow,[N]),np.append(tcol,[N]);


        icol=np.append(icol,[N]);
        irow=np.append(irow,[N]);

        edata*=2;
        edata=np.append(edata,[cR_2pi]);
        
        tcol=np.append(C.col,irow);
        trow=np.append(C.row,icol);

        nnz=A.nnz;

        tdata=np.append(A.data[0:nnz],edata);

        Aex=coo_scheme((N1,N1),row=trow,col=tcol,data=tdata);

        self.AC=(Aex,Cex);

class tri_grad_t:
    def __init__(self,tri,vxs,index_base=0):

        pass


        

    
def test():
    '''
    A=coo_matrix([[1,2,3,4],[10,20,30,40],[100.,200,300,400],[1000,2000,3000,4000]])
    C=-A;
    V=coo_matrix([[1,0,0,0],[0,1,0,0],[0.,0,1.,0],[0,0,0,1]])
    cm=[0,1,0,0];
    R=777;

    coil=coil_t((A,C,V),cm,45,R);
    (a,c)=coil.AC;

    #(a,c)=coil_t((A,C,V),cm,88,R).AC;
    '''
    print('test');
    A=77*np.ones((7,7));
    C=33*np.ones((7,7));
    V=np.eye(7);

    b1=[0,1,0,0,0,0,0];
    b2=[0,0,1,0,0,0,0];
    b3=[0,0,0,0,0,1,0];
    jcc=[
        {
            'C':0.707,
            'R':11,'coils':[
            {'Nw':100,'mask':b3}
            ]
         },
            {'R':55,'coils':[
                {'Nw':10,'mask':b1},
                {'Nw':-10,'mask':b2}
                ]
             }
        ]

    sACV= [  sp.coo_matrix(m) for m in (A,C,V)]

    cc=coil_chain_t(sACV);

    a,c,c0,fuJ,nn=cc.make_cicuit(jcc);

    cc=c.tocsc();
    x=np.array([0,1,0,0,0,0,0,0,0]);
    y=cc*x;


    c=c.todense()
    c0=c0.todense()
    a=a.todense()

    print('A');
    print(a);
    print('C');
    print(c);
    print('C0');
    print(c0);
    pass

if __name__=='__main__':
    test()

