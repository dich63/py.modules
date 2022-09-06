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

def tolist(l):
    if (type(l)==tuple) or (type(l)==list):
        return l;
    else:
        return (l,);

def nzmask(a,eps=1e-15):    
    return np.abs(a)>eps;

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
            mCR=+cc*R;

            dCCs=[mCR*dA for dA in dAs]            
            dCs+=dCCs;
            iCs+=iAs;
            jCs+=jAs;



        iAs+=[anj];
        jAs+=[anj];
        dAs+=[aR];

        return ((iAs,jAs,dAs),(iCs,jCs,dCs))


        

    def make_cicuit(self,jcc,eps=1e-15):

        def to_coo_scheme(Nf,col,row,data):
            col=np.concatenate(col);
            row=np.concatenate(row);
            data=np.concatenate(data);
            return coo_scheme((Nf,Nf),row=row,col=col,data=data);

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


        for k in range(Nj):
            ((iA,jA,dA),(iC,jC,dC))=self.make_ecj(jcc[k],N+k,eps);

            iAs+=iA;
            jAs+=jA;
            dAs+=dA;

            iCs+=iC;
            jCs+=jC;
            dCs+=dC;

        

        Aex=to_coo_scheme(Nf,row=iAs,col=jAs,data=dAs);
        Cex=to_coo_scheme(Nf,row=iCs,col=jCs,data=dCs);

        
        (ua,ue,ur)=self.getunits()

        #bur=1./ur;
        #jcurr = lambda x: bur*x[:,N:];
        jcurr = lambda x: ur*x[:,N:];

        return (Aex,Cex,jcurr,N);



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
            'C':0.7070707,
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

    a,c,fuJ,nn=cc.make_cicuit(jcc);


    c=c.todense()
    a=a.todense()

    print('A');
    print(a);
    print('C');
    print(c);
    pass

if __name__=='__main__':
    test()

