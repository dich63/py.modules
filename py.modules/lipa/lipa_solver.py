
import numpy as np
import ctypes
from scipy import sparse as sp
from scipy.sparse import csc_matrix,coo_matrix
from scipy.sparse import  linalg as sla




class lipa_node:

    def __init__(self,n,poles,res,A,C,x0):
        self.n=n;
        self.C=C;
        self.poles=poles;
        self.res=res;
        self.x0=x0;
        #bres=np.complex(1)/res;
        bres=1/res;
        Az=bres*(A+poles*C);
        options=dict(Equil=False, IterRefine='SINGLE')
        options=dict(Equil=False)
        self.lu=sla.splu(Az)#,options=options)

    def tic(self,x):
        return self.lu.solve(x);
    def itic(self):
        #for i in range(10):
        self.yy=self.lu.solve(self.x0);
        return True


class lu_polus:
    def __init__(self,Deltaz,Cz):
        #options=dict(Equil=False, IterRefine='SINGLE')
        #options=dict(Equil=False)
        self.lu=sla.splu(Deltaz)#,options=options)

    def step(self,xin,xout):
        xout[:]=self.lu.solve(xin);





class lu_poles_node:

    def zero_like(C):
        return csc_matrix(C.shape,dtype = C.dtype);
    def eye_like(C):
        return sp.identity(C.shape[0],dtype = C.dtype);

    def calc_poly(C,z):
        zn=1;
        cz=C[0];
        C=C[1:];
        for c in C:
            zn*=z;
            cz+=c*zn;
        return cz;

    def reset_J(z,J):

        if J:
            jj=zero_like(self.c0);
            z=self.polus;
            for d in J:
                g=d.get('g',0);
                m=-d.get('m',0);
                j=d['j'];
                j*=(z+g)**(m+1);
                jj+=j;
        else:
            jj=False;
        return jj;

    def __init__(self,AC,polus,res,xz,c0,cz0,freal=True,lu_solver_class=lu_polus):
        #AC=
        A=AC[0];
        C=AC[1:];
        z=polus;

        if len(C)==0:
            Cz=eye_like(C);
        else:
            Cz=calc_poly(C,z);

        br=1/res;
        Az =br*(A+z*Cz);
        self.Cz = Cz;
        self.polus=polus;
        self.J=np.zeros_like(cz0);
        (selt.xz,selt.c0,selt.cz0)=(xz,c0,cz0);
        #options=dict(Equil=False, IterRefine='SINGLE')
        #options=dict(Equil=False)
        self.lu=lu_solver_class(Az,Cz)#,options=options)
        self.dJ=False;
        """
        if freal:
            self.step=self._step_real;
        else:
            self.step=self._step_cmplx;
        """

    def step(self):
        xb=self.c0;
        self.lu.step(self,xb,xz);
        if not self.xz0 is self.cz0:
            self.cz0[:]=self.Cz*x;
        return True;


class lipa_solver_real_st:

    def __init__(self,n,pade,A,C):
        l=pade.count;


