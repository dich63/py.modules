# -*- coding: utf-8 -*-
"""
Created on Sat May  7 02:05:08 2022

@author: wwww
"""

import lipa.pade_exp_poles_res
from utils import *
from utils.c_structs import *
#from LIPA2.qp_solver import *
from LIPA2.tools import *
#from klu_lipa import *
from jsonrpc.json_io import *
import os
#import scipy
#import scipy.sparse as sp
    
import sparse_batch.sp_jet as sj
from sparse_batch.sp_jet import *



norm=np.linalg.norm
normm= lambda x:  norm(x,ord=np.inf)
randnc= lambda *lp:  np.random.randn(*lp)+1j*np.random.randn(*lp)

def sp_mulsum(CC,xx,yout):
        yout[:]=0;
        for k in range(len(CC)):
            yout+=CC[k].dot(xx[k])
        return yout;

N=4;

A=sp.eye(N,dtype=complex)
B=sp.eye(N,dtype=complex)
C=sp.eye(N,dtype=complex)

w=A+B
fn=r'O:\__ss\matrix\sfFEM1024k.json'
#fn=r'O:\__ss\matrix\sfFEM128k.json'
#fn=r'O:\__ss\matrix\sfFEM24k.json'
d=decode(fn,1)

A,B,C=[sp.csr_matrix(m,dtype=complex) for m in [d.K,d.G,d.M]]

N=A.shape[0];

z=0+0j

Jet=(A,B,C)
Jet=(A,B,B,C)
D=len(Jet)
xx=randnc(D,N);
xx[:]=1+0j


for k in range(D):
    xx[k,:]=range(1,N+1);

y=randnc(N);
y2=randnc(N);
y2[:]=np.nan
pr=lipa.pade_exp_poles_res.poles_res(8,8)
ml=10;
#zz=np.array([ml*p[0] for p in pr])
zz=np.array([p[0] for p in pr])
res=np.array([p[1] for p in pr])
zzD=np.array([zz**k for k in range(1,D+1)])


zzDt=np.ascontiguousarray(zzD)
sm=make_spm(*Jet);
tic()
smzz=new_spm_like(sm,nbatch=len(zz))
zzDt=np.array([zz**k for k in range(1,D+1)])
ppzz=ptr_ptr_array(zzDt);
toc('spm_like')
print('pid=',os.getpid())
tic();sj._csr_mul_zz(sm.ptr,ppzz.ptr,smzz.ptr,0,1111);toc('_csr_mul_zz')
#def make_crs_smzz(sm,z,chunk=11111,nthread=0,smzz_out=None,dtype=np.complex128)
tic();
smz2=None
for k in range(1):
    #smz2=None
    smz2=make_crs_smzz(sm,zz,1111,0,smz2);
toc('make_crs_smzz')
    
raise SystemExit(0);

ppx=ptr_ptr_array(xx);z
py=ptr_array(y2);
py=y2.ctypes.data

print('pid=',os.getpid())
    
tic(); sp_mulsum(Jet,xx,y);t1=toc('sp_mulsum:')   


nt,chunk=0,290
tic(); st=sj._csr_mulsum(sm.ptr,ppx.ptr,py,nt,chunk);t2=toc('csr_mulsum:')

print("nt,chunk=",nt,chunk)
err=2*norm(y-y2)/(norm(y)+norm(y2))
print('err',err)
print('perf',t1/t2)

xz=xxz_buffer_t(pr,xx)

'''
print(xz.xxz.real)
xz()
print(xz.xxz.real)
'''
tic();xz();toc('jetz')

sh=(len(zz),)+xx.shape
xxz_out=np.empty(sh,dtype=complex);
yyz_out=np.zeros((len(zz),xx.shape[1]),dtype=complex);
yyz_outp=np.zeros((len(zz),xx.shape[1]),dtype=complex);

ppyz=ptr_ptr_array(yyz_out);

tic();st=xz.gaxpyjet(sm,ppyz,1111);toc('gaxpyjet')
yyz_out[:]=0
tic();st=xz.gaxpyjet(sm,ppyz,1111);t1=toc('gaxpyjet')
#print('xz.xxz.real',xz.xxz.real)
#print('np.sum(xz.xxz,axis=1)')
#print(np.sum(xz.xxz,axis=1))
tic()
for k in range(xz.xxz.shape[0]):
    sp_mulsum(Jet[1:],xz.xxz[k],yyz_outp[k]);
t2=toc('gaxpyjet_python')    
    
err=2*norm(yyz_outp-yyz_out)/(norm(yyz_outp)+norm(yyz_out))
print('err',err)
print('perf',t2/t1)    

xxs=copy.deepcopy(xxz_out);
ppxxs=ptr_ptr_array(xxs);

print('pid=',os.getpid())
tic();st=xz.assembly(ppxxs.ptr,1111,0);toc('xz.assembly')



'''
print('yyz_out=')
print(yyz_out.real)
    
print('yyz_outp=')
print(yyz_outp.real)
'''

'''
tic()
JetA=(222*A,B,C)
AzC0(xx,JetA,yyz_outp);
for k in range(len(zz)):
    xxz_out[k][0]=xx[0]
    AzC1(xx,JetA,zz[k],xxz_out[k],yyz_outp);
    
toc('jet:')


print('yp=',yyz_outp)


'''






    