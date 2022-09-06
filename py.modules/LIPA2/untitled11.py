# -*- coding: utf-8 -*-
"""
Created on Wed Apr 20 21:45:20 2022

@author: wwww
"""
import numpy as np
import pypardiso
from utils import *
from parallel.sparse import *
from jsonrpc.json_io import *
from utils import *
from LIPA2.qp_solver import *
from scipy.sparse.linalg import splu

norm=np.linalg.norm

import scipy.sparse as sp


def sp_LU_f(A):    
    #lu=sp.linalg.splu(A
    options=dict(Equil=False,PivotGrowth=True
                                     ,PrintStat=False
                                     ,DiagPivotThresh=0.1
                                     ,ColPerm='MMD AT PLUS A'#'MMD_ATA'#
                                     ,ConditionNumber=False,IterRefine='NOREFINE')
    #    options={}
    #options=dict(Equil=False, IterRefine='SINGLE')
    lu=sp.linalg.splu(A,permc_spec='MMD AT PLUS A',options=options)
    
    return lambda x: lu.solve(x) ;


N=100#00
density=0.0001
tic()
A = sp.eye(N)+sp.rand(N, N, density=density, format='csr')+1j*sp.rand(N, N, density=density, format='csr')
toc('cr:')
b =np.random.rand(N) + 1j*np.random.rand(N)

cc=A.tocsc()


tic()
y = sp_LU_f(A)(b)
toc('scipy:')
err=2*norm(A*y-b)/(norm(b)+norm(A*y))
print('err=',err)


tic()
si=sparse_solver_invoker(A)
st=si.analyze()
#toc('analy:')
#tic()
st=si.factorize()
#toc('lu ctx')
si.y=b; 
si.solve()
y=si.y
toc('lipa ctx ')
err=2*norm(A*y-b)/(norm(b)+norm(A*y))
print('err=',err)





#pp=pypardiso.PyPardisoSolver()
tic()
y = pypardiso.spsolve(A, b)
toc('pypardiso:')
err=2*norm(A*y-b)/(norm(b)+norm(A*y))
print('err=',err)
