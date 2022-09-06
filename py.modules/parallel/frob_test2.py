#import sys;sys.path.append('v:/ipc/py.modules')


import numpy as np
import env;

import  parallel.sparse as psp 



A=[
    [0,  2, 3,4 ],
    [5, 1., 7,8 ],
    [9, 10, 1,12],
    [13,14, 15,16]
    ]


nA=np.array(A)

sol=psp.sparse_solver_invoker(nA,opts={"common":{"numerical_rank":3}})
sol.factorize()

b=np.array([1,2,3,4+1j]);

x=np.dot(nA,b);



sol.y=x[:];
y=sol.y;
for k in range(1,2):
    sol.solve()
    print('[',k,'] y=',y)


t=11;