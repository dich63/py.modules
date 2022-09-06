# -*- coding: utf-8 -*-
"""
Created on Fri Apr  1 22:50:40 2022

@author: wwww
"""

from LIPA2.ilaplace import *
from LIPA2.qp_solver import *

from utils import *
from LIPA2.qp import *

norm=np.linalg.norm

'''
import matplotlib.pyplot as plt
fig=plt.figure(figsize=(12,8))
#fig.suptitle('sensors data + pure response ', fontsize=16)
plt.grid(True,which='major')
plt.minorticks_on()
plt.grid(True,which='minor',alpha=0.2)
'''


def ilaplace_jet(rps,datas={}):    
    
    frl=sympy.factorial
    exp=sympy.exp
    rp=rps
    #rp=subs_keys(rps,datas)
    
    class jet(object):
        def __call__(self):
            return rp
    
        
    '''
    
    '''
    
    return jet();
        
        




FF=np.ones([1,1],dtype=complex)
G=np.ones([1,1],dtype=complex)

qp=np.array([[1,1,1,1]],dtype=complex)

QP=np.empty_like(qp,dtype=complex)


QP[:]=qp;

g=-.0
G[:]=g

print(QP.shape)
print(G.shape)
print(FF.shape)




DC =[0, 1]
LM=[6,8]


dt=1

lqp=lipa_qp_number(DC,LM=LM,FF=FF,qp=QP,g=G).reset(dt);

nr=111

y=lqp(nr)

t=nr*dt
#print((1-np.exp(-g*t*nr))/g)
ye=np.sum(qp*[t,t**2/2,t**3/6,t**4/24])
#ye=(np.exp(g*t)-1)/g

#ye=np.sum(qp*[t,t**2/2])
print('ex=',ye)
print('sol=',y)
print('err=',(y-ye))
print('rerr=',100.0*norm((y-ye).reshape(-1))/norm((ye).reshape(-1)),'%')


#print(np.sum(qp*[t,t**2/2]))
