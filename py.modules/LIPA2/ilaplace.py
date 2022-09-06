# -*- coding: utf-8 -*-
"""
Created on Tue Mar 29 15:35:29 2022

@author: wwww
"""
import sympy 
from sympy import *
from sympy.matrices import Matrix, eye, zeros, ones, diag, GramSchmidt
from sympy.abc import x,t
#from sympy.functions.special.hyper import *
from sympy.functions.special.hyper import *
from sympy.integrals.transforms import inverse_laplace_transform
from mpmath import *
import numpy as np
sexp=sympy.exp;

from jsobj import *

mp.dps = 55; 
mp.pretty = True

#mpmathify('1+4j',strings=True)

def isiter(o):
    try:
        i = iter(o)
        return True
    except TypeError:
        return False

def toiter(o):
    return o if isiter(o) else [o];
    
def issym(ex):
    return isinstance(ex, sympy.Expr)

def issymconst(ex,t):    
    return ex.subs(t,0)==ex;

'''
def ilaplace(F,z,plane=None):
    
    F,z=[sympify(s) for s in (F,z)]
    
    iF=inverse_laplace_transform(F,z,t,plane=plane);
    
    
    def fiF(r):
        if isiter(r):            
            return [ fiF(v) for v in r];
            
        iFt=iF.subs(t,r);
        return complex(iFt.evalf())
    
    return fiF,iF;
'''



def remequ(o):
    return list(dict.fromkeys(o))    

def repare_roots(*rs):
    rr=()
    for r in rs:
        if isiter(r):
            if len(r)>1:
                rr+=( ( r[0] ,int(r[1]) ), )
            else:
                rr+=((r[0],1),)
        else:
            rr+=((r,1),)          
            
    gg=remequ([r[0] for r in rr ])
    
    rrg=()
    
    for g in gg:
        rg=[g,0];
        for r in rr:
            if g==r[0]:
                rg[1]+=r[1]
        rrg+=( tuple(rg)  ,)
      
            
    return rrg;

def repare_qp(*rs):
    
    def addL(x,y):
        nx=len(x)
        ny=len(y);
        if nx<ny:
            m=nx
            x,y=y,x;
        else:
            m=ny;
            
        s=[v for v in x]        
        
        for n in range(m):
            s[n]+=y[n]
            
        return s;
    
    def tolist(q):
        return list(q) if isiter(q) else [q]
        
    rr=();
    for r in rs:
        if isiter(r):
            if len(r)>1:
                q,g=r;
                rr+=( (tolist(q),g ), )
            else:                
                rr+=( (tolist(r[0]),0),)
        else:
            rr+=((tolist(r),0),)
            
    gg=remequ([r[1] for r in rr ])
    
    rrg=()
    
    for g in gg:
        rg=[[0],g];
        for r in rr:
            if g==r[1]:
                rg[0]=addL(rg[0], r[0])
                
        rrg+=( tuple(rg) ,)
      
            
    return rrg;
    '''        
    return rr;
    '''

def make_qp(qp,z=Symbol('z')):
    qp=repare_qp(*qp);
    ex=0;
    rs=()
    for r in qp:
        q,g=r
        m=len(q)
        rs+=((g,m),)
        for k in range(m):
            ex=ex+q[k]/(z-g)**(k+1);
    return ex,rs,qp        

def coeff_series(ex,deg=6,z0=0,z=Symbol('z')):
    
    frl=sympy.factorial
    
    cc=[];
    ex=ex.subs(z,z+z0);
    ss=ex.series(z,0,n=deg+1)
    
    for m in range(deg+1):
        v=ss.diff(z,m).subs(z,0);
        v=v/frl(m);
        cc+=[v];          
    
    return cc,ss;

def prod_roots(rs,z=Symbol('z')):
    ex=1
    rs=repare_roots(*rs);
    deg=0;    
    
    for r in rs:        
        g,m=r;
        deg+=m;
        ex=ex*(z-g)**m;
        
    cc,ss=coeff_series(ex,deg,z0=0,z=z);    
    '''
    eex=expand(ex)
    
    ct=poly(eex).terms();
    cc=[None for k in range(len(ct))];
    
    for c in ct:
        i=int(c[0][0])
        cc[i]=c[1]
    '''
    
    return [ex,rs,cc,ss]


def res_expand(exp,roots,z=Symbol('z'),t=Symbol('t')):
    
    
    rr=()
    ex=0;
    roots=repare_roots(*roots)
    
    for root in roots:
        g,m=root
        Bs=[]
        zm=(z-g)**m;
        expzm=exp*zm;
        print('g=',g,'m=',m)
        for k in range(m):
           bk=expzm.diff(z,m-k-1).limit(z,g)/factorial(m-k-1);
           bk=complex(bk)
           Bs+=[bk];
        rr+=( (g,Bs),)
        
    return rr;
'''
def res(ex,m,k,z,z0):    
    exzm=simplify(ex*(z-z0)**m)
    return exzm.diff(z,m-k).limit(z,z0)/factorial(m-k)
'''
def res(ex,m,k,z,z0):    
    
    exzm=simplify(ex*(z-z0)**m)
    
    return exzm.diff(z,m-k).limit(z,z0)/factorial(m-k)

def resm(ex,m,z,z0):    
    exzm=simplify(ex*(z-z0)**m)
    bm=[]
    for k in range(m):        
        k1=k+1
        print(k1,m)
        r=exzm.limit(z,z0);
        print(r)
        bm+=[r]
        if k1<m:
            exzm=exzm.diff(z);
    return bm;
            
        
        
    return exzm.diff(z,m-k).limit(z,z0)/factorial(m-k)


def sym_res(ex,g,mr=1,z=Symbol('z')):      
    
    frl=sympy.factorial
    
    ex=ex.subs(z,z+g);
    exzm=ex*z**mr;
    ss=exzm.series(z,0,n=mr+1).removeO();
    bm=[];
    for m in range(mr):
        #print(mr-m,mr)
        b=ss.diff(z,m).subs(z,0)/frl(m);
        bm+=[b]
    bm.reverse();
    return bm;

def subs_keys(ex,o):
    
    '''
    if type(ex) is int:
        print('ex=',ex)
        return ex;
    '''    
    
    
    if type(ex) is str:
        ex=eval(ex);        
    elif len(o)==0:
        return ex;                
    elif isiter(ex):
        return tuple([ subs_keys(i,o)  for i in ex])
        
    for k in o:
        if hasattr(ex,'subs'):
            sk=Symbol(str(k))
            ex=ex.subs(sk,o[k]);
        else:
            break
            
    return ex;

def sym2complex(r):
    if isiter(r):
        return tuple([ sym2complex(i)  for i in r])
    else:
        return r if type(r) is int else complex(r);

def sym2numeric(ex,o):
    return sym2complex(subs_keys(ex,o))

def sym_repare(r):
    if type(r) is str:
        return Symbol(r)
    elif isiter(r):
        return tuple([ sym_repare(i)  for i in r])
    else:
        return r;    
        
    
def sym_reduce(DC,qp,z=Symbol('z')):
    
    DC,qp=[ sym_repare(v) for v in (DC,qp) ]
    
    exDC,rootsDC,DC_list,exDCss=prod_roots(DC,z=z)
    exF,rootsF,qp_list=make_qp(qp,z=z)
    
    exRF=exF/exDC;
    
    allroots=repare_roots(*(rootsDC+rootsF));
    
    res_poles=()
    
    for r,m in allroots:
        
        bm=sym_res(exRF,r,m,z);        
        res_poles+=((bm,r),);
    
    #qlmax=
        
    return res_poles,DC_list,qp_list;


def sym_reduce_ex(DC,qp=[],sq=[],z=Symbol('z')):
    def make_sq(sq,z):
        ex=0;       
        zn=1
        for c in sq:
            ex=ex+c*zn
            zn=zn*z;
        return ex;
            
            
    
    DC,qp,sq=[ sym_repare(v) for v in (DC,qp,sq) ]
    
    exDC,rootsDC,DC_list,exDCss=prod_roots(DC,z=z)
    exF,rootsF,qp_list=make_qp(qp,z=z)
    
    exF=exF+make_sq(sq,z)
    exRF=exF/exDC;
    
    allroots=repare_roots(*(rootsDC+rootsF));
    
    res_poles=()
    
    for r,m in allroots:
        
        bm=sym_res(exRF,r,m,z);        
        res_poles+=((bm,r),);
    
    #qlmax=
        
    return res_poles,DC_list,qp_list,sq;
        

def maxlen(qp):
    return max([len(q[0]) for q in qp ])

def make_shape(qp):
    return ( len(qp),maxlen(qp));


def make_numeric_data(d,q,o,dtype=np.complex128):
    
    [d,q]=sym2numeric([d,q],o)
    
    dn=np.array(d,dtype=dtype)
    
    if not q:
        return dn,np.array(q),np.array(q),None
    
    shape=make_shape(q);
    gn=np.zeros(shape[0],dtype=dtype)
    qpn=np.zeros(shape,dtype=dtype)
    FF=np.ones([1,shape[0]],dtype=dtype)
    for k in range(shape[0]):
        Q,gn[k]=q[k]
        Qp= qpn[k]
        for m in range(len(Q)):
            Qp[m]=Q[m]
    
    return dn,qpn,gn,FF;

def make_numeric_data_ex(d,q=[],sing=[],o={},dtype=np.complex128):    
    sing=sym2numeric(sing,o)
    dn,qpn,gn,FF=make_lipa_data(d,q,o,dtype);
    FFS=np.ones([len(sing),1],dtype=dtype);
    return dn,qpn,gn,sing,FF,FFS;

make_lipa_data=make_numeric_data
make_lipa_data_ex=make_numeric_data_ex

    
def syms_def(gl,cc):
    cc=cc.split(',');
    for c in cc:
        sc=str(c).strip()
        gl[sc]=Symbol(sc);
    pass

def toargs(d):
    r={};
    for k in d:
        r[str(k)]=d[k]
        
    return r;


def ilaplace_functor(rp,dp,nd=0):
    
    frl=sympy.factorial
    exp=sympy.exp
    
    rp=subs_keys(rp,dp)
    
    t=Symbol('t')
    
    ex=0;
    for q,r in rp:
        Q=0;
        tk=1
        for k in range(len(q)):
            Q=Q+q[k]/frl(k)*tk;
            tk=tk*t
        ex=ex+Q*exp(r*t);
        #print(r,':',ex)
    
    if nd:
        ex=ex.diff(t,nd)    
        
    ex=Heaviside(t,1)*ex;
    
    fu=lambdify(t,ex,'numpy') 
    
    def fua(x):
         return fu(np.array(x)) if isiter(x) else fu(x);
     
    def fuac(x):
        return np.full_like(x,fu(0),dtype=np.complex128) if isiter(x) else fu(x);
   
    return [fuac,ex] if issymconst(ex,t) else [fua,ex];          
    

def ilaplace_jet(rps,datas):
    
    
    frl=sympy.factorial
    exp=sympy.exp
    rp=subs_keys(rps,datas)
    
    t=Symbol('t')
    
    def getex(rp,nd=0): 
        ex=0;
        for q,r in rp:
            Q=0;
            tk=1
            for k in range(len(q)):
                Q=Q+q[k]/frl(k)*tk;
                tk=tk*t
            ex=ex+Q*exp(r*t);
        
        if nd>0:
            ex=ex.diff(t,nd)
        return ex;
    
    ex=getex(rp);
        
    def get_njet(nd):
        
        e=ex.diff(t,nd)        
        fu=lambdify(t,e,'numpy')         
        
        def fa(x):
             return fu(np.array(x)) if isiter(x) else fu(x);
         
        def fc(x):
            return np.full_like(x,fu(0),dtype=np.complex128) if isiter(x) else fu(x);
   
        return fc if issymconst(e,t) else fa;          
        
        
    
    class _inner_jet_t(object):
        def __init__(self):
            self.rp,self.rps,self.ex=rp,rps,ex
            self.getex=lambda nd=0 : getex(rps,nd);
            pass
        def __call__(self,ijet=0):
            
            if isiter(ijet):
                
                funcs=np.array([ get_njet(nd) for nd in ijet])
                nnd=funcs.shape[0];
                def ffjet(t):
                    return np.array([f(t) for f in funcs])                                
                return ffjet;
            else:
                return get_njet(ijet)
            
            
    return _inner_jet_t();


if __name__=='__main__':
    
    from utils import *
    #from jsobj import *
    '''
    c0=Symbol('c0');c1=Symbol('c1');c2=Symbol('c2');c3=Symbol('c3');c4=Symbol('c4');    
    b1=Symbol('b1');b2=Symbol('b2');b3=Symbol('b3');b4=Symbol('b4');b5=Symbol('b5');
    a1=Symbol('a1');a2=Symbol('a2');a3=Symbol('a3');a4=Symbol('a4');
    g1=Symbol('g1');g2=Symbol('g2');g3=Symbol('g3');
    
    y=Symbol('y');
    z=Symbol('z');
    g=Symbol('g');
    g1=Symbol('g1');
    
    G=(z-c1)*(z-c2)*(z-c3)*(z-c4)
    F=b1/(z-g)**1+b2/(z-g)**2+b3/(z-g)**3+b4/(z-g)**4+b5/(z-g)**5
    F1=a1/(z-g1)**1+a2/(z-g1)**2+a3/(z-g1)**3+a4/(z-g1)**4
    '''
    
    
    tic();
    #[p,d,q]=sym_reduce([c0,c1,c2,c3],[( [b1,b2,b3] ,g  )]);
    
    syms_def(locals(),' c0, c1 ,c0,c1,c2,c3,c4,b1, b2, b3,b4, g, a1, a2, a3, a4, g1')
    
    
    
    
    #[p,d,q]=sym_reduce([c0 ,c1,c0,c1,c2,c3],[( [0,b1,b2,b3] ,g  ),([0,0,0,a4],g1)]);
    
    #[p,d,q]=sym_reduce([(c0,2),(c1,2),c3,c3],[( [0,b1,b2,b3,b4] ,g  ),([0,0,0,a4],g1)]);
    [p,d,q]=sym_reduce([c0,c1,c2,c3],[( [b1,b2,b3,b4] ,g  )]);
    #[p,d,q]=sym_reduce([c0 ,c1,c0,c1,c2,c3],[( [0,b1,b2,b3] ,g  )]);
    
    toc('reduce:')
    #jo=arg2jso(c1=1.1,c2=2.3,c3=1.11,c4=21,c0=10j,g=3.11,g1=2j,b1=2.1,b2=1.2,b3=1.1j,a4=2)
    jo={
        c1:-0.015+7j,
        c2:-3.2,
        c3:-0.13+2j,
        c4:-21,
        c0:-0.01+2.1j,
        g:-1.051+7j,        
        b1:2.1e-3,
        b2:1.2,
        b3:10+1.1j,
        g1:-1+2j,
        a4:2.0e-2}
    
    jo={
    
    c0:-.5+250j,
    c1:-.5-265.5j,
    c2:-1.5+258j,
    c3:-1.5-268j,
    c4:-21,
    
    g:-57.5+260j,
    #g:0,                
    b1:10,
    b2:-1j,
    b3:3,
    b4:0,
    
    g1:-1+2j,
    a4:2.0e-12
    }


    
    
    '''
    tic();
    
    [pc,dc,qc]=sym2complex(subs_keys([p,d,q],toargs(jo)))
    toc('subs:')
    
    tic();
    
    GF=ilaplace_functor(pc,jo)
    F=ilaplace_functor(qc,jo)
    toc('functors:')
    '''
    
    tic();
    
    [GF,exg]=ilaplace_functor(p,jo)
    [F,exf]=ilaplace_functor(q,jo)
    toc('functors:')
    
    
    
    import matplotlib.pyplot as plt
    fig=plt.figure(figsize=(12,8))
    #fig.suptitle('sensors data + pure response ', fontsize=16)
    plt.grid(True,which='major')
    plt.minorticks_on()
    plt.grid(True,which='minor',alpha=0.2)
    #tt=np.linspace(0,300,2000)
    tt=np.linspace(0,10,7200)
    rn=lambda x: x.real/np.max(np.abs(x))
    #plt.plot(tt,rn(F(tt)),tt,rn(GF(tt)))
    plt.plot(tt,rn(F(tt)))
    plt.plot(tt,rn(GF(tt)))
    
    '''
    raise SystemExit()
    bb=sym_res(F/G,g,5)
    toc(':')
    raise SystemExit()
    
    F=1/(z-g)**4
    m=3;
    
    k=3
    tic()
    bm=[]
    for k in range(m):        
        print(m-k);
        bmk=res(F/G,m,m-k,z,g)       
        bm+=[bmk]
        toc(':')
    print(bm)
    
    
    raise SystemExit()
    
    
    ll=prod_roots((1,(9,1),4j,-4j,(7,1)));
    r=make_qp((1,[[2,3,4],3j],[[7],3j]));
    print('ll=',ll)
    print('=============')
    print('r=',r)
    
    F=r[0]/ll[3]
    rs=r[1]+ll[1]
    print('len=',len(rs))
    rp=res_expand(r[0],r[1]);rp
    
    rp=res_expand(F,rs,z)
    print(rp)
    
    
    
    a=(y**2+1*y+8)
    a1=(y**2+1*y+9)
    b=1/(a*a1)
    
    tic()
    #fu,iF=ilaplace(b/(y+3)**4,y)
    fu,iF=ilaplace(b*(2/(y+5)+1/(y+2)**2),y)
    toc(':')
    w=fu(1.2)
    '''
    
    
    
    
        
