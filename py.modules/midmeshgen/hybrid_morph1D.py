import numpy as np

def hybrid_morph1D_(seg=None,q=None,*args,**kwargs):

    morphmesh=seg[0]

    for k in np.arange(0,(len(seg) - 1)):
        dx1=seg[k][-1] - seg[k][-2]
        l = seg[k + 1][0] - seg[k][-1]
        if len(seg[k + 1]) > 1:
            dx2=seg[k + 1][1] - seg[k + 1][0]
            xgv0=morph_(dx1,dx2,l,q)
        else:
            xgv0=xgv_(dx1,l,q,1)
            xgv0=np.delete(xgv0, 0)
        xgv0=xgv0 + seg[k][-1]
        morphmesh=np.append(morphmesh,xgv0)
        morphmesh = np.append(morphmesh,seg[k + 1])
    morphmesh=np.delete(morphmesh, -1)
    return morphmesh

def morph_(dx1=None,dx2=None,l=None,q=None,*args,**kwargs):

    if (l <= 0):
        print ('geometry error')
        return
    lmin=(1 + q + q ** 2 + q ** 3) * (dx1 + dx2)
    if (l > lmin):
        xgv0=exp_morph_(dx1,dx2,l,q)
    else:
        if (l > (dx1 + dx2)):
            if ((l - dx1 - dx2) > (dx1 + dx2)):
                xgv0=np.append(dx1,prb_morph_(dx1,dx2,l))
                xgv0=np.append(xgv0,(l - dx2))
            else:
                xgv0=np.append(dx1, l - dx2)
                if ((l - dx1 - dx2) < 0.3 * np.max([dx1,dx2])):
                    pass
        else:
            if (l < 0.4 * np.min([dx1,dx2])):
                print ('morphing is not possible')
            xgv0=np.array([])
    return xgv0

def exp_morph_(dx1=None,dx2=None,l=None,q=None,*args,**kwargs):

    k=np.log(dx1 / dx2) / np.log(q)
    k=np.round(k)
    c=l * (q - 1) + dx2 + dx1
    f=c / (dx1 + dx2 * q ** k)
    n1 = np.log(f) / np.log(q)
    n1=np.floor(n1)
    l1=dx1 * (q ** n1 - 1) / (q - 1)
    l2=l - l1
    if (l1 < (1 + q) * dx1):
        l1=np.copy(dx1)
        l2=l - l1
    if (l2 < (1 + q) * dx2):
        l2=np.copy(dx2)
        l1=l - l2
    xgv1=xgv_(dx1,l1,q,1)
    xgv2=xgv_(dx2,l2,q,- 1)[0]
    xgv0=np.append(np.array(xgv1),xgv2[1:len(xgv2)] + l1)
    xgv0=np.delete(xgv0, 0)
    xgv0=np.delete(xgv0, -1)
    return xgv0

def xgv_(dx=None,l=None,q0=None,flag=None,*args,**kwargs):

    if (l > dx):
        a=l * (q0 - 1) / dx + 1
        n0=np.log(a) / np.log(q0)
        n=int(np.floor(n0))
        k=- l
        m=l - dx
        h=np.array(dx)
        if (n > 2):
            p=np.append(h,np.zeros(n - 2))
            p = np.append(p, k)
            p = np.append(p,m)
        else:
            p=np.append(h,k)
            p = np.append(p,m)
        r=np.roots(p)
        f=np.copy(True)
        while f:

            c,mq= (np.abs(r.real - q0)).min(),(np.abs(r.real - q0)).argmin()
            #q=float(r[mq])
            q=r[mq].real
            s=dx * (q ** n - 1) / (q - 1)
            if (np.abs(s - l) / l) < 1e-06:
                f=np.copy(False)
            else:
                r = np.delete(r,mq)

        xgv=np.zeros(n + 1)
        for i in np.arange(1,(n + 1)):
            xgv[i]=xgv[i - 1] + dx * q ** (i - 1)
    else:
        xgv=np.append(0,dx)
    if (flag == - 1):
        xgv=np.fliplr([xgv])
        xgv=np.abs(xgv - l)
    return xgv

def prb_morph_(dx1=None,dx2=None,l=None,*args,**kwargs):
    xgv0=morphgrid1_(dx1,l - dx2,dx1,dx2)
    return xgv0

def morphgrid1_(L=None,R=None,hL=None,hR=None,*args,**kwargs):

    dX=R - L
    n,q=stepcalc_(dX,hL,hR,nargout=2)
    n=int(n)
    G=np.zeros(n - 1)
    g=np.copy(L)
    for k in np.arange(n - 1):
        h=q[0] * k ** 2 + q[1] * k + q[2]
        g=g + h
        G[k]=g
    return G

def stepcalc_(dX=None,hL=None,hR=None,*args,**kwargs):

    n=2 * dX / (hR + hL)
    n=np.floor(n)
    q=step_co_(n,dX,hL,hR)
    return n,q

def step_co_(n=None,dX=None,hL=None,hR=None,*args,**kwargs):

    nn=np.arange(0,n)
    sn1=np.sum(nn)
    sn2=np.sum(nn ** 2)
    A=np.array([[sn2,sn1,n],[(n + 1) ** 2,(n + 1),1],[0,0,1]])
    q=np.linalg.solve(A,[[dX],[hR],[hL]])
    return q