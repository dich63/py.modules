import sys
#sys.path.append('v:/ipc/py.modules')

if __name__=='__main__':
    import os
    import numpy as np
    import matplotlib.pyplot as plt
    from matplotlib.pyplot import cm,bar,legend
    from midmeshgen.MID_decay_3B import make_barriers_3B,get_MID_decay_3B

    N=10
    color=iter(cm.rainbow(np.linspace(0,1,N)))
    color=iter(cm.spectral(np.linspace(0,1,N)))
    color=iter(0.9*cm.hot(np.linspace(0,1,N)))
    
    dt=0.001
    count=250
    th=np.linspace(1,20,N)
    for n in range(N):
        #th=0.2+float(n)/float(N)*20

        #b=make_barriers_3B(csg={"th":th[n]},sensor='LS')
        #b=make_barriers_3B(csg={"th":th[n]},sensor='LS')
        #b=make_barriers_3B(csg2={"th":th[n]},sensor='LS')
        
        #b["tbg"]["th"]=th[n];
        #b["csg"]["th"]=th[n];
        #b["csg2"]["th"]=th[n];
        smbg=(np.double(n>0)*0.000001*np.power(10,n) ,1)
        smbg=(np.double(n<N-1)*5000.0*np.power(10.0,-n) ,1)
        b=make_barriers_3B()
        #smbg=(0,1)
        si=b["tbg"]["sigma"];
        #si=(0.1+np.float(n)*si)
        #        si=0.01*si
        #
        si=2.1e6 
        
        #si=smbg[0]
        b["tbg"]["sigma"]=1*si;
        b["csg"]["sigma"]=1*si;
        b["csg2"]["sigma"]=si;
        mu=b["tbg"]["mu"]
        #
        mu=10.
        #mu=100.
        b["tbg"]["mu"]=mu;
        b["csg"]["mu"]=mu;
        b["csg2"]["mu"]=mu;

        b['zmax']=2500
        
        #fbc_mask=[n==0,1,1,1]
        #fbc_mask=[n==0,n==0,n==0,n==0]
        fbc_mask=[1,1,1,1]
        
        #kk=n-3
        #smbg=(np.power(10,kk) ,1)
        jst1=get_MID_decay_3B(b,count=count,dt=dt,fbc_mask=fbc_mask,sigma_mu_bg=smbg,mt=1,bag=0)

        print('end[',n+1,'] th=',th[n]);
        c=next(color)
        jj=np.clip(-np.array(jst1),0,np.Inf)
        jj[jj==0]=1;
        tt=dt*np.array(range((count+1)),dtype='d');
        #plt.plot(jj, c=c,label='sb=10^'+str(kk))
        #plt.plot(np.abs(jst1), c=c,label=str(th[n]))
        print(smbg)
        #plt.plot(jj, c=c,marker='.',ls='',label='{0:2.2g}'.format(smbg[0]))
        #
        plt.plot(tt,jj, c=c,label='{0:2.2g}'.format(smbg[0]))
        #plt.plot(jj, c=c,label=str(fbc_mask))
    plt.legend()
    plt.yscale('log')
    #    plt.xscale('log')
    plt.title('b')
    print('show...');
    #plt.show() 
    plt.savefig('zzz.svg')
    os.system('start zzz.svg')

"""
    #jst1=make_barriers_3B()
    jst1=get_MID_decay_3B()
    #plt.plot(np.abs(jst1), 'r-', np.abs(jst2), "-g")
    plt.plot(np.abs(jst1), 'r-')
    plt.yscale('log')
    plt.show()
"""

