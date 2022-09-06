# -*- coding: utf-8 -*-
"""
Created on Mon Apr  5 19:56:39 2021

@author: wwww
"""

import csv
import numpy as np;



_parsecmlx=lambda r: eval(r.replace('i','j'))


def _cmlx2str(v):
    s=str(v.real);
    vi=v.imag;
    if vi>=0:
        s+="+";
    s+=str(vi)+"i";
    return s
    

def csv2array(fn,with_header=False):    
    
    
    ls=[]    
    with open(fn) as cf:
        sr=csv.reader(cf,delimiter=',')
        head=next(sr)
        for r in sr:
            pp=[_parsecmlx(s) for s in r ]
            ls+=[pp]
        
    
    f=np.array(ls);
    return (f,head) if with_header else f;


def strip_dict0(s):
    d={};
    for k,v in s.items():
        d[k.strip()]=s[k];
        
    return d;



def csv2dict(fn,farray=True,ftrimkey=True):
    
    def strip_dict(s):
        d={};
        for k,v in s.items():
            name=k.split(';')[0];
            d[name.strip()]=s[k];       
            
        return d;

    d={};
    with open(fn) as cf:
        sr=csv.DictReader(cf,delimiter=',')        
        names=sr.fieldnames;
        #names=[s.strip() for s in sr.fieldnames ];
        for n in names:            
            d[n]=[];
        
        for r in sr:
            for k,v in r.items():
                d[k]+=[_parsecmlx(v)]
            
    if farray:
        for k,v in d.items():
            d[k]=np.array(d[k]);
    
    return strip_dict(d) if ftrimkey else d;

def dict2csv(d,fn):
    
    kn=[k for k in d.keys()]
    skn=[str(k) for k in kn];
    
    ds={};
    for k in kn:
        ds[k]=np.array(d[k]).reshape(-1);
        
    #x=ds[kn[0]];
    #N=np.array(x).size;
    
    N=ds[kn[0]].size;
    
    with open(fn,'w',newline='') as cf:      
        
        writer = csv.DictWriter(cf, fieldnames=skn,delimiter=',')
        writer.writeheader();
        
        for n in range(N):
            sdn={}
            for k in kn:
                sdn[str(k)]=_cmlx2str(ds[k][n])
                
            writer.writerow(sdn);
            
        

def jso2csv(o,fn):
    from  jsobj import to_dict
    dict2csv(to_dict(o),fn);
        
def csv2jso(fn,farray=True):
    from  jsobj import jso;
    return jso(csv2dict(fn,farray));
    
    


if __name__=='__main__':
    pass
    #with open("V:/HW/dSigXY.csv") as cf:    
    d=csv2dict("V:/HW/d.csv")
    d=csv2jso("V:/HW/dSigXY.csv")
    jso2csv(d,"V:/HW/d2.csv")
    #csv2array("V:/HW/d.csv")