# -*- coding: utf-8 -*-
"""
Created on Sat Jan 23 21:39:52 2016

@author: dich
"""
import numpy as np
from .jsonclass import config,class_decode,class_encode
from .ndarray_marshal import ndarray_base64_marshal,ndarray_base64_unmarshal,complex_json_marshal
from scipy.sparse import csc_matrix,coo_matrix


def sparse_base64_marshal(a):
    a=a.tocoo();
    col=ndarray_base64_marshal(a.col,'buffer',False)
    row=ndarray_base64_marshal(a.row,'buffer',False)
    data=ndarray_base64_marshal(a.data,'buffer',False)
    r={"data":data,"col":col,"row":row,"format":"COO","nnz":a.nnz,"lbound":0,"size":a.shape}
    return r

def sparse_json_marshal(a,json_class_name):
    t=a.dtype.type
    if (t==np.complex128) or (t==np.complex64):
        r=complex_json_marshal(a)
    else:
        r={'__jsonclass__':[json_class_name,[sparse_base64_marshal(a)]]}
    return r




def sparse_base64_unmarshal(o):
    fmt=o['format'].upper()
    if not fmt=='COO':
        raise Exception('sparse format not support');
    data=ndarray_base64_unmarshal(o['data'],'buffer')
    col=ndarray_base64_unmarshal(o['col'],'buffer')
    row=ndarray_base64_unmarshal(o['row'],'buffer')
    shape=o['size'];
    lb=np.int(o.get('lbound',1));
    if lb!=0:
        col=col-lb;
        row=row-lb;
    r=coo_matrix((data,(row,col)),shape=shape);
    return r





def sparse_json_unmarshal(d):
    #o=d['__jsonclass__'][1][0];
    return sparse_base64_unmarshal(d[1][0])
    
def register(cc):
    for c in cc:
        config.classes.register(c,'SparseMatrix',sparse_json_marshal,sparse_json_unmarshal);
    
       
register((coo_matrix,csc_matrix));

#config.classes.register(coo_matrix,'SparseMatrix',sparse_json_marshal,sparse_json_unmarshal)
#config.classes.register(csc_matrix,'SparseMatrix',sparse_json_marshal,sparse_json_unmarshal)