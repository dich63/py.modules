from  ltx.mm  import mm_buffer
import numpy as np
from scipy import sparse as sp
from scipy.sparse import csc_matrix,coo_matrix


def ltx2sparse(om):
    if om.format!='COO':
        raise Exception('Bad ltx sparse format!');
    size=om.size;    
    shape=(size[0],size[1])
    row=mm_buffer(om.row);
    col=mm_buffer(om.col);
    data=mm_buffer(om.data);
    r=np.array(row.value)
    c=np.array(col.value);
    d=np.array(data.value);

    try:
        l=om.lbound;
    except Exception:
        l=1

    if l>0:
        c-=l
        r-=l

    return coo_matrix((d,(r,c)),shape=shape,copy=True);


