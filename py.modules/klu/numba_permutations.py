#
import numba  
import numba as nb
import numpy as np


@numba.jit("int32(int32,complex128[:],complex128[:],int32[:])"
           ,nopython=True,
           parallel=True,
           nogil=True,
           fastmath=True           
           )
def permute(n,s,d,ip):
    for i in numba.prange(n):
        d[i]=s[ip[i]];
    return n;


@numba.jit("int32(int32,int32,complex128[:,:],complex128[:,:],int32[:])"
           ,nopython=True,
           parallel=True,
           nogil=True,
           fastmath=True           
           )
def permute2(m,n,s,d,ip):
    for i in numba.prange(n):
        for k in range(m):
            d[k,i]=s[k,ip[i]];
    return n;

@numba.jit("int32(int32,int32,complex128[:,:],complex128[:,:],int32[:,:])"
           ,nopython=True,
           parallel=True,
           nogil=True,
           fastmath=True           
           )
def permute3(m,n,s,d,ip):
    for i in numba.prange(n):
        for k in range(m):
            d[k,i]=s[k,ip[k,i]];
    return n;

@numba.jit("int32(int32,int32,complex128[:,:],complex128[:,:],int32[:,:],float64[:,:])"
           ,nopython=True,
           parallel=True,
           nogil=True,
           fastmath=True           
           )
def permute3_rescale(m,n,s,d,ip,rs):
    for i in numba.prange(n):
        for k in range(m):
            j=ip[k,i]
            d[k,i]=s[k,j]/rs[k,j];
    return n;

