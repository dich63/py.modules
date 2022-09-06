import numba as nb
import numpy as np

@nb.njit
def TDMA(a,b,c,d):
    n = len(d)
    x = np.zeros(n)
    bc = np.zeros(len(b))
    bc[0] = b[0]
    dc = np.zeros(len(d))
    dc[0] = d[0]
    
    for i in range(1, n):
        w = a[i - 1] / bc[i - 1]
        bc[i] = b[i] - w * c[i - 1]
        dc[i] = d[i] - w * dc[i - 1]

    x[n - 1] = dc[n - 1] / bc[n - 1]
    for k in range(n - 2, -1, -1):
        x[k] = (dc[k] - c[k] * x[k + 1]) / bc[k]
    return x