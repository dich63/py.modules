import numpy as np
import scipy
#import matcompat

# if available import pylab (from matlibplot)
try:
    import matplotlib.pylab as plt
except ImportError:
    pass

def minmorph(reg, dmin):

    # Local Variables: a, k, flag, dmin, reg, dreg
    # Function calls: sort, false, length, minmorph, diff, remove_closest, true
    #% dmin = 0.25; %0.2 0.125
    reg = np.sort(reg)
    reg = remove_closest(reg, dmin)
    flag = True
    while flag:
        flag = False
        dreg = np.diff(reg)
        for k in np.arange(1., (len(dreg))+1):
            if dreg[int(k)-1] > dmin:
                flag = True
                a = (reg[int((k+1.))-1]-reg[int(k)-1])/2.+reg[int(k)-1]
                reg = np.array(np.hstack((reg[0:int(k)], a, reg[int(k+1.)-1:])))
                break
    return reg

def remove_closest(reg, dmin):

    # Local Variables: k, flag, dmin, reg, dreg
    # Function calls: false, length, abs, diff, remove_closest, true
    flag = True
    while flag:
        flag = False
        dreg = np.diff(reg)
        for k in np.arange(1., (len(dreg))+1):
            if np.abs(np.divide(dreg[int(k)-1], dmin))<1e-2:
                flag = True
                reg[int((k+1.))-1] = np.array([])
            break
    return reg