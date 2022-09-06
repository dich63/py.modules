import numpy as np


def minmorph(reg, dmin):
    reg = np.sort(np.array(reg).flatten())
    reg = remove_closest(reg, dmin)

    out = np.array([])
    for i in range(len(reg) - 1):
        if reg[i+1] - reg[i] > 2 * dmin:
            steps = np.round((reg[i+1] - reg[i]) / dmin)
            out = np.r_[out, np.linspace(reg[i], reg[i+1], int(steps))[:-1]]
        else:
            out = np.r_[out, reg[i]]
    out = np.r_[out, reg[-1]]
    return out


def remove_closest(reg, dmin):
    regdf = np.diff(reg)
    ind = []

    for i in range(len(regdf)):
        if regdf[i] < 1e-2:
            ind.append(i+1)
    out = np.delete(np.array(reg), ind)
    return out.tolist()
