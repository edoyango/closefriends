import numpy as np


def setupc(dim, npoints, maxnpair_perp):
    x = np.random.rand(npoints, dim)
    cutoff = 2.0 * npoints ** (-1.0 / dim)
    maxnpair = maxnpair_perp * npoints
    return x, cutoff, maxnpair
