import numpy as np


def setupf(dim, npoints, maxnpair_perp):
    x = np.random.rand(dim * npoints).reshape((dim, npoints), order="F")
    cutoff = 2.0 * npoints ** (-1.0 / dim)
    maxnpair = maxnpair_perp * npoints
    return x, cutoff, maxnpair
