import numpy as np


def setupc(dim, npoints, maxnpair_perp):
    rng = np.random.default_rng(12345)
    x = rng.random(npoints*dim, dtype=np.float64).reshape((npoints, dim))
    cutoff = 2.0 * npoints ** (-1.0 / dim)
    maxnpair = maxnpair_perp * npoints
    return x, cutoff, maxnpair
