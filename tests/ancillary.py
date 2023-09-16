import numpy as np

def setupf(dim, npoints, maxnpair_perp):
    x = np.random.rand(dim*npoints).reshape((dim, npoints), order='F')
    cutoff = 2.*npoints**(-1./dim)
    maxnpair = maxnpair_perp*npoints
    return x, cutoff, maxnpair