def setupf(dim, npoints, maxnpair_perp):
    import numpy as np
    x = np.random.rand(dim*npoints).reshape((dim, npoints), order='F')
    cutoff = npoints**(-1./dim)
    maxnpair = maxnpair_perp*npoints
    return x, cutoff, maxnpair
