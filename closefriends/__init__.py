from .closefriends import closefriends as _closefriends
import numpy as _np

def celllist(x, cutoff, maxnpair):
    n, pair_i, pair_j =  _closefriends.celllist(x, cutoff, maxnpair)
    pairs = _np.stack((pair_i[:n], pair_j[:n]), axis=0)
    return _np.asfortranarray(pairs)
