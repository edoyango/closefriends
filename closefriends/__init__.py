from .closefriends import closefriends as _closefriends
import numpy as _np


def celllist(x, cutoff, maxnpair):
    n, pairs = _closefriends.celllist_noreorder(x, cutoff, maxnpair)
    return pairs[0:2, :n]
