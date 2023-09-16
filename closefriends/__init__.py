from .closefriends import closefriends as _closefriends

def query_pairs(x, cutoff, maxnpair):
    n, pairs = _closefriends.query_pairs_noreorder(x, cutoff, maxnpair)
    return pairs[0:2, :n]
