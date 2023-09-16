import ancillary, closefriends, scipy, numpy as np

def test_1D():
    x, cutoff, maxnpair = ancillary.setupf(1, 100, 4)

    '''
    output is (100, 2) ndarray. To compare, need to:
        1. sort each pair so the lower index is first
        2. sort entire list by first index in each pair
    '''
    tree = scipy.spatial.cKDTree(x.transpose())
    pairs_scipy = tree.query_pairs(cutoff, output_type="ndarray")
    pairs_scipy = sorted([sorted(xi) for xi in pairs_scipy])
    
    '''
    output is (1, 100) ndarray. To compare, need to:
        1. transpose for c-ordering
        2. sort each pair so the lower index is first
        3. sort entire list by first index in each pair
    '''
    pairs_closefriends = closefriends.celllist(x, cutoff, maxnpair).transpose()
    pairs_closefriends = sorted([sorted(xi) for xi in pairs_closefriends])
    
    assert pairs_closefriends == pairs_scipy
