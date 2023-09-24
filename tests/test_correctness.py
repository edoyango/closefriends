import ancillary, closefriends
from scipy.spatial import cKDTree


def compare_ND_set(dim, maxPairsPerPoint):
    x, cutoff, maxnpair = ancillary.setupc(dim, 100, maxPairsPerPoint)

    # output is set of tuples.
    tree = cKDTree(x)
    pairs_scipy = tree.query_pairs(cutoff, output_type="set")

    # output is set of tuples.
    pairs_closefriends = closefriends.query_pairs(x, cutoff, maxnpair, output_type="set", retain_order=True)
    
    # pairs_closefriends = sorted([sorted(xi) for xi in pairs_closefriends])
    return pairs_closefriends == pairs_scipy

def compare_ND_ndarray(dim, maxPairsPerPoint):
    x, cutoff, maxnpair = ancillary.setupc(dim, 100, maxPairsPerPoint)

    """
    output is (100, dim) ndarray. To compare, need to:
        1. sort entire list by first index in each pair
    """

    tree = cKDTree(x)
    pairs_scipy = tree.query_pairs(cutoff, output_type="set")
    pairs_scipy = sorted([sorted(list(xi)) for xi in pairs_scipy])
    print(pairs_scipy)

    """
    output is (dim, 100) ndarray. To compare, need to:
        1. return to original ordering
        2. sort each pair so the lower index is first
        3. sort entire list by first index in each pair
    """
    pairs_closefriends, idx = closefriends.query_pairs(x, cutoff, maxnpair, output_type="ndarray", retain_order=False)
    
    pairs_closefriends = [[idx[i], idx[j]] for i, j in pairs_closefriends] # return to original ordering
    pairs_closefriends = [sorted(pairi) for pairi in pairs_closefriends] # sort each pair so lower index is first
    pairs_closefriends = sorted(pairs_closefriends) # sort pairs

    return pairs_closefriends == pairs_scipy

def test_1D():
    assert compare_ND_set(1, 4) and compare_ND_ndarray(1, 4)


def test_2D():
    assert compare_ND_set(2, 13) and compare_ND_ndarray(2, 13)


def test_3D():
    assert compare_ND_set(3, 34) and compare_ND_ndarray(3, 34)


def test_5D():
    assert compare_ND_set(5, 166) and compare_ND_ndarray(5, 166)
