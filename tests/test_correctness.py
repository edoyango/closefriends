import ancillary, closefriends, scipy


def compare_ND(dim, maxPairsPerPoint):
    x, cutoff, maxnpair = ancillary.setupf(dim, 100, maxPairsPerPoint)

    """
    output is (100, dim) ndarray. To compare, need to:
        1. sort each pair so the lower index is first
        2. sort entire list by first index in each pair
    """

    tree = scipy.spatial.cKDTree(x.transpose())
    pairs_scipy = tree.query_pairs(cutoff, output_type="ndarray")
    pairs_scipy = sorted([sorted(xi) for xi in pairs_scipy])

    """
    output is (dim, 100) ndarray. To compare, need to:
        1. transpose for c-ordering
        2. sort each pair so the lower index is first
        3. sort entire list by first index in each pair
    """
    pairs_closefriends = closefriends.celllist(x, cutoff, maxnpair).transpose()
    pairs_closefriends = sorted([sorted(xi) for xi in pairs_closefriends])

    return pairs_closefriends == pairs_scipy


def test_1D():
    assert compare_ND(1, 4)


def test_2D():
    assert compare_ND(2, 13)


def test_3D():
    assert compare_ND(3, 34)


def test_5D():
    assert compare_ND(5, 166)
