def test_import():
    try:
        import closefriends
    except ImportError:
        assert False, "Failed to import closefriends"
    else:
        assert True

# Checking that output_type="seperate-ndarrays" returns two seperate ndarrays
def test_returnseperatearrays():
    import closefriends
    from ancillary import setupc
    x, cutoff, maxnpair = setupc(3, 100, 60)
    pair_i, pair_j, idx = closefriends.query_pairs(x, cutoff, maxnpair, output_type="seperate-ndarrays")
    pair_i_shape = pair_i.shape
    pair_j_shape = pair_j.shape
    assert (len(pair_i_shape) == 1) and (len(pair_j_shape) == 1) and (len(idx.shape) == 1)

# Checking that output_type="ndarray" returns a 2D ndarray
def test_returnarray():
    import closefriends
    from ancillary import setupc

    x, cutoff, maxnpair = setupc(3, 100, 60)
    pairs, idx = closefriends.query_pairs(x, cutoff, maxnpair, output_type="ndarray")
    pairshape = pairs.shape
    assert (len(pairshape) == 2) and (pairshape[1] == 2) and (len(idx.shape) == 1)

# Checking that specifying output_type="set" returns a set of tuples
def test_returnset():
    import closefriends
    from ancillary import setupc

    x, cutoff, maxnpair = setupc(3, 100, 60)
    pairs, idx = closefriends.query_pairs(x, cutoff, maxnpair, output_type="set")
    
    istuples = [isinstance(ipair, tuple) for ipair in pairs]

    assert isinstance(pairs, set) and all(istuples)

# Testing scipy compatibility api, which corresponds to retain_order=True.
# This should ensure the x array order is unmodified, as well as indices in each pair being ordered such that i < j
def test_retainorder():
    import closefriends
    from ancillary import setupc

    x, cutoff, maxnpair = setupc(3, 100, 60)
    x_before = x.copy()
    pairs = closefriends.query_pairs(x, cutoff, maxnpair, retain_order = True, output_type="ndarray")

    isiltj = [i < j for i, j in pairs]

    assert all(isiltj) and (x == x_before).all()

# Check that the default behaviour is to return an ndarray where x and idx is returned
def test_defaultreturn():
    import closefriends
    from ancillary import setupc

    x, cutoff, maxnpair = setupc(3, 100, 60)

    x_before = x.copy()

    pairs_idx = closefriends.query_pairs(x, cutoff, maxnpair)

    returnstuple = isinstance(pairs_idx, tuple)

    pair_i = pairs_idx[0]
    pair_j = pairs_idx[1]
    idx = pairs_idx[2]
    
    pair_i_shape = pair_i.shape
    pair_j_shape = pair_j.shape
    pairs_is_seperatendarrays = (len(pair_i_shape) == 1) and (len(pair_j_shape) == 1)
    idx_is_1dndarray = (len(idx.shape) == 1)

    x_is_modified = (x != x_before).any()

    any_i_gt_j = any(pair_i>pair_j)

    assert returnstuple and pairs_is_seperatendarrays and idx_is_1dndarray and x_is_modified and any_i_gt_j