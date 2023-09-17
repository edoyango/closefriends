def test_import():
    try:
        import closefriends
    except ImportError:
        assert False, "Failed to import closefriends"
    else:
        assert True


def test_returnarray():
    import closefriends
    from ancillary import setupc

    x, cutoff, maxnpair = setupc(3, 100, 60)
    pairs = closefriends.query_pairs(x, cutoff, maxnpair)
    pairshape = pairs.shape
    assert (len(pairshape) == 2) & (pairshape[1] == 2)
