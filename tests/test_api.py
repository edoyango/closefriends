def test_import():
    try:
        import closefriends
    except ImportError:
        assert False, "Failed to import closefriends"
    else:
        assert True


def test_returnarray():
    import closefriends
    from ancillary import setupf

    x, cutoff, maxnpair = setupf(3, 100, 60)
    pairs = closefriends.celllist(x, cutoff, maxnpair)
    pairshape = pairs.shape
    assert (len(pairshape) == 2) & (pairshape[0] == 2)
