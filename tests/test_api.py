def test_import():
    try:
        import closefriends
    except ImportError:
        assert False, "Failed to import closefriends"
    else:
        assert True