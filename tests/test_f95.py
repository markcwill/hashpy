#
"""
"""
 
# Import will fail if build/linking to Fortran ext is bad
def test_f95():
    import hashpy.libhashpy
    assert True

