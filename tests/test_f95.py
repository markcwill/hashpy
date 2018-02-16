#
import hashpy

 
# Import will fail if build/linking to Fortran ext is bad
def test_f95():
    assert True
    script = '''
    hp = hashpy.HashPype()
    load_pf(hp, 'dbhash')
    hp.input('/data/2013/005/reno', 'ANTELOPE', orid=979567)
    hp.driver2(False, False)
    print hp.output()
    '''


