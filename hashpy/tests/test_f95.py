
import hashpy
#from hashpy.io.antelopeIO import load_pf
# 
def test_f95():
    assert True
    script = '''
    hp = hashpy.HashPype()
    load_pf(hp, 'dbhash')
    hp.input('/data/2013/005/reno', 'ANTELOPE', orid=979567)
    hp.driver2(False, False)
    print hp.output()
    '''


