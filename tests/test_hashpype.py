#/usr/bin/env python
#
# test_hashpype.py -MCW 2013
#

import unittest
from hashpy.hashpype import HashPype

def test_hashpype():
    hp = HashPype()
    s = hp.settings_str
    # TODO: check substring for defaults?
