#!/usr/bin/env python

from distutils.core import setup
import os

kwargs = {'name' : 'HASHpy',
          'version' : '1.0',
          'description' : 'Routines for running HASH algorithms',
          'author' : 'Mark Williams',
          'url' : 'https//github.com/markcwill',
          'packages' : ['hashpy', 'hashpy.db'],
          'package_data' : {'hashpy':['libhashpy.so','src/*.inc','Makefile','data/*', 'src/*.f']},
}

if 'ANTELOPE' in os.environ:
    ant_bin = os.path.join(os.environ['ANTELOPE'], 'bin')
    kwargs['data_files'] = [(ant_bin, ['hashpy/db/dbhash'])]

setup(**kwargs)
