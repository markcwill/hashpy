#!/usr/bin/env python

from numpy.distutils.core import setup, Extension
import os

# Build extension based on the fucntion numpy.f2py.f2py2e.run_compile
srcf = ['fmech_subs.f', 'uncert_subs.f', 'util_subs.f',
        'pol_subs.f', 'vel_subs.f', 'station_subs.f', 'vel_subs2.f' ]

src_list = ['hashpy/src/' + src for src in srcf]

ext = Extension('hashpy.libhashpy', sources=src_list)

# Regular setup stuff
s_args = {'name'         : 'HASHpy',
          'version'      : '2.0',
          'description'  : 'Routines for running HASH algorithms',
          'author'       : 'Mark Williams',
          'url'          : 'https//github.com/markcwill',
          'packages'     : ['hashpy', 'hashpy.db'],
          'package_data' : {'hashpy':['src/*.inc','src/Makefile','data/*', 'src/*.f']},
          'ext_modules'  : [ext],
}

# copy pf and bins from hashpy.db to antelope if available
if 'ANTELOPE' in os.environ:
    ant_bin = os.path.join(os.environ['ANTELOPE'], 'bin')
    ant_pf  = os.path.join(os.environ['ANTELOPE'], 'data', 'pf')
    s_args['data_files'] = [(ant_bin, ['hashpy/db/dbhash']),
                            (ant_pf,  ['hashpy/db/dbhash.pf'])]
# Go
setup(**s_args)
