#
# setup.py file for compiling HASH and installing hashpy
#
import sys
import os
import inspect

from numpy.distutils.core import setup, Extension

#--- libhashpy Fortran extension --------------------------------------------#
#
# Build extension from FORTRAN source of HASH subroutines
# (based on the fucntion numpy.f2py.f2py2e.run_compile)
#
def hash_extension_args():
    srcdir = os.path.join('hashpy', 'src')
    srcf = ['fmamp_subs.f95', 'fmech_subs.f95', 'uncert_subs.f95', 'util_subs.f95',
            'pol_subs.f95', 'vel_subs.f95', 'station_subs.f95', 'vel_subs2.f95']
    src_list = [os.path.join(srcdir, src) for src in srcf]
    ext_args = {'sources': src_list}
    return ext_args

libhashpy_extension = Extension('hashpy.libhashpy', **hash_extension_args())


### SETUP ####################################################################
#
setup_args = {
    'name': 'hashpy',
    'version': '1.0.0-beta',
    'description': 'Python wrapper for HASH first-motion focal mech lib',
    'author': 'Mark Williams',
    'url': 'https//github.com/markcwill/hashpy',
    'packages': [
        'hashpy', 
        'hashpy.io', 
        'hashpy.io.obspy', 
    ],
    'package_data': {
        'hashpy': [
            'src/*.f95', 'src/*.inc', 'src/Makefile', 'data/*',
        ],
    },
    'ext_modules': [libhashpy_extension],
}

setup(**setup_args)
