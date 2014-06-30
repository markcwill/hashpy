#
# setup.py file for compiling HASH and installing hashpy
#
import sys
import os
import inspect

from numpy.distutils.core import setup, Extension


SETUP_DIRECTORY = os.path.dirname(os.path.abspath(inspect.getfile(
    inspect.currentframe())))

# Import the version string.
UTIL_PATH = os.path.join(SETUP_DIRECTORY, "hashpy", "util")
sys.path.insert(0, UTIL_PATH)
from version import get_git_version  # @UnresolvedImport
sys.path.pop(0)

#--- libhashpy Fortran extension --------------------------------------------#
#
# Build extension from FORTRAN source of HASH subroutines
# (based on the fucntion numpy.f2py.f2py2e.run_compile)
srcdir = os.path.join('hashpy', 'src')
srcf = ['fmamp_subs.f95', 'fmech_subs.f95', 'uncert_subs.f95', 'util_subs.f95',
        'pol_subs.f95', 'vel_subs.f95', 'station_subs.f95', 'vel_subs2.f95']
src_list = [os.path.join(srcdir, src) for src in srcf]
ext_args = {'sources': src_list}

#
# Use as a template for non-standard (non-distro) python installls...
#
def get_linker_args_for_virtualenv(virtualenv=None):
    """Return linker args relative to a virtual env"""
    np_inc = os.path.join('lib', 'python2.7', 'site-packages', 'numpy', 'core',
                          'include')
    inc_dirs = [os.path.join(virtualenv, inc) for inc in ('include', np_inc)]
    lib_dirs = [os.path.join(virtualenv, lib) for lib in ('lib',)]
    return  {'include_dirs': inc_dirs,
             'library_dirs': lib_dirs,
             }

# Have to link against antelope libs if installing to Antelope python
#
# TODO: get python version/path automagically.
if 'antelope' in sys.executable:
    python_folder = '/opt/antelope/python2.7.2-64'
    ANT_EXT_ARGS = get_linker_args_for_virtualenv(python_folder)
    ext_args.update(ANT_EXT_ARGS)
#----------------------------------------------------------------------------#


### SETUP ####################################################################
s_args = {'name'         : 'HASHpy',
          'version'      : get_git_version(),
          'description'  : 'Routines for running HASH algorithms',
          'author'       : 'Mark Williams',
          'url'          : 'https//github.com/markcwill/hashpy',
          'packages'     : ['hashpy', 'hashpy.io', 'hashpy.plotting'],
          'package_data' : {'hashpy': [
                                'RELEASE-VERSION',
                                'src/*.inc','src/Makefile','data/*',
                                'scripts/*', 'src/*.f95'
                                ]
                            },
          'ext_modules'  : [Extension('hashpy.libhashpy', **ext_args)],
}
##############################################################################


# hashpy.db -----------------------------------------------------------------#
# TODO: OBSELETED - break out to separate module/package so hashpy can
# stand alone.
#
# copy pf and bins from hashpy.db to antelope if available
if 'ANTELOPE' in os.environ:
    ant_bin = os.path.join(os.environ['ANTELOPE'], 'bin')
    ant_pf = os.path.join(os.environ['ANTELOPE'], 'data', 'pf')
    s_args['data_files'] = [(ant_bin, ['hashpy/scripts/dbhash']),
                            (ant_pf,  ['hashpy/data/dbhash.pf'])]
#----------------------------------------------------------------------------#


setup(**s_args)
