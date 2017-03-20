#
# setup.py file for compiling HASH and installing hashpy
#
import sys
import os
import inspect
import setuptools #enables develop

from numpy.distutils.core import setup, Extension


SETUP_DIRECTORY = os.path.dirname(os.path.abspath(inspect.getfile(
    inspect.currentframe())))

# Import the version string.
UTIL_PATH = os.path.join(SETUP_DIRECTORY, "hashpy", "util")
sys.path.insert(0, UTIL_PATH)
from version import get_git_version  # @UnresolvedImport
sys.path.pop(0)

#
# For non-system python installs, or non-activated virtualenvs.
#
def get_linker_args_for_virtualenv(virtualenv=None):
    """Return linker args relative to a virtual env"""
    np_inc = os.path.join('lib', 'python2.7', 'site-packages', 'numpy', 'core',
                          'include')
    inc_dirs = [os.path.join(virtualenv, inc) for inc in ('include', np_inc)]
    lib_dirs = [os.path.join(virtualenv, lib) for lib in ('lib',)]
    return  {
        'include_dirs': inc_dirs,
        'library_dirs': lib_dirs,
        }

#--- libhashpy Fortran extension --------------------------------------------#
#
# Build extension from FORTRAN source of HASH subroutines
# (based on the fucntion numpy.f2py.f2py2e.run_compile)
#
def get_extension_args():
    srcdir = os.path.join('hashpy', 'src')
    srcf = ['fmamp_subs.f95', 'fmech_subs.f95', 'uncert_subs.f95', 'util_subs.f95',
            'pol_subs.f95', 'vel_subs.f95', 'station_subs.f95', 'vel_subs2.f95']
    src_list = [os.path.join(srcdir, src) for src in srcf]
    ext_args = {'sources': src_list}

    # Have to link against antelope libs if installing to Antelope python
    #
    if 'antelope' in sys.executable:
        python_folder = sys.executable.rstrip("/bin/python")
        ANT_EXT_ARGS = get_linker_args_for_virtualenv(python_folder)
        ext_args.update(ANT_EXT_ARGS)

    return ext_args
#----------------------------------------------------------------------------#


### SETUP ####################################################################
ext_args = get_extension_args()

setup_args = {
    'name'         : 'HASHpy',
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

# dbhash --------------------------------------------------------------------#
# Copy pf and bins to antelope if available
if 'ANTELOPE' in os.environ:
    ant_bin = os.path.join(os.environ['ANTELOPE'], 'bin')
    ant_pf = os.path.join(os.environ['ANTELOPE'], 'data', 'pf')
    setup_args['data_files'] = [(ant_bin, ['hashpy/scripts/dbhash']),
                                (ant_pf,  ['hashpy/data/dbhash.pf'])]
#----------------------------------------------------------------------------#

##############################################################################

setup(**setup_args)
