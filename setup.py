#
# setup.py file for compiling HASH and installing hashpy
#
from numpy.distutils.core import setup, Extension
import sys, os

#--- libhashpy SETUP ---------------------------------------------------------
#
# Build extension from FORTRAN source
# (based on the fucntion numpy.f2py.f2py2e.run_compile)
srcdir = 'hashpy/src'

srcf = ['fmech_subs.f', 'uncert_subs.f', 'util_subs.f',
        'pol_subs.f', 'vel_subs.f', 'station_subs.f', 'vel_subs2.f' ]

src_list = [ os.path.join(srcdir,src) for src in srcf ]

ext_args = { 'sources' : src_list }

#
# Use as a template for non-standard (non-distro) python installls...
#
def get_linker_args_for_virtualenv(virtualenv=None):

     return  {'include_dirs' : [ virtualenv + '/include',
                                 virtualenv + '/lib/python2.7/site-packages/numpy/core/include'],

              'library_dirs' : [ virtualenv + '/lib'], 
              }

# Have to link against antelope libs if installing to Antelope python
if 'antelope' in sys.executable:
    python_folder = '/opt/antelope/python2.7.2-64'
    ANT_EXT_ARGS  = get_linker_args_for_virtualenv(python_folder)
    ext_args.update(ANT_EXT_ARGS)

ext = Extension('hashpy.libhashpy', **ext_args)
#-----------------------------------------------------------------------------


### Regular setup stuff ######################################################

s_args = {'name'         : 'HASHpy',
          'version'      : '0.5.1',
          'description'  : 'Routines for running HASH algorithms',
          'author'       : 'Mark Williams',
          'url'          : 'https//github.com/markcwill',
          'packages'     : ['hashpy', 'hashpy.io', 'hashpy.plotting'],
          'package_data' : {'hashpy': ['src/*.inc','src/Makefile','data/*','scripts/*', 'src/*.f']},
          'ext_modules'  : [ext],
}

# hashpy.db ------------------------------------------------------------------
# TO BE OBSELETED - break out to separate module/package so hashpy can
# stand alone.
#
# copy pf and bins from hashpy.db to antelope if available
if 'ANTELOPE' in os.environ:
    ant_bin = os.path.join(os.environ['ANTELOPE'], 'bin')
    ant_pf  = os.path.join(os.environ['ANTELOPE'], 'data', 'pf')
    s_args['data_files'] = [(ant_bin, ['hashpy/scripts/dbhash']),
                            (ant_pf,  ['hashpy/data/dbhash.pf'])]
#-----------------------------------------------------------------------------

# Go
setup(**s_args)
