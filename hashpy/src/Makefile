# HASHpy Makefile for Fortran subroutines
#
# - Mark Williams 2012.320
#
# Based on Makefile from HASH v1.2
#
# This is NOT used by the setup.py file, the subs are compiled using
# numpy.distutils as a Python extension. However, this file makes it
# possible to change the default values in the include files and compile
# a new module.
#
# This may be necessary until the code is migrated to using allocatable
# arrays, and/or passing program-wide lengths via common blocks or in the
# subroutines.
#
# This requires the 'f2py' script in your path.

SRC= fmamp_subs.f95 fmech_subs.f95 uncert_subs.f95 util_subs.f95 \
     pol_subs.f95 vel_subs.f95 station_subs.f95 vel_subs2.f95

# HASHpype library module:
libhashpy.so: Makefile $(SRC)
	f2py -c -m libhashpy $(SRC)     

.PHONY : install clean
install : libhashpy.so
	@mv libhashpy.so ..

clean :
	@rm -f libhashpy.so
