# HASHpype Makefile
#
# - Mark Williams 2012.320
# Based on HASH Makefile

SRC=    fmamp_subs.f95 fmech_subs.f95 uncert_subs.f95 util_subs.f95 \
        pol_subs.f95 vel_subs.f95 station_subs.f95 vel_subs2.f95

# HASHpype library module:
hashpy: libhashpy.so

libhashpy.so: Makefile $(SRC)
	f2py -c -m libhashpy $(SRC)     

#.PHONY : install
#install : libhashpy.so
#	@mv libhashpy.so ..

