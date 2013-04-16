# HASHpype Makefile
#
# - Mark Williams 2012.320
# Based on HASH Makefile

SRC=    fmech_subs.f uncert_subs.f util_subs.f \
        pol_subs.f vel_subs.f station_subs.f vel_subs2.f

# HASHpype library module:
hashpy: libhashpy.so

libhashpy.so: Makefile $(SRC)
	f2py -c -m libhashpy $(SRC)     

#.PHONY : install
#install : libhashpy.so
#	@mv libhashpy.so ..

