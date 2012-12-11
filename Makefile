# HASHpype Makefile
#
# - Mark Williams 2012.320
# Based on HASH Makefile

SRCF=    fmech_subs.f uncert_subs.f util_subs.f \
        pol_subs.f vel_subs.f station_subs.f vel_subs2.f

SRC= $(addprefix hashpy/src/, $(SRCF))

# HASHpype library module:
hashpy: libhashpy.so

libhashpy.so: Makefile $(SRC)
	f2py -c -m libhashpy $(SRC)
	mv libhashpy.so hashpy
    
.PHONY : install
install : libhashpy.so
	python setup.py install

