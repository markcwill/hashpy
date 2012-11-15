#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# utils.py


import numpy as np

# following two subs are from GMT utlimeca.c
def zero_360(str1):
	'''Stub'''
	#put an angle between 0 and 360 degrees
	#Genevieve Patau
	if str1 >= 360:
		str1 -= 360
	elif str1 < 0:
		str1 += 360
	else:
		pass
	return str1

# FORTRAN routines of Anne Deschamps :
# compute azimuth and plungement of P-T axis
# from nodal plane strikes, dips and rakes.
def ps_pt_axis(str1,da1,sa1,str2,da2,sa2):
	'''
	my $dd1 = shift @_;	 # dip direction principle plane
	my $da1 = shift @_;	 # dip angle principle plane
	my $sa1 = shift @_;	 # rake principle plane
	my $dd2 = shift @_;	 # dip direction aux plane
	my $da2 = shift @_;	 # dip angle aux plane
	my $sa2 = shift @_;	 # rake aux plane
	str1 = dd1 - 90;
	str2 = dd2 - 90;
	'''
	# Constants, mostly unnecessary, fix later:
	# e.g. M_PI = np.pi
	EPSIL   = .0001
	M_PI    = 3.14159265358979323846
	M_PI_2  = 1.57079632679489661923
	M_PI_4  = 0.78539816339744830962
	M_SQRT2 = 1.41421356237309504880
	TWO_PI  = 6.28318530717958647692
	D2R     = M_PI / 180.0

	#my ($pp, $dp, $pt, $dt, $xp, $yp);

	im = 0
	pure_strike_slip = 0

	if abs(np.sin(sa1 * D2R)) > EPSIL:
		im = sa1 / abs(sa1)
	elif abs(np.sin(sa2 * D2R)) > EPSIL:
		im = sa2 / abs(sa2)
	else:
		pure_strike_slip = 1

	if pure_strike_slip:
		if np.cos(sa1 * D2R) < 0:
			pp = zero_360(str1 + 45)
			pt = zero_360(str1 - 45)
		else:
			pp = zero_360(str1 - 45);
			pt = zero_360(str1 + 45);
		dp = 0
		dt = 0
	else:
		cd1 = np.cos(da1 * D2R) *  M_SQRT2
		sd1 = np.sin(da1 * D2R) *  M_SQRT2
		cd2 = np.cos(da2 * D2R) *  M_SQRT2
		sd2 = np.sin(da2 * D2R) *  M_SQRT2
		cp1 = -(np.cos(str1 * D2R)) * sd1
		sp1 = np.sin(str1 * D2R) * sd1
		cp2 = -(np.cos(str2 * D2R)) * sd2
		sp2 = np.sin(str2 * D2R) * sd2

		amz = -(cd1 + cd2)
		amx = -(sp1 + sp2)
		amy =   cp1 + cp2
		dx  = np.arctan2(np.sqrt(amx * amx + amy * amy), amz) - M_PI_2
		px  = np.arctan2(amy, -amx)

		if px < 0:
			px += TWO_PI
		
		amz   = cd1 - cd2
		amx   = sp1 - sp2
		amy   = - cp1 + cp2
		dy = np.arctan2(np.sqrt(amx * amx + amy * amy), -abs(amz)) - M_PI_2
		py = np.arctan2(amy, -amx)

		if amz > 0:
			py -= M_PI
		if py < 0:
			py += TWO_PI

		if im == 1:
			dp = dy
			pp = py
			dt = dx
			pt = px
		else:
			dp = dx
			pp = px
			dt = dy
			pt = py

		pp *= 180 / M_PI
		dp *= 180 / M_PI
		pt *= 180 / M_PI
		dt *= 180 / M_PI
	
	# I added this line b/c the names are confusing - MCW
	dip_p, dip_t, azi_p, azi_t = dp, dt, pp, pt
	return dip_p, dip_t, azi_p, azi_t

