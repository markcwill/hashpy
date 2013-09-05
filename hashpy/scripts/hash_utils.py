# -*- coding: utf-8 -*-
"""
  hash_utils.py

 -by Mark Williams (2013), NSL
 
 utilities for running HASH using the python version of hash_driver2
 (hash_driver2.py)
"""

def parameter(**kwargs):
	'''returns variables inside a fortran 'parameter' call'''
	# FUTURE: could just make them all globals and import from namespace
	return [kwargs[key] for key in kwargs]

def fortran_include(fname):
	'''functions similar to a fortran include'''
	vars = []
	f_inc = open(fname)
	for line in f_inc:
		if 'c' in line[0]:
			pass
		elif 'parameter' in line:
			vars.extend(eval(line))
		else:
			pass
	f_inc.close()
	return vars
	
def get_sta_coords(stfile):
	'''Stub which reads in the HASH example station file'''
	f = open(stfile)
	sta_coords = {}
	for s in f:
		sta = s[0:4]
		lat = float(s[42:50])
		lon = float(s[51:61])
		elv = float(s[62:67])/1000.
		if sta in sta_coords:
			pass
		else:
			sta_coords[sta] = [lat,lon,elv]
	f.close()
	return sta_coords

def test_stereo(azimuths,takeoffs,polarities,sdr=[]):
	'''
	Plots points with given azimuths, takeoff angles, and 
	polarities on a stereonet. Will also plot both planes
	of a double-couple given a strike/dip/rake
	'''
	import matplotlib.pyplot as plt
	import mplstereonet
	from obspy.imaging.beachball import AuxPlane
	
	fig = plt.figure()
	ax = fig.add_subplot(111, projection='stereonet')
	up = polarities > 0
	dn = polarities < 0
	h_rk = ax.rake(azimuths[up]-90.,takeoffs[up],90, 'ro')
	h_rk = ax.rake(azimuths[dn]-90.,takeoffs[dn],90, 'b+')
	if sdr:
		s2,d2,r2 = AuxPlane(*sdr)
		h_rk = ax.plane(sdr[0],sdr[1],'g')
		h_rk = ax.rake(sdr[0],sdr[1],-sdr[2], 'go')
		h_rk = ax.plane(s2,d2, 'g')
	plt.show()
	


