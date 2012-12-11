#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  focalmech.py
#
# by Mark Williams 2012.313
# First motion focal mechanism classes

import numpy as np
from obspy.core.util import AttribDict
from obspy.imaging.beachball import AuxPlane

# following two subs are from GMT utlimeca.c
#
# ported to python by Mark
# from perl port by Gabe
#
def zero_360(str1):
	'''Put an angle between 0 and 360 degrees
	Genevieve Patau
	'''
	if str1 >= 360:
		str1 -= 360
	elif str1 < 0:
		str1 += 360
	else:
		pass
	return str1

def nodal2pt(str1,da1,sa1,str2,da2,sa2):
	'''Compute azimuth and plungement of P-T axis 
	(from nodal plane strikes, dips and rakes.)
	
	Mark's python port from Gabe's perl port from:
	FORTRAN routines of Anne Deschamps ::
	
	Inputs
	------
	*args == (str1,da1,sa1,str2,da2,sa2)
	For each plane:
	str : strike angle in degrees
	da  : dip angle in degrees
	sa  : rake (slip angle) in degrees
	
	Returns
	-------
	Dips and azimuths of P and T axes
	(dip_p, dip_t, azi_p, azi_t)
	
	(Original fxn used azimuth of dip plane, not strike)
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


class NodalPlane(list):
	'''List to hold strike, dip, rake of a nodal plane
	
	Overview
	--------
	Basically, a list consisting of:
	[strike, dip, rake]
	with each element accessible by name as well as index.
	
	Construct with sequence, list, or named keyword args, see
	constructor doc for deatails.
	
	Attributes
	----------
	strike	:	int or float of degrees
	dip		:	int or float of degrees
	rake	:	int or float of degrees
	'''
	_imap = {'strike': 0,
			 'dip'   : 1,
			 'rake'  : 2}
	
	def __init__(self, *args, **kwargs):
		'''
		NodalPlane(strk, dp, rk)
		NodalPlane([strk,dp,rk])
		NodalPlane(strike=strk, dip=dp, rake=rk)
		
		strike	:	int or float of degrees
		dip		:	int or float of degrees
		rake	:	int or float of degrees
		
		Examples
		--------
		>>> l = [123, 45, 67]
		>>> p = NodalPlane(l)
		>>> p = NodalPlane(145, 45, 67)
		>>> p = NodalPlane(strike=145, dip=45, rake=67)
		>>> p.dip = 30
		'''
		super(NodalPlane,self).__init__([None,None,None])
		
		if args:
			if isinstance(args[0], list) or isinstance(args[0], tuple) and len(args[0]) == 3 :
				self[:] = args[0]
			elif len(args) == 3:
				self[:] = args
			else:
				pass
		if kwargs:
			for key in kwargs:
				if key in self._imap:
					self.__setattr__(key, kwargs[key])
	
	def __getindex(self, key):
		'''Maps an attribute key to a list index'''
		if key in self._imap:
			index = self._imap[key]
		else:
			index = None
		return index
	
	def __getattr__(self, key):
		'''Look for attribute in list'''
		index = self.__getindex(key)
		if index is not None:
			return self[index]
		else:
			raise AttributeError("Attribute must be 'strike', 'dip', or 'rake'")
	
	def __setattr__(self, key, value):
		'''Set attribute in list'''
		index = self.__getindex(key)
		if index is not None:
			self[index] = value
		else:
			raise AttributeError("Attribute must be 'strike', 'dip', or 'rake'")


class DoubleCouple(object):
	'''Holds nodal planes and P and T axes of a double couple focal mech
	
	The attributes are set up to calulate everything on the fly from the
	initial plane (strike, dip, rake), so you can change something (like
	a rake in your primary plane), and calling for a 'P' axis, e.g. will
	give you the new answer...
	
	Attributes
	----------
	plane1	:	NodalPlane of primary plane
	plane2	:	NodalPlane of auxiliary plane
	axis	:	AttribDict of axes ('P' and 'T')
					containing list of [azi, dip]
	'''
	
	plane1 = None
	
	@property
	def plane2(self):
		'''Return Auxiliary plane'''
		auxplane = AuxPlane(*self.plane1)
		return NodalPlane(*auxplane)
	
	@property
	def axis(self):
		'''return direction and dip for P and T axes'''
		dipP, dipT, aziP, aziT = nodal2pt(*self.plane1+self.plane2)
		return AttribDict({'P': [aziP, dipP], 'T': [aziT, dipT] })
	
	def __init__(self, nodal_plane=None):
		self.plane1 = nodal_plane

#
# Works in progress...
#

class FirstMotionFM(object):
	'''Stub'''
	planes = None # AttribDict of 2 obspy NodalPlanes
	picks  = None # array of station/azimuth/takeoff/polarities
	data   = None # a HashRun object, e.g.
	_dt = np.dtype([('station', 'a6'), ('azimuth', float), ('takeoff',float), ('polarity','a2')])
	
	def __init__(self):
		'''Build'''
		planes = AttribDict({'primary': [], 'aux': []})
		
	def load_hash(self, hro=HashRun() ):
		'''Map HASH variables to data for methods to use
		
		If you don't specify a HashRun object these fields will just be junk
		'''
		
		assert isinstance(hro, HashRun), "Must pass a HashRun in here to load data!"
		
		n = hro.npol
		picks = np.empty(n, dtype=_dt)
		picks['station'] = hro.sname[:n]
		picks['azimuth'] = hro.p_azi_mc[:n,0]
		picks['takeoff'] = hro.p_the_mc[:n,0]
		picks['polarity'] = hro.p_pol[:n]
		dc = DoubleCouple(NodalPlane(hro.str_avg[0], hro.dip_avg[0], hro.rak_avg[0]))
		plane1 = dc.plane1
		plane2 = dc.plane2
		
		self.data = hro
		self.picks = picks
		self.planes = {'main':plane1, 'aux':plane2}
		
	def load_fplane(self, orid=[]):
		'''stub'''
	
	def plot(self, labels = False):
		'''stub'''
		import matplotlib.pyplot as plt
		import mplstereonet
		
		fig = plt.figure()
		ax = fig.add_subplot(111, projection='stereonet')
		# pull out variables from mechanism
		azimuths = self.picks['azimuth']
		# HASH takeoffs are 0-180 from vertical UP!!
		# Stereonet angles 0-90 inward (dip)
		# Classic FM's are toa from center???
		takeoffs = abs(self.picks['takeoff'] - 90)
		polarities = self.picks['polarity']
		strike1,dip1,rake1 = self.planes['main']
		strike2,dip2,rake2 = self.planes['aux']
		up = polarities > 0
		dn = polarities < 0
		if False:
			# plot HASH trial planes (nout2) OR avg planes (nmult)
			for n in range(self.data.nout2):
				s1, d1, r1 = self.data.strike2[n], self.data.dip2[n], self.data.rake2[n]
				s2, d2, r2 = AuxPlane(s1, d1, r1)
				h_rk = ax.plane(s1,d1, color='#999999')
				h_rk = ax.plane(s2,d2,'#888888')
		# plot best fit average plane
		h_rk = ax.plane(strike1, dip1, color='black', linewidth=3)
		h_rk = ax.rake( strike1, dip1, -rake1, 'k^', markersize=8)
		h_rk = ax.plane(strike2, dip2, color='black', linewidth=3)
		# plot station takeoffs
		h_rk = ax.rake(azimuths[up]-90.,takeoffs[up],90, 'wo', markersize=8, markeredgewidth=2, markerfacecolor=None)
		h_rk = ax.rake(azimuths[dn]-90.,takeoffs[dn],90, 'ko', markersize=8, markeredgewidth=2)
		# hack to throw in station names for temp debugging...
		if labels:
			for i in range(len(self.picks)):
				h_rk = ax.rake(azimuths[i]-90,takeoffs[i]+5,90, marker='$   {0}$'.format(self.picks[i]['station']), color='black',markersize=20)
		# and go.
		plt.show()
