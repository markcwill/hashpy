#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# dbhashpype.py
#
# by Mark Williams 2012.313
#
# Class to run HASHpy using Antelope database I/O
#
#import os.path
#import numpy as np
#from hashpy import HashPype
#from antelope.stock import pfget
#from antelope.datascope import dblookup, dbprocess
#from aug.contrib import AttribDbptr, open_db_or_string, DbrecordList
#from obspy.imaging.beachball import AuxPlane

from numpy import arange
from matplotlib import pyplot as plt
import mplstereonet
from obspy.core import read, UTCDateTime as utc
from obspy.imaging.beachball import AuxPlane
from aug.contrib import open_db_or_string
from antelope.datascope import *

def get_waveform_from_arid(database, arid, window=4):
	db,o = open_db_or_string(self.dbin)
	db = dbprocess(db, ['dbopen arrival', 'dbsubset arid=={0}'.format(arid)])
	db.record = 0
	time, sta, chan = dbgetv(db,'time','sta','chan')
	db = dbprocess(db,['dbopen wfdisc',
	'dbsubset time >= {0} && endtime <= {0}'.format(time),
	'dbsubset sta=={0} && chan=={1}'.format(sta,chan)])
	db.close()
	
class Plotter(object):
	
	fig = None
	hro = None
	h_up = None
	h_down = None
	iup = None
	idn = None
	
	def onpick(self, event):
		N = len(event.ind)
		if not N: return True
		
		if event.artist is self.h_up:
			inds = self.iup
		elif event.artist is self.h_down:
			inds = self.idn
		else:
			return True
		
		ax = self.fig.add_subplot(212)
		ax.clear()
		for subplotnum, dataind in enumerate(event.ind):
			k = inds[dataind]
			ax.set_title('{0}'.format(self.hro.sname[k]))
			ax.plot(0.5, abs(self.hro.p_pol[k]),'bo')
		self.fig.show()
		return True
	
	def __init__(self,hro):
		self.hro = hro
		fig = plt.figure()
		ax = fig.add_subplot(211, projection='stereonet')
		ax.set_title('Click on station to plot time series')
		tlab  = ax.set_azimuth_ticklabels([])
		
		# pull out variables from mechanism
		azimuths = hro.p_azi_mc[:hro.npol,0]
		# HASH takeoffs are 0-180 from vertical UP!!
		# Stereonet angles 0-90 inward (dip)
		# Classic FM's are toa from center???
		takeoffs = abs(hro.p_the_mc[:hro.npol,0] - 90)
		polarities = hro.p_pol[:hro.npol]
		strike1,dip1,rake1 = hro.str_avg[0], hro.dip_avg[0], hro.rak_avg[0]
		strike2,dip2,rake2 = AuxPlane(strike1, dip1, rake1)
		up = polarities > 0
		dn = polarities < 0
		if False:
			# plot trial planes (nout2) OR avg planes (nmult)
			for n in range(hro.nout2):
				s1, d1, r1 = hro.strike2[n], hro.dip2[n], hro.rake2[n]
				s2, d2, r2 = AuxPlane(s1, d1, r1)
				h_rk = ax.plane(s1,d1, color='#999999')
				h_rk = ax.plane(s2,d2,'#888888')
		# plot best fit average plane
		h_rk = ax.plane(strike1, dip1, color='black', linewidth=3)
		h_rk = ax.rake( strike1, dip1, -rake1, 'k^', markersize=8)
		h_rk = ax.plane(strike2, dip2, color='black', linewidth=3)
		# plot station takeoffs
		h_rk_up, = ax.rake(azimuths[up]-90.,takeoffs[up],90, 'ko', picker=5, markersize=8, markeredgewidth=2, markerfacecolor=None)
		h_rk_dn, = ax.rake(azimuths[dn]-90.,takeoffs[dn],90, 'wo', picker=5, markersize=8, markeredgewidth=2)
		self.h_up = h_rk_up
		self.h_down = h_rk_dn
		self.iup = arange(hro.npol)[up]
		self.idn = arange(hro.npol)[dn]
		#h_t  = ax.set_title("ORID: {0}".format(hro.icusp))
		# hack to throw in station names for temp debugging...
		if True:
			for i in range(hro.npol):
				h_rk = ax.rake(azimuths[i]-90,takeoffs[i]+5,90, marker='$   {0}$'.format(hro.sname[i]), color='black',markersize=20)
		# and go.
		fig.canvas.mpl_connect('pick_event', self.onpick)
		self.fig = fig
		plt.show()
	
