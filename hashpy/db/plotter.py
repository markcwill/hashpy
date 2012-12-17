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
from matplotlib.widgets import Button
import mplstereonet
from obspy.core import read, UTCDateTime as utc
from obspy.imaging.beachball import AuxPlane
from aug.contrib import open_db_or_string
from antelope.datascope import *
from obspy_ext.antelope import readANTELOPE

def get_waveform_from_arid(database, arid, window=4.):
	db,o = open_db_or_string(database)
	dbv = dbprocess(db, ['dbopen arrival', 'dbsubset arid=={0}'.format(arid)])
	dbv.record = 0
	time, sta, chan = dbgetv(dbv,'time','sta','chan')
	t0 = utc(time)-(window/2.)
	t1 = t0 + window
	st = readANTELOPE(database, station=sta, channel=chan, 
		starttime=t0, endtime=t1)
	db.close()
	return st
	
	
class Plotter(object):
	
	fig = None
	ax = None
	hro = None
	h_up = None
	h_down = None
	iup = None
	idn = None
	rax = None
	l = None
	pol = None
	
	picks_list = ('X', 'U', 'D')
	wf_color = { 'U' : 'red', 'D' : 'blue' }
	
			
	def onpick(self, event):
		N = len(event.ind)
		if not N: return True
		
		if event.artist is self.h_up:
			inds = self.iup
			pick = 'U'
		elif event.artist is self.h_down:
			inds = self.idn
			pick = 'D'
		else:
			return True
		
		self.pol = pick
		fig = plt.figure()
		ax = fig.add_subplot(111) #212
		#plt.subplots_adjust(left=0.3)
		ax.clear()
		for subplotnum, dataind in enumerate(event.ind):
			k = inds[dataind]
			ax.set_title("{0} -- {1}".format(self.hro.sname[k],self.hro.arid[k]))
			ax.set_xlabel("{0}".format(pick))
			st = get_waveform_from_arid(self.hro.dbin, self.hro.arid[k], window=0.5)
			#ax.set_position([0.3, 0, 1, 1])
			l, = ax.plot(st[0].data, color=self.wf_color[pick], lw=2)
			self.ax = ax
			self.l = l
		
		#rax = fig.add_axes([0.05, 0.7, 0.15, 0.15], axisbg='white')
		#radio = Button(rax, 'Flip')
		#radio.on_clicked(self.fmpick)
		#self.rax = rax
		fig.canvas.mpl_connect('axes_enter_event', self.enter_axes)
		fig.canvas.mpl_connect('axes_leave_event', self.leave_axes)
		fig.canvas.mpl_connect('button_press_event', self.switch_polarity)
		plt.show()
		return True
	
	def enter_axes(self, event):
		#print 'enter_axes', event.inaxes
		event.inaxes.patch.set_facecolor('yellow')
		event.canvas.draw()
	
	def leave_axes(self, event):
		#print 'leave_axes', event.inaxes
		event.inaxes.patch.set_facecolor('white')
		event.canvas.draw()
		
	def switch_polarity(self, event):
		pick = self.pol
		if pick == 'U':
			npick = 'D'
		else:
			npick = 'U'
		self.pol = npick
		self.l.set_color(self.wf_color[npick])
		self.ax.set_xlabel("{0}".format(npick))
		event.canvas.draw()
	
	def __init__(self,hro):
		self.hro = hro
		fig = plt.figure()
		ax = fig.add_subplot(111, projection='stereonet')
		ax.set_title('{0} - click to plot station time series'.format(hro.icusp))
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
		h_rk_dn, = ax.rake(azimuths[dn]-90.,takeoffs[dn],90, 'wo', picker=5, markersize=8, markeredgewidth=2, markerfacecolor=None)
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
	
