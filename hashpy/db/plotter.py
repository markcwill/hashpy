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
from matplotlib.gridspec import GridSpec
import mplstereonet
from obspy.core import read, UTCDateTime as utc
from obspy.imaging.beachball import AuxPlane
from aug.contrib import open_db_or_string
from antelope.datascope import *
from obspy_ext.antelope import readANTELOPE

def get_waveform_from_arid(database, arid, window=4.):
	'''Return an ObsPy stream containing traces which match a given
	arrival ID from Datascope database
	
	Uses Mark's readANTELOPE function to load from wfdisc to Stream
	'''
	db,o = open_db_or_string(database)
	# Pull out needed params from 'arrival'
	dbv = dbprocess(db, ['dbopen arrival', 'dbsubset arid=={0}'.format(arid)])
	dbv.record = 0
	time, sta, chan = dbgetv(dbv,'time','sta','chan')
	# Use sta-chan-teim-endtime to pull out waveform
	t0 = utc(time)-(window/2.)
	t1 = t0 + window
	st = readANTELOPE(database, station=sta, channel=chan, 
		starttime=t0, endtime=t1)
	db.close()
	return st
	
	
class Plotter(object):
	'''Class to create interactive stereonet plot of a HASH first motion'''
	fig = None
	ax = None
	ax2 = None
	hro = None
	h_up = None
	h_down = None
	iup = None
	idn = None
	gs = None
	l = None
	pol = None
	
	picks_list = ('X', 'U', 'D')
	wf_color = { 'U' : 'red', 'D' : 'blue' }
	
			
	def onpick(self, event):
		'''Plot waveform of Z channel when a station is clicked
		
		When a first motion is clicked on the stereonet, get the
		waveform from the database and plot it. This will be the data
		associated with the first motion P pick.
		
		See the 'get_waveform_from_arid' fuction for more info 
		'''
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
		fig = self.fig 
		ax = self.ax2
		ax.clear()
		dataind = event.ind[0]
		k = inds[dataind]
		ax.set_title("{0} -- {1}".format(self.hro.sname[k],self.hro.arid[k]))
		ax.set_xlabel("{0}".format(pick))
		st = get_waveform_from_arid(self.hro.dbin, self.hro.arid[k], window=0.5)
		l, = ax.plot(st[0].data, color=self.wf_color[pick], lw=2)
		self.l = l
		plt.show()
		return True
	
	def enter_axes(self, event):
		'''When mouse enters waveform window'''
		if event.inaxes is self.ax2:
			event.inaxes.patch.set_facecolor('gray')
			event.canvas.draw()
	
	def leave_axes(self, event):
		'''When mouse leaves waveform window'''
		if event.inaxes is self.ax2:
			event.inaxes.patch.set_facecolor('white')
			event.canvas.draw()
		
	def switch_polarity(self, event):
		'''Change the first motion designation of a pick
		
		When waveform axes is clicked on, change the polarity to its
		opposite, e.g., from 'UP' to 'DOWN'. The waveform will change
		color correspondingly.
		'''
		if event.inaxes is self.ax2:
			pick = self.pol
			if pick == 'U':
				npick = 'D'
			else:
				npick = 'U'
			self.pol = npick
			self.l.set_color(self.wf_color[npick])
			self.ax2.set_xlabel("{0}".format(npick))
			event.canvas.draw()
	
	def plot_on_stereonet(self):
		hro = self.hro
		ax = self.ax
		# pull out variables from mechanism
		azimuths = hro.p_azi_mc[:hro.npol,0]
		#--- HASH takeoffs are 0-180 from vertical UP!!
		#--- Stereonet angles 0-90 inward (dip)
		#--- Classic FM's are toa from center???
		takeoffs = abs(hro.p_the_mc[:hro.npol,0] - 90)
		polarities = hro.p_pol[:hro.npol]
		
		# Planes
		strike1,dip1,rake1 = hro.str_avg[0], hro.dip_avg[0], hro.rak_avg[0]
		strike2,dip2,rake2 = AuxPlane(strike1, dip1, rake1)
		
		# Indices
		up = polarities > 0
		dn = polarities < 0
		
		# Plotting --------------------------------------#
		# plot trial planes (nout2) OR avg planes (nmult)
		if False:
			for n in range(hro.nout2):
				s1, d1, r1 = hro.strike2[n], hro.dip2[n], hro.rake2[n]
				s2, d2, r2 = AuxPlane(s1, d1, r1)
				h_rk = ax.plane(s1,d1, color='#999999')
				h_rk = ax.plane(s2,d2,'#888888')
		# plot best fit average plane, aux plane
		h_rk = ax.plane(strike1, dip1, color='black', linewidth=3)
		h_rk = ax.rake( strike1, dip1, -rake1, 'k^', markersize=8)
		h_rk = ax.plane(strike2, dip2, color='black', linewidth=3)
		# plot station takeoffs
		h_rk_up, = ax.rake(azimuths[up]-90.,takeoffs[up],90, 'o', picker=5, markersize=8, markeredgewidth=2, markeredgecolor='black', markerfacecolor='red')
		h_rk_dn, = ax.rake(azimuths[dn]-90.,takeoffs[dn],90, 'o', picker=5, markersize=8, markeredgewidth=2, markeredgecolor='blue',  markerfacecolor='white')
		# save to instance so other functions can change the plot
		self.h_up = h_rk_up
		self.h_down = h_rk_dn
		self.iup = arange(hro.npol)[up]
		self.idn = arange(hro.npol)[dn]
		# hack to throw in station names for temp debugging...
		if True:
			for i in range(hro.npol):
				h_rk = ax.rake(azimuths[i]-90,takeoffs[i]+5,90, marker='$   {0}$'.format(hro.sname[i]), color='black',markersize=20)
		#------------------------------------------------#
	
	def __init__(self,hro):
		'''Create a plot for focal mechanisms'''
		self.hro = hro
		
		# Draw figure and set up
		fig = plt.figure()
		gs = GridSpec(4,1)
		ax = fig.add_subplot(gs[:-1,:], projection='stereonet')
		ax.clear() 
		ax.set_title('{0} - click to plot station time series'.format(hro.icusp))
		tlab  = ax.set_azimuth_ticklabels([])
		
		# Save to plotting object
		self.fig = fig
		self.ax = ax
		self.gs = gs
		self.ax2 = fig.add_subplot(self.gs[-1,0])
		# Plot FM stuff
		self.plot_on_stereonet()
		
		# Set mouse events to method functions.
		fig.canvas.mpl_connect('pick_event', self.onpick)
		fig.canvas.mpl_connect('axes_enter_event', self.enter_axes)
		fig.canvas.mpl_connect('axes_leave_event', self.leave_axes)
		fig.canvas.mpl_connect('button_press_event', self.switch_polarity)
		# Save and draw
		plt.show()
	
