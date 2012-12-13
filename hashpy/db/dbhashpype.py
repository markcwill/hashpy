#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# dbhashpype.py
#
# by Mark Williams 2012.313
#
# Class to run HASHpy using Antelope database I/O
#
import os.path
import numpy as np
from hashpy import HashPype
from antelope.stock import pfget
from antelope.datascope import dblookup, dbprocess
from aug.contrib import AttribDbptr, open_db_or_string, DbrecordList
from obspy.imaging.beachball import AuxPlane

# for HASH compatiblity, change later.
degrad = 180. / np.pi
rad = 1. / degrad


class DbHashPype(HashPype):
	'''Object which holds all data from a HASH instance for one event
	
	The variables are named and structured as close to the original code
	as possible. Methods act on the variables within this namespace and
	most are accessible as attributes.
	
	One can make a "HASH" driver program by calling the methods of this
	class on the data held within it.
	'''
	
	dbin = None
	dbout = None
	
	def load_pf(self, pffile='dbhash.pf'):
		'''Update runtime settings from a pf file
		
		One can also specify names of velocity  model files in the pf.
		
		Right now these settings are inherited from the HashPype class,
		 and are not instance attributes.
		 
		 Input
		 -----
		 pffile : string of full path to pf file
		 
		'''
		pf_settings = pfget(pffile)
		
		# Little hack to do type conversions 
		for key in pf_settings:
			pfi = pf_settings[key]
			if key in ['badfrac','prob_max']:
				pfi = float(pfi)
			elif key in ['npolmin','max_agap','max_pgap','dang','nmc','maxout', 'delmax','cangle']:
				pfi = int(pfi)
			else:
				pass
			self.__setattr__(key, pfi)
		
		if 'vmodel_dir' in pf_settings and 'vmodels' in pf_settings:
			self.vmodels = [os.path.join(self.vmodel_dir, table) for table in self.vmodels]
	
	def get_phases_from_db(self, dbname, evid=None, orid=None):
		'''Input HASH data from Antelope database
		
		This will accept a database name OR Antelope Dbptr, and either
		an ORID, or an EVID, in which case the 'prefor' ORID is looked
		up and used.
		
		Inputs
		------
		dbname	:	str or antelope.datascope.Dbptr
		orid	:	int of ORID
		evid	:	int of EVID
		'''
		self.dbin = dbname
		
		db, oflag = open_db_or_string(dbname)
		if orid is None:
			dbv = dbprocess(db,['dbopen event', 'dbsubset evid == '+str(evid)])
			dbv.record = 0
			orid = dbv.getv('prefor')[0]
		db = dbprocess(db,[ 'dbopen origin', 'dbsubset orid == '+str(orid),
						'dbjoin origerr', 'dbjoin assoc',  'dbjoin arrival',
						'dbjoin affiliation', 'dbjoin site',
						'dbsubset iphase =~ /.*[Pp].*/',
						'dbsubset (ondate <= time)',
						'dbsubset (time <= offdate) || (offdate == -1)']
						)
		
		phases = AttribDbptr(db)
		self.nrecs = len(phases)
		assert len(phases) > 0, "No picks for this ORID: {0}".format(orid)
		ph = phases[0]
		self.tstamp = ph['origin.time']
		self.qlat   = ph['origin.lat']
		self.qlon   = ph['origin.lon']
		self.qdep   = ph['origin.depth']
		self.qmag   = ph['origin.ml']
		self.icusp  = ph['origin.orid']
		self.seh	= ph['origerr.smajax']
		self.sez	= ph['origerr.sdepth']
		
		aspect = np.cos(self.qlat / degrad) # convert using python later.
		
		# The index 'k' is deliberately non-Pythonic to deal with the fortran
		# subroutines which need to be called and the structure of the original HASH code.
		# May be able to update with a rewrite... YMMV
		k = 0
		for ph in phases:
			# load up params
			# in future, could use the acol() method?
			self.sname[k]	 = ph['sta']
			self.snet[k]	  = ph['net']
			self.scomp[k]	 = ph['chan']
			self.pickonset[k] = 'I'
			self.pickpol[k]   = ph['fm']
			self.arid[k]	  = ph['arid']
			
			flat, flon, felv = ph['site.lat'],ph['site.lon'],ph['site.elev']
			self.esaz[k] = ph['esaz']

			#print '{0} {1} {2} {3} {4} {5} {6} {7}'.format(k, sname[k], snet[k], scomp[k], pickonset[k], pickpol[k], flat, flon)
			
			# dist @ azi, get from db OR obspy or another python mod (antelope) could do this on WGS84
			dx = (flon - self.qlon) * 111.2 * aspect
			dy = (flat - self.qlat) * 111.2
			dist = np.sqrt(dx**2 + dy**2)
			qazi = 90. - np.arctan2(dy,dx) * degrad
			
			if (qazi < 0.):
				qazi = qazi + 360.
			if (dist > self.delmax):
				continue
			if (self.pickpol[k] in 'CcUu'):
				self.p_pol[k] = 1
			elif (self.pickpol[k] in 'RrDd'):
				self.p_pol[k] = -1
			else:
				continue
			
			# save them for other functions -MCW
			self.dist[k] = dist
			self.qazi[k] = qazi
			self.flat[k] = flat
			self.flon[k] = flon
			self.felv[k] = felv
			
			if (self.pickonset[k] in 'Ii'):
				self.p_qual[k] = 0
			else:
				self.p_qual[k] = 1
			
			# polarity check in original code... doesn't work here
			#self.p_pol[k] = self.p_pol[k] * self.spol
			k += 1
		#npol = k - 1
		self.npol = k # k is zero indexed in THIS loop
		db.close()
	
	def save_result_to_db(self, dbout=None, solution=0):
		'''Write the preferred HASH solution to Datascope database.
		
		This writes the strike, dip, rakes to 'fplane', arids used for a
		given mech in 'predmech' and the takeoffs in 'predarr'.
		
		Currently depends on a HashPype hack function which keeps this
		internally in a list of dictionaries... (for fplane, 
		others are from the object attr, I'm working on it)
		
		Input
		-----
		dbout	: str or antelope.datascope.Dbptr to database
		solution : <STUB> int of desired solution.
		
		'''
		
		fp = self.fplane[solution]
		
		db, oflag = open_db_or_string(dbout, perm='r+')
		mechid = db.nextid('mechid')
		
		dbfpln = dblookup(db,table='fplane')
		
		dbfpln.record = dbfpln.addnull()
		dbfpln.putv('orid', fp['orid'],
			   'str1', round(fp['str1'],1) ,
			   'dip1', round(fp['dip1'],1) ,
			   'rake1',round(fp['rake1'],1),
			   'algorithm', fp['algorithm'],
			   'mechid', mechid
			   )
		if True:
			fp['str2'],fp['dip2'],fp['rake2'] = AuxPlane(fp['str1'],fp['dip1'],fp['rake1'])
			dbfpln.putv('str2', round(fp['str2'],1) ,
			   'dip2', round(fp['dip2'],1) ,
			   'rake2',round(fp['rake2'],1))
		
		dbpmec = dblookup(db,table='predmech')
		dbparr = dblookup(db,table='predarr')
		for k in range(self.npol):
			if self.p_pol[k] > 0:
				fm = 'U'
			else:
				fm = 'D'
			dbpmec.record = dbpmec.addnull()
			dbpmec.putv('arid', int(self.arid[k]) ,
						'orid', self.icusp,
						'mechid', mechid,
						'fm', fm,
						)
			dbparr.record = dbparr.addnull()
			dbparr.putv('arid', int(self.arid[k]),
						'orid', self.icusp, 
						'esaz', self.qazi[k], 
						'dip' , self.p_the_mc[k,0],
						)
		db.close()
	
	
	# This one can probably be booted, really should be abstracted out
	def read_result_from_db(self, dbin=None, orid=None):
		'''Read in a mechanism from the fplane table'''
		if orid is None:
			orid = self.icusp
		d = open_db_or_string(dbin)
		d = d.lookup(table='fplane')
		d = d.subset('orid == {0}'.format(self.icusp))
		assert d.nrecs() is not 0, 'No solution for this ORID: {0}'.format(orid)
		self.fplane = DbrecordList(d)
		


def plotter(hro):
	from matplotlib import pyplot as plt
	import mplstereonet
	from obspy.imaging.beachball import AuxPlane
	
	def onpick(event):
	
		if event.artist!=line: return True
	
		N = len(event.ind)
		if not N: return True
	
	
		figi = plt.figure()
		for subplotnum, dataind in enumerate(event.ind):
			ax = figi.add_subplot(N,1,subplotnum+1)
			ax.plot(hro.p_pol[dataind])
			ax.text(0.05, 0.9, 'Station {0}'.format(hro.sname[dataind]),
					transform=ax.transAxes, va='top')
			ax.set_ylim(-0.5, 1.5)
		figi.show()
		return True
	
	fig = plt.figure()
	ax = fig.add_subplot(111, projection='stereonet')
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
	h_rk_up = ax.rake(azimuths[up]-90.,takeoffs[up],90, 'ko', picker=5, markersize=8, markeredgewidth=2, markerfacecolor=None)
	h_rk_dn = ax.rake(azimuths[dn]-90.,takeoffs[dn],90, 'wo', picker=5, markersize=8, markeredgewidth=2)
	#h_t  = ax.set_title("ORID: {0}".format(hro.icusp))
	# hack to throw in station names for temp debugging...
	if False:
		for i in range(hro.npol):
			h_rk = ax.rake(azimuths[i]-90,takeoffs[i]+5,90, marker='$   {0}$'.format(hro.sname[i]), color='black',markersize=20)
	# and go.
	fig.canvas.mpl_connect('pick_event', onpick)
	plt.show()
