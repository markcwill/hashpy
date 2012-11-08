#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  hashpype.py

from hashpy import *
import numpy as np
from hash_utils import fortran_include, get_sta_coords, test_stereo
from ant_utils import add_antelope_path, open_db_or_string
from obspy_ext.antelope import *
add_antelope_path()

# for HASH compatiblity, change later.
degrad = 180. / np.pi
rad = 1. / degrad

class FocalMech(object):
	'''First motion focal mechanisms'''
	# These MUST be the same as the fortran includes!!
	# (They are compiled into the Fortran subroutines)
	npick0, nmc0, nmax0 = None, None, None
	dang0, ncoor        = None, None
	
	# initialize arrays
	
	# Input arrays
	sname     = None
	scomp     = None
	snet      = None
	pickpol   = None
	pickonset = None
	p_pol     = None
	p_qual    = None
	spol      = None
	p_azi_mc  = None
	p_the_mc  = None
	index     = None
	qdep2     = None
	
	# Output arrays
	f1norm  = None
	f2norm  = None
	strike2 = None
	dip2    = None
	rake2   = None
	str_avg = None
	dip_avg = None
	rak_avg = None
	var_est = None
	var_avg = None
	mfrac   = None
	stdr    = None
	prob    = None
	qual    = None
	
	# Chosen fm run setings: (pf file should look like this:)
	# dbhash.pf -----------------------------------------------------
	# Defaults
	npolmin  = 8    # Enter mininum number of polarities (e.g., 8)
	max_agap = 90   # Enter maximum azimuthal gap (e.g., 90)
	max_pgap = 60   # Enter maximum takeoff angle gap (e.g., 60)
	dang     = 5    # Enter grid angle for focal mech search, in degrees 5(min {0})
	nmc      = 30   # Enter number of trials (e.g., 30)
	maxout   = 500  # Enter maxout for focal mech. output (e.g., 500)
	badfrac  = 0.1  # Enter fraction of picks presumed bad (e.g., 0.10)
	delmax   = 500  # Enter maximum allowed source-station distance, in km (e.g., 120)
	cangle   = 45   # Enter angle for computing mechanisms probability, in degrees (e.g., 45)
	prob_max = 0.1 # Enter probability threshold for multiples (e.g., 0.1)
	#----------------------------------------------------------------
	
	# Mark's spec'd object variables
	vmodels = []
	vtables = None
	ntab = 0
	
	dist = None
	qazi = None
	flat = None
	flon = None
	felv = None
	
	# polarity reversals, [-1,1] stub for now
	spol = 1
	
	# defaults, set for now...
	npol = 0 
	
	def __init__(self, **kwargs):
		'''Make one, empty or pass to other fxns.'''
		npick0, nmc0, nmax0 = fortran_include('param.inc')
		dang0, ncoor        = fortran_include('rot.inc')
		
		# initialize arrays
		self.dang2=max(dang0, self.dang)
		
		# Input arrays
		self.sname     = np.empty(npick0, 'a4', 'F')
		self.scomp     = np.empty(npick0, 'a3', 'F')
		self.snet      = np.empty(npick0, 'a2', 'F')
		self.pickpol   = np.empty(npick0, 'a1', 'F')
		self.pickonset = np.empty(npick0, 'a1', 'F')
		self.p_pol     = np.empty(npick0, int, 'F')
		self.p_qual    = np.empty(npick0, int, 'F')
		self.spol      = np.empty(npick0, int, 'F')
		self.p_azi_mc  = np.empty((npick0,nmc0), float, 'F')
		self.p_the_mc  = np.empty((npick0,nmc0), float, 'F')
		self.index     = np.empty(nmc0, int, 'F')
		self.qdep2     = np.empty(nmc0, float, 'F')
		
		# Output arrays
		self.f1norm  = np.empty((3,nmax0), float, 'F')
		self.f2norm  = np.empty((3,nmax0), float, 'F')
		self.strike2 = np.empty(nmax0, float, 'F')
		self.dip2    = np.empty(nmax0, float, 'F')
		self.rake2   = np.empty(nmax0, float, 'F')
		self.str_avg = np.empty(5, float, 'F')
		self.dip_avg = np.empty(5, float, 'F')
		self.rak_avg = np.empty(5, float, 'F')
		self.var_est = np.empty((2,5), float, 'F')
		self.var_avg = np.empty(5, float, 'F')
		self.mfrac   = np.empty(5, float, 'F')
		self.stdr    = np.empty(5, float, 'F')
		self.prob    = np.empty(5, float, 'F')
		self.qual    = np.empty(5, 'a', 'F')
		
		# Other (added for Python classiness)
		self.dist    = np.empty(npick0, float)
		self.qazi    = np.empty(npick0, float)
		self.flat    = np.empty(npick0, float)
		self.flon    = np.empty(npick0, float)
		self.felv    = np.empty(npick0, float)
		self.esaz    = np.empty(npick0, float)
		
		# Save include vars for other fucntions to access
		self.npick0 = npick0
		self.nmc0   = nmc0
		self.nmax0  = nmax0
		self.dang0  = dang0
		self.ncoor  = ncoor
		
			
		if kwargs:
			self.__dict__.update(kwargs)
		
	def load_pf(self, pffile='dbhash.pf'):
		'''update some run settings from a pf file'''
		pass
		
	def load_velocity_models(self, model_list=None):
		'''load velocity model data'''
		# Future -- allow adding on fly, check and append to existing
		#
		# THIS LOADS INTO hashpy.angtable.table (nx,nd,nindex) !!!
		if model_list:
			self.vmodels = model_list
		for n,v in enumerate(self.vmodels):
			self.ntab = mk_table_add(n+1,v)
			self.vtable = angtable.table

	def get_phases_from_db(self, dbname, evid=None, orid=None, pf=None):
		db, oflag = open_db_or_string(dbname)
		if orid is None:
			dbv = db.process(['dbopen event', 'dbsubset evid == '+str(evid)])
			dbv.record = 0
			orid = dbv.getv('prefor')[0]
		db = db.process([ 'dbopen origin', 'dbsubset orid == '+str(orid),
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
		self.seh    = ph['origerr.smajax']
		self.sez    = ph['origerr.sdepth']
		
		aspect = np.cos(self.qlat / degrad) # convert using python later.
		
		# The index 'k' is deliberately non-Pythonic to deal with the fortran
		# subroutines which need to be called and the structure of the original HASH code.
		# May be able to update with a rewrite... YMMV
		k = 0
		for ph in phases:
			# load up params
			# in future, could use the acol() method?
			self.sname[k]     = ph.sta
			self.snet[k]      = ph.net
			self.scomp[k]     = ph.chan
			self.pickonset[k] = 'I'
			self.pickpol[k]   = ph.fm
			
			flat,flon,felv = ph['site.lat'],ph['site.lon'],ph['site.elev']
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
		
	def generate_trial_data(self):
		'''Make data for running trials MUST have loaded data and vel mods alreday
		
		Algorithm NOT written by me (From HASH driver script)
		'''
		# choose a new source location and velocity model for each trial
		self.qdep2[0] = self.qdep
		self.index[0] = 1
		for nm in range(1,self.nmc):
			val = ran_norm()
			self.qdep2[nm] = abs(self.qdep + self.sez * val) # randomly perturbed source depth
			self.index[nm] = (nm % self.ntab) + 1  # index used to choose velocity model
			
	def calculate_takeoff_angles(self):
		'''Use HASH fortran subroutine to calulate takeoff angles for each trial'''
		# loop over k picks b/c I broke it out -- NOTE: what does iflag do?
		#
		# find azimuth and takeoff angle for each trial
		for k in range(self.npol):
			for nm in range(self.nmc):
				self.p_azi_mc[k,nm] = self.qazi[k]
				self.p_the_mc[k,nm], iflag = get_tts(self.index[nm],self.dist[k],self.qdep2[nm])

	def view_polarity_data(self):
		'''Print out a list of polarity data for interactive runs'''
		for k in range(self.npol):
			print '{0}   {1} {2} {3} {4}'.format(k,self.sname[k],self.p_azi_mc[k,0],self.p_the_mc[k,0],self.p_pol[k])
	
	def check_minimum_polaritiy(self):
		if self.npol >= self.npolmin:
			return True
		else:
			return False
	
	def check_maximum_gap(self):
		magap,mpgap = get_gap(self.p_azi_mc[:self.npol,0],self.p_the_mc[:self.npol,0],self.npol)
		if ((magap > self.max_agap) or (mpgap > self.max_pgap)):
			return False
		else:
			return True

	def calculate_hash_focalmech(self):
		# determine maximum acceptable number misfit polarities
		nmismax = max(int(self.npol * self.badfrac),2)        # nint
		nextra  = max(int(self.npol * self.badfrac * 0.5),2)  # nint
		
		# find the set of acceptable focal mechanisms for all trials
		nf2, self.strike2, self.dip2, self.rake2, self.f1norm, self.f2norm = focalmc(self.p_azi_mc, self.p_the_mc, self.p_pol[:self.npol], self.p_qual[:self.npol], self.nmc, self.dang2, self.nmax0, nextra, nmismax, self.npol)
		self.nout2 = min(self.nmax0,nf2)  # number mechs returned from sub
		self.nout1 = min(self.maxout,nf2) # number mechs to return
		
		# find the probable mechanism from the set of acceptable solutions
		nmult, self.str_avg, self.dip_avg, self.rak_avg, self.prob, self.var_est = mech_prob(self.f1norm[:,:self.nout2], self.f2norm[:,:self.nout2], self.cangle, self.prob_max, self.nout2) # nout2
		
		for imult in range(nmult):
			self.var_avg[imult] = (self.var_est[0,imult] + self.var_est[1,imult]) / 2.
			#print 'cid = {0} {1}  mech = {2} {3} {4}'.format(self.icusp, imult, self.str_avg[imult], self.dip_avg[imult], self.rak_avg[imult])
			# find misfit for prefered solution
			self.mfrac[imult], self.stdr[imult] = get_misf(self.p_azi_mc[:self.npol,0], self.p_the_mc[:self.npol,0], self.p_pol[:self.npol], self.p_qual[:self.npol], self.str_avg[imult], self.dip_avg[imult], self.rak_avg[imult], self.npol) # npol
			
			# HASH default solution quality rating
			if ((self.prob[imult] > 0.8) and (self.var_avg[imult] < 25) and (self.mfrac[imult] <= 0.15) and (self.stdr[imult] >= 0.5)):
				self.qual[imult]='A'
			elif ((self.prob[imult] > 0.6) and (self.var_avg[imult] <= 35) and (self.mfrac[imult] <= 0.2) and (self.stdr[imult] >= 0.4)):
				self.qual[imult]='B'
			elif ((self.prob[imult] > 0.5) and (self.var_avg[imult] <= 45) and (self.mfrac[imult] <= 0.3) and (self.stdr[imult] >= 0.3)):
				self.qual[imult]='C'
			else:
				self.qual[imult]='D'


	def plot_beachball(self):
		'''
		test_stereo(self.p_azi_mc[:self.npol,0], self.p_the_mc[:self.npol,0], self.p_pol[:self.npol], sdr=[self.str_avg[0], self.dip_avg[0], self.rak_avg[0]])
		azimuths,takeoffs,polarities,sdr=[]
		'''
		import matplotlib.pyplot as plt
		import mplstereonet
		from obspy.imaging.beachball import AuxPlane
		
		# pull out variables from mechanism
		azimuths = self.p_azi_mc[:self.npol,0]
		takeoffs = self.p_the_mc[:self.npol,0]
		polarities = self.p_pol[:self.npol]
		strike1,dip1,rake1 = self.str_avg[0], self.dip_avg[0], self.rak_avg[0]
		strike2,dip2,rake2 = AuxPlane(strike1, dip1, rake1)
		up = polarities > 0
		dn = polarities < 0
		# plot stations, trial planes, avg planes
		fig = plt.figure()
		ax = fig.add_subplot(111, projection='stereonet')
		for n in range(self.nout1):
			s1, d1, r1 = self.strike2[n], self.dip2[n], self.rake2[n]
			s2, d2, r2 = AuxPlane(s1, d1, r1)
			h_rk = ax.plane(s1,d1, color='#999999')
			#h_rk = ax.plane(s2,d2,'#888888')
		h_rk = ax.plane(strike1, dip1, color='black', linewidth=3)
		h_rk = ax.rake( strike1, dip1, -rake1, 'ko')
		h_rk = ax.plane(strike2, dip2, color='black', linewidth=3)
		#h_rk = ax.rake(azimuths[up]+90.,takeoffs[up],90, 'wo', markersize=8, markeredgewidth=2, markerfacecolor=None)
		#h_rk = ax.rake(azimuths[dn]+90.,takeoffs[dn],90, 'r+', markersize=8, markeredgewidth=2)
		#for i in range(self.npol):
		#	h_rk = ax.rake(azimuths[i]+90,takeoffs[i]+5,90, marker='$   {0}$'.format(self.sname[i]), color='black',markersize=20)
		plt.show()
#------------------------------------------
# for testing:

test_orid = 933757 # 2012 206 (shot SPE3)

c_evid = '''
379439 2012167
380226 2012174
380653 2012175
380828 2012175
'''

test_table_list =  ['vz.socal',
					'vz.north',
					'vz.lab1',
					'vz.sgm1',
					'vz.vb1']

test_table_list = ['vz.pickema1', 'vz.pickema2', 'vz.pickema3']

foo = FocalMech(maxout=100)
foo.get_phases_from_db('/data/2012/206/reno', orid=test_orid)
foo.load_velocity_models(test_table_list)
foo.generate_trial_data()
foo.calculate_takeoff_angles()
foo.calculate_hash_focalmech()
foo.plot_beachball()

# Code to make a quick and dirty station map
if False:
	figure()
	plot(foo.flon[:foo.npol],foo.flat[:foo.npol],'o')
	for i in range(foo.npol):
		text(foo.flon[i],foo.flat[i],foo.sname[i])
	plot(foo.qlon, foo.qlat, 'ro')
