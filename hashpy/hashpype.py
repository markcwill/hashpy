# -*- coding: utf-8 -*-
#
#  hashpype.py
#
# by Mark Williams 2012.313
# First motion focal mechanism classes for running HASH

# import all Fortran subroutines and common blocks from the mod
# plus my custom utils and numpy arrays

import numpy as np
import os
from pwd import getpwuid
from libhashpy import (mk_table_add, angtable, ran_norm, get_tts, get_gap,
	focalmc, mech_prob, get_misf,)

def parameter(**kwargs):
	"""
	Returns variables inside a fortran 'parameter' call
	"""
	return kwargs

def fortran_include(fname):
	"""
	Fortran include
	
	Quickie for access to 'include' files, needed for access to
	parameters useful for interacting with f2py wrapped functions
	
	Input
	-----
	fname : str of fortran include filename
	
	Returns : dict of parameters described in any 'parameter' call
	"""
	inc_params = {}
	f = open(fname)
	for line in f:
		if 'c' in line[0]:
			continue
		elif 'parameter' in line:
			# evaluate the line 'parameter(foo=999)' as python expr
			inc_params.update(eval(line))
		else:
			pass
	f.close()
	return inc_params
	
	
class HashPype(object):
	"""
	Object which holds all data from a HASH instance for one event
	
	Methods are based on the 'hash_driver2.f' program from H & S
	
	The variables are named and structured as close to the original code
	as possible. Methods act on the variables within this namespace and
	most are accessible as attributes.
	
	One can make a "HASH" driver program by calling the methods of this
	class on the data held within it, but the input data must be entered
	somehow. There are no input methods, and one very simple output
	method. Why?
	
	The idea is to use this as a metaclass which inherits all these
	methods and attributes to another class which defines various
	input/output for whatever system is used.
	
	*See the EventHashPype or DbHashPype class in the program for details.
	
	"""
	# These MUST be the same as the fortran includes!!
	# (They are compiled into the Fortran subroutines)
	#--- param.inc ---#
	npick0, nmc0, nmax0 = None, None, None
	#--- rot.inc ---#
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
	prob_max = 0.1  # Enter probability threshold for multiples (e.g., 0.1)
	#----------------------------------------------------------------
	
	# Mark's spec'd object variables
	vmodel_dir = None # string of directory containing velocity models
	vmodels = []      # list of string names of velocity model files
	vtables = None    # actual travel time tables in common block
	arid    = None    # unique number for each pick...
	author  = None    # logged in user
	
	
	# Variables HASH keeps internally, for ref and passing to fxns
	ntab = 0     # number of tables loaded
	npol = 0     # number of polarity picks
	dist = None  # distance from source im km
	qazi = None  # azimuth from event to station
	flat = None  # pick station lat
	flon = None  # pick station lon
	felv = None  # pick station elv
	nf2 = None   # number of fm's returned
	nmult = None # number of fm solutions returned
	magap = None # calulated max azi gap
	mpgap = None # calculated max plunge gap (I think)
	
	tstamp  = None
	qlat    = None
	qlon    = None
	qdep    = None
	qmag    = None
	icusp   = None
	seh     = None
	sez     = None
	
	
	# polarity reversals, [-1,1] stub for now
	spol = 1
	
	def __init__(self, **kwargs):
		"""
		Make an empty HASH run object.
		
		Initialize all the arrays needed for a HASH run. They are based
		on the maximum size of the arrays passed to the FORTRAN
		subroutines. 
		
		ANYTHING initialized here can be changed by passing
		a keyword arg at the Constructor call. Useful for using an
		input parameter file for HASH.
		
		Example
		-------
		>>> h = HashPype()
		"""
		
		# INCLUDES ------------------------------------------------
		# contain array sizes, get em
		directory = os.path.dirname(__file__)
		param_inc_file = os.path.join(directory,'src','param.inc')
		rot_inc_file   = os.path.join(directory,'src','rot.inc')
		
		param_inc = fortran_include(param_inc_file) # npick0, nmc0, nmax0 
		rot_inc = fortran_include(rot_inc_file) # dang0, ncoor
		
		# Save include vars for other fucntions to access
		self.__dict__.update(param_inc)
		self.__dict__.update(rot_inc)
		
		# For convenience, throw in fxn namespace
		npick0 = self.npick0
		nmc0   = self.nmc0
		nmax0  = self.nmax0
		dang0  = self.dang0
		ncoor  = self.ncoor
		
		# ARRAYS ---------------------------------------------------
		# initialize arrays for HASH
		self.dang2=max(self.dang0, self.dang)
		
		# Input arrays
		self.sname     = np.empty(npick0, 'a6', 'F')
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
		
		# Internals (added for Python classiness)
		self.dist    = np.empty(npick0, float)
		self.qazi    = np.empty(npick0, float)
		self.flat    = np.empty(npick0, float)
		self.flon    = np.empty(npick0, float)
		self.felv    = np.empty(npick0, float)
		self.esaz    = np.empty(npick0, float)
		self.arid  	 = np.empty(npick0, int) * 0
		
		self.author = getpwuid(os.getuid()).pw_name
		
		if kwargs:
			self.__dict__.update(kwargs)
	
	def __repr__(self):
		"""
		String saying what you are
		"""
		return '{0}(icusp={1})'.format(self.__class__.__name__, self.icusp)
	
	def load_velocity_models(self, model_list=None):
		"""
		"Load velocity model data
		
		Inputs
		------
		model_list : list of str of velocity model filenames
		
		** if None, will use the list in 'self.vmodels'
		"""
		# Future -- allow adding on fly, check and append to existing
		#
		# THIS LOADS INTO hashpy.angtable.table (nx,nd,nindex) !!!
		
		# take care of 
		if model_list:
			models = model_list
		else:
			models = self.vmodels
		
		for n,v in enumerate(models):
			self.ntab = mk_table_add(n+1,v)
			self.vtable = angtable.table
		
	def generate_trial_data(self):
		"""
		Make data for running trials
		(MUST have loaded data and vel mods already)
		
		Algorithm by H & S (From HASH driver script)
		"""
		# choose a new source location and velocity model for each trial
		self.qdep2[0] = self.qdep
		self.index[0] = 1
		for nm in range(1,self.nmc):
			val = ran_norm()
			self.qdep2[nm] = abs(self.qdep + self.sez * val) # randomly perturbed source depth
			self.index[nm] = (nm % self.ntab) + 1  # index used to choose velocity model
			
	def calculate_takeoff_angles(self):
		"""
		Use HASH fortran subroutine to calulate takeoff angles for each trial
		"""
		# loop over k picks b/c I broke it out -- NOTE: what does iflag do?
		#
		# find azimuth and takeoff angle for each trial
		for k in range(self.npol):
			for nm in range(self.nmc):
				self.p_azi_mc[k,nm] = self.qazi[k]
				self.p_the_mc[k,nm], iflag = get_tts(self.index[nm],self.dist[k],self.qdep2[nm])
		self.magap, self.mpgap = get_gap(self.p_azi_mc[:self.npol,0],self.p_the_mc[:self.npol,0],self.npol)
	
	def view_polarity_data(self):
		"""
		Print out a list of polarity data for interactive runs
		"""
		for k in range(self.npol):
			print '{0}   {1} {2} {3} {4}'.format(k,self.sname[k],self.p_azi_mc[k,0],self.p_the_mc[k,0],self.p_pol[k])
	
	def check_minimum_polarity(self):
		"""Polarity check"""
		if self.npol >= self.npolmin:
			return True
		else:
			return False
	
	def check_maximum_gap(self):
		"""Gap check"""
		if ((self.magap > self.max_agap) or (self.mpgap > self.max_pgap)):
			return False
		else:
			return True
	
	def calculate_hash_focalmech(self):
		"""
		Run the actual focal mech calculation, and find the probable mech
		"""
		# determine maximum acceptable number misfit polarities
		nmismax = max(int(self.npol * self.badfrac),2)        # nint
		nextra  = max(int(self.npol * self.badfrac * 0.5),2)  # nint
		
		# find the set of acceptable focal mechanisms for all trials
		self.nf2, self.strike2, self.dip2, self.rake2, self.f1norm, self.f2norm = focalmc(self.p_azi_mc, self.p_the_mc, self.p_pol[:self.npol], self.p_qual[:self.npol], self.nmc, self.dang2, self.nmax0, nextra, nmismax, self.npol)
		self.nout2 = min(self.nmax0, self.nf2)  # number mechs returned from sub
		self.nout1 = min(self.maxout, self.nf2) # number mechs to return
		
		# find the probable mechanism from the set of acceptable solutions
		self.nmult, self.str_avg, self.dip_avg, self.rak_avg, self.prob, self.var_est = mech_prob(self.f1norm[:,:self.nout2], self.f2norm[:,:self.nout2], self.cangle, self.prob_max, self.nout2) # nout2
	
	def calculate_quality(self):
		"""
		Do the quality calulations
		"""
		for imult in range(self.nmult):
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

class HashError(Exception):
	"""Throw this if something happens while running"""
	pass





