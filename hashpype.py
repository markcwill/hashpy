#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  hashpype.py

from sys import argv
from hashpy import *
import numpy as np
from hash_utils import fortran_include, get_sta_coords, test_stereo
from ant_utils import add_antelope_path, open_db_or_string
from obspy_ext.antelope import *
add_antelope_path()

test_orid = 946268


class FocalMech(object):
	'''First motion focal mechanisms'''
	# These MUST be the same as the fortran includes!!
	# (They are compiled into the Fortran subroutines)
	npick0, nmc0, nmax0 = None, None, None
	dang0, ncoor        = None, None
	
	# initialize arrays
	
	# Input arrays
	sname     = np.empty(npick0, 'a4', 'F')
	scomp     = np.empty(npick0, 'a3', 'F')
	snet      = np.empty(npick0, 'a2', 'F')
	pickpol   = np.empty(npick0, 'a1', 'F')
	pickonset = np.empty(npick0, 'a1', 'F')
	p_pol     = np.empty(npick0, int, 'F')
	p_qual    = np.empty(npick0, int, 'F')
	spol      = np.empty(npick0, int, 'F')
	p_azi_mc  = np.empty((npick0,nmc0), float, 'F')
	p_the_mc  = np.empty((npick0,nmc0), float, 'F')
	index     = np.empty(nmc0, int, 'F')
	qdep2     = np.empty(nmc0, float, 'F')
	
	# Output arrays
	f1norm  = np.empty((3,nmax0), float, 'F')
	f2norm  = np.empty((3,nmax0), float, 'F')
	strike2 = np.empty(nmax0, float, 'F')
	dip2    = np.empty(nmax0, float, 'F')
	rake2   = np.empty(nmax0, float, 'F')
	str_avg = np.empty(5, float, 'F')
	dip_avg = np.empty(5, float, 'F')
	rak_avg = np.empty(5, float, 'F')
	var_est = np.empty((2,5), float, 'F')
	var_avg = np.empty(5, float, 'F')
	mfrac   = np.empty(5, float, 'F')
	stdr    = np.empty(5, float, 'F')
	prob    = np.empty(5, float, 'F')
	qual    = np.empty(5, 'a', 'F')
	
	degrad = 180. / np.pi
	rad = 1. / degrad
	spol = 1
	
	# defaults, set for now...
	npol = 0
	
	def __init__(self, **kwargs):
		'''Make one, empty or pass to other fxns.'''
		npick0, nmc0, nmax0 = fortran_include('param.inc')
		dang0, ncoor        = fortran_include('rot.inc')
	
	def get_data(self):
		'''Import the data from a db given some params'''

	def get_data_from_db(self, dbname, evid=None, orid=None, pf=None):
		db, oflag = open_db_or_string(dbname)
		if orid is None:
			dbv = db.process(['dbopen event', 'dbsubset evid == '+str(evid)])
			orid = dbv.getv('prefor')[0]
		db = db.process([ 'dbopen origin', 'dbsubset orid == '+str(orid),
						'dbjoin origerr', 'dbjoin assoc',  'dbjoin arrival',
						'dbjoin affiliation', 'dbjoin site',
						'dbsubset iphase =~ /.*[Pp].*/',
						'dbsubset (ondate <= time)',
						'dbsubset (time <= offdate) || (offdate == -1)']
						)
		
		phases = AttribDbptr(db)
		ph = phases[0]
		tstamp = ph['origin.time']
		qlat   = ph['origin.lat']
		qlon   = ph['origin.lon']
		qdep   = ph['origin.depth']
		qmag   = ph['origin.ml']
		icusp  = ph['origin.orid']
		seh    = ph['origerr.smajax']
		sez    = ph['origerr.sdepth']
		
		aspect = np.cos(qlat / degrad) # convert using python later.
		
		# choose a new source location and velocity model for each trial
		qdep2[0] = qdep
		index[0] = 1
		for nm in range(1,nmc):
			val = ran_norm()
			qdep2[nm] = qdep + sez * val # randomly perturbed source depth
			index[nm] = (nm % ntab) + 1  # index used to choose velocity model
		
		for k, ph in enumerate(phases):
			# load up params
			# in future, could use the acol() method?
			sname[k]     = ph.sta
			snet[k]      = ph.net
			scomp[k]     = ph.chan
			pickonset[k] = 'I'
			pickpol[k]   = ph.fm
			
			flat,flon,felv = ph['site.lat'],ph['site.lon'],ph['site.elev']
			#print '{0} {1} {2} {3} {4} {5} {6} {7}'.format(k, sname[k], snet[k], scomp[k], pickonset[k], pickpol[k], flat, flon)
			
			# dist @ azi, get from db OR obspy or another python mod (antelope) could do this on WGS84
			dx = (flon - qlon) * 111.2 * aspect
			dy = (flat - qlat) * 111.2
			dist = np.sqrt(dx**2 + dy**2)
			qazi = 90. - np.arctan2(dy,dx) * degrad
			
			if (qazi < 0.):
				qazi = qazi + 360.
			if (dist > delmax):
				continue
			if (pickpol[k] in 'CcUu+'):
				p_pol[k] = 1
			elif (pickpol[k] in 'RrDd-'):
				p_pol[k] = -1
			else:
				continue
			
			if (pickonset[k] in 'Ii'):
				p_qual[k] = 0
			else:
				p_qual[k] = 1
			
			#spol = check_pol(plfile,sname[k],iyr,imon,idy,ihr)  # SCSN station polarity reversal information - ** YOU MAY NEED TO USE YOUR OWN SUBROUTINE **
			p_pol[k] *= spol
			
			# find azimuth and takeoff angle for each trial
			for nm in range(nmc):
				p_azi_mc[k,nm] = qazi
				p_the_mc[k,nm], iflag = get_tts(index[nm],dist,qdep2[nm])
				
			#k += 1
			#continue
		#npol = k - 1
		self.npol = k
	
	def __repr__(self):
		for k in range(npol):
			print '{0}   {1} {2} {3} {4}'.format(k,sname[k],p_azi_mc[k,0],p_the_mc[k,0],p_pol[k])


