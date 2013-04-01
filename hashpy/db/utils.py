# -*- coding: utf-8 -*-
#
# utils.py
#
# by Mark Williams 2012.313
#
# Utilities for using Antelope db with focal mechanisms and obspy
#


import sys, os
import numpy as np
from obspy.core.util import AttribDict
from obspy.imaging.beachball import AuxPlane
from obspy.core import read, Stream, UTCDateTime
utc = UTCDateTime
array = np.array

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
	
	_plane = None
	
	@property
	def plane1(self):
		'''Return Preferred plane'''
		return NodalPlane(*self._plane)
	
	@property
	def plane2(self):
		'''Return Auxiliary plane'''
		auxplane = AuxPlane(*self._plane)
		return NodalPlane(*auxplane)
	
	@property
	def axis(self):
		'''return direction and dip for P and T axes'''
		dipP, dipT, aziP, aziT = nodal2pt(*self.plane1+self.plane2)
		return AttribDict({'P': AttribDict({'azi': aziP, 'dip': dipP}),
						   'T': AttribDict({'azi': aziT, 'dip': dipT}) })
	
	def __init__(self, nodal_plane=None):
		self._plane = nodal_plane



def add_antelope_path():
	_version_string = os.environ['ANTELOPE'].split('/')[-1]
	_pydirs = ['data','python']
	if float(_version_string[:3]) < 5.2:
		_pydirs = ['local'] + _pydirs
	_pypath = os.path.join(os.environ['ANTELOPE'], *_pydirs)
	if _pypath not in sys.path:
		sys.path.append(_pypath)
add_antelope_path()

from aug.contrib import *
from antelope.datascope import *
from antelope.stock import *

def load_pf(pffile='dbhash.pf'):
	'''Load HASH runtime settings from a pf file to a dictionary
	
	One can also specify names of velocity model files in the pf.
	
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
		pf_settings[key] = pfi
	
	if 'vmodel_dir' in pf_settings and 'vmodels' in pf_settings:
		pf_settings['vmodels'] = [os.path.join(pf_settings['vmodel_dir'], table) for table in pf_settings['vmodels']]
	return pf_settings

def db2object(dbv):
	'''
	Port of Antelope MATLAB toolbox 'db2struct' function.
		
	Returns a list-like object, this is the function version of calling
	DbrecordList() directly.
	
	:type dbv: antelope.datascope.Dbptr
	:param dbv: Open pointer to an Antelope database view or table
	:rtype: :class:`~obspy.antelope.Dbview`
	:return: Dbview of Dbrecord objeccts
	'''
	if isinstance(dbv, Dbptr):
		db = Dbptr(dbv)
	else:
		raise TypeError("'{0}' is not a Dbptr object".format(dbv))
	return DbrecordList(db)

	
def readANTELOPE(database, station=None, channel=None, starttime=None, endtime=None):
	'''
	Reads a portion of a Antelope wfdisc table to a Stream.
	
	Attempts to return one Trace per line of the 'wfdisc' view passed.	
	Additionally, will filter and cut with respect to any of the fields
	in the primary key IF specified. (sta chan time::endtime)
	
	NOTE: Currently MUST have both times (start/end) or neither.
	the returned Traces will have a new attribute, 'db'

	:type database: string or antelope.datascope.Dbptr
	:param database: Antelope database name or pointer
	:type station: string
	:param station: Station expression to subset
	:type channel: string
	:param channel: Channel expression to subset
	:type starttime: :class: `~obspy.core.utcdatetime.UTCDateTime`
	:param starttime: Desired start time
	:type endtime: :class: `~obspy.core.utcdatetime.UTCDateTime`
	:param endtime: Desired end time
		
	:rtype: :class: `~obspy.core.stream.Stream'
	:return: Stream with one Trace for each row of the database view
	
	.. rubric:: Example
	
	>>> st = readANTELOPE('/Volumes/colza_HD/dbs/land', station='TOL0', channel='LH.',
						starttime=UTCDateTime(2008,6,13), endtime=UTCDateTime(2008,6,14))
	>>> print(st)
	6 Trace(s) in Stream:
	XA.TOL0..LHE | 2008-06-12T23:59:59.640000Z - 2008-06-13T00:04:11.640000Z | 1.0 Hz, 253 samples
	XA.TOL0..LHE | 2008-06-13T00:04:12.640000Z - 2008-06-13T23:59:59.640000Z | 1.0 Hz, 86148 samples
	XA.TOL0..LHN | 2008-06-12T23:59:59.640000Z - 2008-06-13T00:04:11.640000Z | 1.0 Hz, 253 samples
	XA.TOL0..LHN | 2008-06-13T00:04:12.640000Z - 2008-06-13T23:59:59.640000Z | 1.0 Hz, 86148 samples
	XA.TOL0..LHZ | 2008-06-12T23:59:59.640000Z - 2008-06-13T00:04:21.640000Z | 1.0 Hz, 263 samples
	XA.TOL0..LHZ | 2008-06-13T00:04:22.640000Z - 2008-06-13T23:59:59.640000Z | 1.0 Hz, 86138 samples
	
	Also adds a Dbrecord as an attribute of the Trace
	
	>>> st[0].db
	Dbrecord('View43' -> TOL0 LHE 1213229044.64::1213315451.64)
 
	'''
	if isinstance(database,Dbptr):
		db = Dbptr(database)
	elif isinstance(database,str):
		db = dbopen(database, 'r')
		db = dblookup(db,table='wfdisc')
	else:
		raise TypeError("Must input a string or pointer to a valid database")
		
	if station is not None:
		db = dbsubset(db,'sta=~/{0}/'.format(station))
	if channel is not None:
		db = dbsubset(db,'chan=~/{0}/'.format(channel))
	if starttime is not None and endtime is not None:
		ts = starttime.timestamp
		te = endtime.timestamp
		db = dbsubset(db,'endtime > {0} && time < {1}'.format(ts,te) )
	else:
		ts = starttime
		te = endtime
	assert db.nrecs() is not 0, "No records for given time period"
	
	st = Stream()
	for db.record in range(db.nrecs() ):
		fname = db.filename() 
		dbr = DbrecordPtr(db)
		t0 = UTCDateTime(dbr.time)
		t1 = UTCDateTime(dbr.endtime)
		if dbr.time < ts:
			t0 = starttime
		if dbr.endtime > te:
			t1 = endtime
		_st = read(fname, starttime=t0, endtime=t1)		 # add format?
		_st = _st.select(station=dbr.sta, channel=dbr.chan) #not location aware
		#_st[0].db = dbr
		if dbr.calib < 0:
			_st[0].data *= -1
		st += _st
	# Close what we opened, BUT garbage collection may take care of this:
	# if you have an open pointer but pass db name as a string, global
	# use of your pointer won't work if this is uncommented:
	#
	#if isinstance(database,str):
	#	db.close()
	return st

def dbloc_source_db(db):
	'''Checks if you are in a dbloc2 'trial' db and returns the source
	one if you are, otherwise returns the same Dbptr
	
	INPUT: Dbptr
	OUTPUT: Dbptr to database that dbloc2 is using.
	'''
	dbname = db.query(dbDATABASE_NAME)
	pfdef = pfget('dbloc2','Define')
	tempdb = pfdef['Temporary_db']
	workdir = pfdef['Work_dir']
	dblocdb = os.path.join(workdir,tempdb)
	if dbname.endswith(tempdb):
		# path of trial db from dbloc2
		dbcwd = os.path.dirname(dbname)
		# relative name of 1st db in 'trial' database decriptor file
		dbpath0 = db.query(dbDBPATH).split(':')[0].translate(None,'{}')
		# full absolute path database name to source
		realdb = os.path.abspath(os.path.join(dbcwd, dbpath0))
		db.close()
		db = dbopen(realdb, perm='r+')
	return db

def focalmech2db(focalmech):
	'''Write the preferred HASH solution to Datascope database.
	
	Writes to 'fplane', 'predmech' and 'predarr' tables
	'''
	fp = focalmech
	axes = fp.axis
	T = axes['T']
	P = axes['P']
	
	database = focalmech.source
	db, oflag = open_db_or_string(database, perm='r+')
	# Use the original db if in a dbloc2 'tmp/trial' db
	db = dbloc_source_db(db)
	
	mechid = db.nextid('mechid')
	
	dbfpln = dblookup(db,table='fplane')
	dbfpln.record = dbfpln.addnull()
	dbfpln.putv('orid', fp.orid,
		'str1', round(fp.plane1.strike,1) ,
		'dip1', round(fp.plane1.dip,1) ,
		'rake1',round(fp.plane1.rake,1),
		'str2', round(fp.plane2.strike,1) ,
		'dip2', round(fp.plane2.dip,1) ,
		'rake2',round(fp.plane2.rake,1),
		'taxazm',round(T.azi,1),
		'taxplg',round(T.dip,1),
		'paxazm',round(P.azi,1),
		'paxplg',round(P.dip,1),
		'algorithm', fp.algorithm,
		'mechid', mechid
		)
	dbpmec = dblookup(db,table='predmech')
	dbparr = dblookup(db,table='predarr')
	for k in range(len(fp.picks)):
		pk = fp.picks[k]
		if pk['polarity'] > 0:
			fm = 'U'
		else:
			fm = 'D'
		dbpmec.record = dbpmec.addnull()
		dbpmec.putv('arid', int(pk['arid']) ,
					'orid', fp.orid,
					'mechid', mechid,
					'fm', fm,
					)
		dbparr.record = dbparr.addnull()
		dbparr.putv('arid', int(pk['arid']),
					'orid', fp.orid, 
					'esaz', pk['azimuth'], 
					'dip' , pk['takeoff'],
					)
	db.close()

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
	
def change_arrival_fm(database, arid, new_fm):
	'''Change the first motion string of a pick in the arrival table'''
	db,o = open_db_or_string(database, perm='r+')
	dbv = dbprocess(db, ['dbopen arrival', 'dbsubset arid=={0}'.format(arid)])
	dbv.record = 0
	d = DbrecordPtr(dbv)
	lddate = d.lddate
	d.fm = new_fm
	d.lddate = lddate
	db.close()
