# -*- coding: utf-8 -*-
#
# utils.py
#
# by Mark Williams 2012.313
#
# Utilities for using Antelope db with focal mechanisms and obspy
#
import sys, os
from obspy.core import read, Stream, UTCDateTime


class WhiteRumpError(Exception):
	'''Call this if you have a problem you want to catch from here'''
	pass

def add_antelope_path():
	_version_string = os.environ['ANTELOPE'].split('/')[-1]
	_pydirs = ['data','python']
	if float(_version_string[:3]) < 5.2:
		_pydirs = ['local'] + _pydirs
	_pypath = os.path.join(os.environ['ANTELOPE'], *_pydirs)
	if _pypath not in sys.path:
		sys.path.append(_pypath)

add_antelope_path()

from aug.contrib.orm import open_db_or_string, AttribDbptr, DbrecordPtr
from antelope.datascope import (Dbptr, dbopen, dblookup, dbsubset, dblookup,
	dbprocess, dbDATABASE_NAME, dbDBPATH)


def readANTELOPE(database, station=None, channel=None, starttime=None, endtime=None):
	"""
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
 
	"""
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
		if os.path.exists(fname):
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
	"""
	Checks if you are in a dbloc2 'trial' db and returns the source
	one if you are, otherwise returns the same Dbptr
	
	INPUT: Dbptr
	OUTPUT: Dbptr to database that dbloc2 is using.
	"""
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

def eventfocalmech2db(event=None, database=None):
	"""
	Write the preferred HASH solution to Datascope database.
	
	Writes to 'fplane', 'predmech' and 'predarr' tables
	"""
	focm = event.preferred_focal_mechanism()
	o = focm.triggering_origin_id.getReferredObject()
	
	plane1 = focm.nodal_planes.nodal_plane_1
	plane2 = focm.nodal_planes.nodal_plane_2
	T = focm.principal_axes.t_axis
	P = focm.principal_axes.p_axis
	orid = int(o.creation_info.version)
	
	db, oflag = open_db_or_string(database, perm='r+')
	try:
		# Use the original db if in a dbloc2 'tmp/trial' db
		db = dbloc_source_db(db)
		# save solution as a new mechid
		mechid = db.nextid('mechid')
		# in fplane...
		dbfpln = dblookup(db,table='fplane')
		dbfpln.record = dbfpln.addnull()
		dbfpln.putv('orid', orid,
			'str1', round(plane1.strike,1) ,
			'dip1', round(plane1.dip,1) ,
			'rake1',round(plane1.rake,1),
			'str2', round(plane2.strike,1) ,
			'dip2', round(plane2.dip,1) ,
			'rake2',round(plane2.rake,1),
			'taxazm',round(T.azimuth,1),
			'taxplg',round(T.plunge,1),
			'paxazm',round(P.azimuth,1),
			'paxplg',round(P.plunge,1),
			'algorithm', focm.method_id.resource_id,
			'auth', focm.creation_info.author,
			'mechid', mechid,
			)
		dbpmec = dblookup(db,table='predmech')
		dbparr = dblookup(db,table='predarr')
		for av in o.arrivals:
			pk = av.pick_id.getReferredObject()
			if pk.polarity is 'positive':
				fm = 'U'
			elif pk.polarity is 'negative':
				fm = 'D'
			else:
				continue
			
			arid = int(av.creation_info.version)
			
			# ..and predmech
			dbpmec.record = dbpmec.addnull()
			dbpmec.putv('arid', arid,
						'orid', orid,
						'mechid', mechid,
						'fm', fm,
						)
			# if there are entries for this arrival already, write over it...
			dbparr.record = dbparr.find('arid=={0} && orid=={1}'.format(arid, orid))
			if dbparr.record < 0:
				dbparr.record = dbparr.addnull()
			dbparr.putv('arid', arid,
						'orid', orid, 
						'esaz', av.azimuth, 
						'dip' , av.takeoff_angle,
						)
	except Exception as e:
		raise e
		#WhiteRumpError("Couldn't write to Antelope! Problem: " + e.message)
	finally:
		db.close()

#--- Useful but unused in this program, to be replaced ---------------#
#def get_waveform_from_arid(database, arid, window=4.):
	#'''Return an ObsPy stream containing traces which match a given
	#arrival ID from Datascope database
	
	#Uses Mark's readANTELOPE function to load from wfdisc to Stream
	#'''
	#db,o = open_db_or_string(database)
	## Pull out needed params from 'arrival'
	#dbv = dbprocess(db, ['dbopen arrival', 'dbsubset arid=={0}'.format(arid)])
	#dbv.record = 0
	#time, sta, chan = dbgetv(dbv,'time','sta','chan')
	## Use sta-chan-teim-endtime to pull out waveform
	#t0 = utc(time)-(window/2.)
	#t1 = t0 + window
	#st = readANTELOPE(database, station=sta, channel=chan, 
		#starttime=t0, endtime=t1)
	#db.close()
	#return st
	
#def change_arrival_fm(database, arid, new_fm):
	#'''Change the first motion string of a pick in the arrival table'''
	#db,o = open_db_or_string(database, perm='r+')
	#dbv = dbprocess(db, ['dbopen arrival', 'dbsubset arid=={0}'.format(arid)])
	#dbv.record = 0
	#d = DbrecordPtr(dbv)
	#lddate = d.lddate
	#d.fm = new_fm
	#d.lddate = lddate
	#db.close()
