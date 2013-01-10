#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# utils.py
#
# by Mark Williams 2012.313
#
# Utilities for using Antelope db with focal mechanisms and obspy
#


import sys, os
from numpy import array
from obspy.core import read, Stream, UTCDateTime
utc = UTCDateTime

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

def db2object(dbv):
	"""
	Port of Antelope MATLAB toolbox 'db2struct' function.
		
	Returns a list-like object, this is the function version of calling
	DbrecordList() directly.
	
	:type dbv: antelope.datascope.Dbptr
	:param dbv: Open pointer to an Antelope database view or table
	:rtype: :class:`~obspy.antelope.Dbview`
	:return: Dbview of Dbrecord objeccts
	"""
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
		dbr = Dbrecord(db)
		t0 = UTCDateTime(dbr.time)
		t1 = UTCDateTime(dbr.endtime)
		if dbr.time < ts:
			t0 = starttime
		if dbr.endtime > te:
			t1 = endtime
		_st = read(fname, starttime=t0, endtime=t1)		 # add format?
		_st = _st.select(station=dbr.sta, channel=dbr.chan) #not location aware
		_st[0].db = dbr
		st += _st
	# Close what we opened, BUT garbage collection may take care of this:
	# if you have an open pointer but pass db name as a string, global
	# use of your pointer won't work if this is uncommented:
	#
	#if isinstance(database,str):
	#	db.close()
	return st

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
