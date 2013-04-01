# -*- coding: utf-8 -*-
#
#  dbhash.py
#
# fucntions which run HASH using Antelope and hashpy
from obspy.core.utcdatetime import UTCDateTime
from hashpy.eventhashpype import EventHashPype, HashError
from hashpy.db.utils import add_antelope_path, load_pf, readANTELOPE
from hashpy.db.dbhashpype import DbHashPype
from hashpy.db.database import DbConnection, db2event
from hashpy.db.plotter2 import PlotterI
from hashpy.db.plotter3 import FocalMechPlotter
from hashpy.doublecouple import FocalMech
add_antelope_path()
from antelope.datascope import dbopen

# FOR TESTING ---------------------------------------------------------#
class Arg():
	pass

args = Arg()
args.dbin = 'reno'
args.orid = 979567
args.pf = None

def dbhash_2(args):
	#--- Little 'dbhash' script ---#
	# Load settings data from a pf file
	if args.pf:
		pf = load_pf(pffile=args.pf)
	else:
		pf = load_pf()
	
	ev = db2event(args.dbin, orid=args.orid, phase_data=True)
	
	hro = EventHashPype(**pf)
	hro.load_event(ev)
	
	try:
		hro.run(gap_check=False)
	except HashError as e:
		print "Failed! " + e.message
	except:
		raise
	
	if hro.nmult:
		ev = hro.output(event=ev, only_fm_picks=True)
		# little script to get the FM pick waveform data for plotting
		dbc = DbConnection(args.dbin)
		adb = dbc.get_first_motions(orid=args.orid)
		t0 = UTCDateTime(adb[0]['arrival.time']) - 10
		t1 = UTCDateTime(adb[-1]['arrival.time']) + 10
		st = readANTELOPE(adb.Ptr, starttime=t0, endtime=t1)
		dbc.close()
		p = FocalMechPlotter(ev, datastream=st,) #database=args.dbin)
	return hro, ev
#----------------------------------------------------------------------#
def dbhash_cli(args):
	"""
	Perform a HASH run using Database input and command line args
	"""
	# Load settings data from a pf file...
	if args.pf:
		pf = load_pf(pffile=args.pf)
	else:
		pf = load_pf()
	# Grab data from the db...
	ev = db2event(args.dbin, orid=args.orid, phase_data=True)
	# Create a HASH run instance and load the data...
	hro = EventHashPype(**pf)
	hro.load_event(ev)
	# Run and catch errors from the minimum requirements checks
	try:
		hro.run(gap_check=False)
	except HashError as e:
		print "Failed! " + e.message
	except:
		raise
	# Grab waveform data and launch plotter or spit out solution...
	if args.review or args.graph:
		if hro.nmult:
			ev = hro.output(event=ev, only_fm_picks=True)
			# little script to get the FM pick waveform data for plotting
			dbc = DbConnection(args.dbin)
			adb = dbc.get_first_motions(orid=args.orid)
			t0 = UTCDateTime(adb[0]['arrival.time']) - 10
			t1 = UTCDateTime(adb[-1]['arrival.time']) + 10
			st = readANTELOPE(adb.Ptr, starttime=t0, endtime=t1)
			dbc.close()
			p = FocalMechPlotter(ev, datastream=st,) #database=args.dbin)
	else:
		# quick orid/strike/dip/rake line
		print hro

	if args.dbout:
		raise NotImplementedError("No support for event2db yet!!!")
		#hro.save_result_to_db(dbout=args.dbout)
	
	# For interactive scripting and debugging:
	return hro


def dbhash_loc2(args):
	'''Perform a HASH run using Database input from dbloc2 menu'''
	### Parse special command line args from dbloc2 'origin' menu:
	# Go into database
	dbin = args.dbin.rstrip('.origin')
	db = dbopen(dbin).lookup(table='origin')
	# hack, if in loc2 mode, dbout is actually rec # of 'tmp/trial'
	db.record = int(args.dbout)
	orid = db.getv('orid')[0]
	
	# For interactive scripting and debugging:
	ev1 = db2event(dbin, orid=orid, phase_data=True)
	# Load data from a pf file
	if args.pf:
		pf = load_pf(pffile=args.pf)
	else:
		pf = load_pf()
	hro = EventHashPype(**pf)
	hro.load_event(ev1)
	hro.run()
	ev2 = hro.output()
	p = FocalMechPlotter(ev2, database=args.dbin)
	return hro

if False: #__name__ == '__main__':
	from argparse import ArgumentParser
	parser = ArgumentParser()
	parser.add_argument("dbin",   help="Input database")
	parser.add_argument("dbout",  help="Output database", nargs='?')
	parser.add_argument("-g", "--graph", help="Plot result", action='store_true')
	parser.add_argument("-r", "--review", help="Interactive reviewer mode", action='store_true')
	parser.add_argument("-l", "--loc", help="dbloc2 mode", action='store_true')
	parser.add_argument("--pf",   help="Parameter file")
	group = parser.add_mutually_exclusive_group() #required=True)
	group.add_argument("--evid", help="Event ID", type=int)
	group.add_argument("--orid", help="Origin ID", type=int)
	args = parser.parse_args()
	if args.loc:
		hash_prog = dbhash_loc2
	else:
		hash_prog = dbhash_cli
	hro = hash_prog(args)
	# done
