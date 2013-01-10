#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  dbhash.py
#
# fucntions which run HASH using Antelope and hashpy

from hashpy.db.utils import add_antelope_path
from hashpy.db.dbhashpype import DbHashPype
from hashpy.db.plotter2 import PlotterI
from hashpy.focalmech import FocalMech
add_antelope_path()
from antelope.datascope import dbopen


def dbhash_cli(args):
	'''Perform a HASH run using Database input and command line args'''
	
	# Make a blank HASH run object
	hro = DbHashPype()
	
	# Load data from a pf file
	if args.pf:
		hro.load_pf(pffile=args.pf)
	else:
		hro.load_pf()
	
	# Go into database
	hro.get_phases_from_db(args.dbin, evid=args.evid, orid=args.orid)
	
	# Generate preliminary data for run
	hro.load_velocity_models()
	hro.generate_trial_data()
	hro.calculate_takeoff_angles()
	
	check1 = hro.check_minimum_polarity()
	check2 = hro.check_maximum_gap()
	
	# If it passes checks, run HASH
	if check1:
		hro.calculate_hash_focalmech()
	else:
		raise ValueError("Didn't pass check: Min # picks={0} Max gap={1}".format(check1,check2))
	hro.add_solution_to_dict()
	
	# Send to any optional outputs (plot or db 'fplane' table)
	if args.review:
		fmech = FocalMech()
		fmech.load_hash(hro)
		p = PlotterI(fmech)
	elif args.graph:
		hro.plot_stereonet(labels=True)
	else:
		# quick orid/strike/dip/rake line
		hro.print_solution_line()
	
	if args.dbout:
		hro.save_result_to_db(dbout=args.dbout)
	
	# For interactive scripting and debugging:
	return hro


def dbhash_loc2(args):
	'''Perform a HASH run using Database input from dbloc2 menu'''
	
	# Make a blank HASH run object
	hro = DbHashPype()
	
	# Load data from a pf file
	hro.load_pf()
	
	### Parse special command line args from dbloc2 'origin' menu:
	# Go into database
	dbin = args.dbin.rstrip('.origin')
	db = dbopen(dbin).lookup(table='origin')
	# hack, if in loc2 mode, dbout is actually rec # of 'tmp/trial'
	db.record = int(args.dbout)
	orid = db.getv('orid')[0]
		
	hro.get_phases_from_db(dbin, orid=orid)
	
	# Generate preliminary data for run
	hro.load_velocity_models()
	hro.generate_trial_data()
	hro.calculate_takeoff_angles()
	
	check1 = hro.check_minimum_polarity()
	check2 = hro.check_maximum_gap()
	
	# If it passes checks, run HASH
	if check1:
		hro.calculate_hash_focalmech()
	else:
		raise ValueError("Didn't pass check: # picks = {0} | Minimum = {1}".format(hro.npol,hro.npolmin))
	hro.add_solution_to_dict()
	
	fmech = FocalMech()
	fmech.load_hash(hro)
	p = PlotterI(fmech)
	
	# For interactive scripting and debugging:
	return hro

if __name__ == '__main__':
	print "dbhash.py is a module for the dbhash program.\n Import 'dbhash_cli' or 'dbhash_loc2'"

