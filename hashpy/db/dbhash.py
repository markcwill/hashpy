#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  dbhash.py
#  

from hashpy.db.dbhashpype import DbHashPype
from hashpy.db.plotter import Plotter
from hashpy.focalmech import FocalMech
from argparse import ArgumentParser

def dbhash(args):
	'''Perform a HASH run using Database input'''
	
	# Make a blank HASH run object
	hro = DbHashPype()
	
	# Load data from a pf file
	if args.pf:
		hro.load_pf(pffile=args.pf)
	
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
		p = Plotter(fmech)
	elif args.graph:
		hro.plot_stereonet(labels=True)
	else:
		# quick orid/strike/dip/rake line
		hro.print_solution_line()
	
	if args.dbout:
		hro.save_result_to_db(dbout=args.dbout)
	
	# For interactive scripting and debugging:
	return hro

# RUN AS A SCRIPT
if __name__ == '__main__':
	parser = ArgumentParser()
	parser.add_argument("dbin",   help="Input database")
	parser.add_argument("dbout",  help="Output database", nargs='?')
	parser.add_argument("-g", "--graph", help="Plot result", action='store_true')
	parser.add_argument("-r", "--review", help="Interactive reviewer mode", action='store_true')
	parser.add_argument("--pf",   help="Parameter file")
	group = parser.add_mutually_exclusive_group(required=True)
	group.add_argument("--evid", help="Event ID", type=int)
	group.add_argument("--orid", help="Origin ID", type=int)
	args = parser.parse_args()
	hro = dbhash(args)

