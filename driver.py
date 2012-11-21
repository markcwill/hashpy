#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  hashpype.py

from hashpype import HashPype 

def main():
	
	#------------------------------------------
	# for testing:
	
	year, jday, test_orid = 2012, 206, 933757 # 2012 206 (shot SPE3)
	#year, jday, test_orid = 2009, 186, 962482 # 2009 186 7/5
	
	c_evid = '''
	379439 2012167
	380226 2012174
	380653 2012175
	380828 2012175
	'''
	
	ttl =  ['vz.socal',
						'vz.north',
						'vz.lab1',
						'vz.sgm1',
						'vz.vb1']
	
	vel_mod_dir = 'data'
	ttl= ['vz.pickema1', 'vz.pickema2', 'vz.pickema3']
	test_table_list = ['/'.join([vel_mod_dir,table]) for table in ttl]
	
	hro = HashPype(maxout=100)
	hro.get_phases_from_db('/data/{y}/{d}/reno'.format(y=year, d=jday), orid=test_orid)
	hro.load_velocity_models(test_table_list)
	hro.generate_trial_data()
	hro.calculate_takeoff_angles()
	hro.calculate_hash_focalmech()
	#hro.print_solution_line()
	hro.plot_beachball(labels=True)
	return hro

if __name__ == '__main__':
	hro = main()

