#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  hashpype.py
from focalmech import FocalMech

def main():
	
	#------------------------------------------
	# for testing:
	
	year, jday, test_orid = 2012, 206, 933757 # 2012 206 (shot SPE3)
	year, jday, test_orid = 2009, 186, 962482 # 2009 186 7/5
	
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
	foo.get_phases_from_db('/data/{y}/{d}/reno'.format(y=year, d=jday), orid=test_orid)
	foo.load_velocity_models(test_table_list)
	foo.generate_trial_data()
	foo.calculate_takeoff_angles()
	foo.calculate_hash_focalmech()
	foo.plot_beachball()

if __name__ == '__main__':
	main()

