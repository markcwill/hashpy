# -*- coding: utf-8 -*-
#
#
"""
obspy antelope module utilities
"""
import sys, os
def add_antelope_path():
	_version_string = os.environ['ANTELOPE'].split('/')[-1]
	_pydirs = ['data','python']
	if float(_version_string[:3]) < 5.2:
		_pydirs = ['local'] + _pydirs
	_pypath = os.path.join(os.environ['ANTELOPE'], *_pydirs)
	if _pypath not in sys.path:
		sys.path.append(_pypath)


add_antelope_path()
from antelope.datascope import Dbptr, dbopen

def open_db_or_string(database, perm='r'):
	'''
	Check if a variable is a valid db or a string

	Returns a pointer to an open db or throw an error
	'''
	if isinstance(database, Dbptr):
		ptr = Dbptr(database)
	elif isinstance(database, str):
		ptr = dbopen(database, perm)
		opened = True
	else:
		raise TypeError("Input must be a Dbptr or string of a valid database path")
	return ptr, opened


