#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  hash_utils.py


def parameter(**kwargs):
	return [kwargs[key] for key in kwargs]

def fortran_include(fname):
	vars = []
	f_inc = open(fname)
	for line in f_inc:
		if 'c' in line[0]:
			pass
		elif 'parameter' in line:
			vars.extend(eval(line))
		else:
			pass
	f_inc.close()
	return vars


