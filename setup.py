#!/usr/bin/env python

from distutils.core import setup

setup(name='HASHpy',
      version='1.0',
      description='Routines for running HASH algorithms',
      author='Mark Williams',
      url='https//github.com/markcwill',
      packages=['hashpy'],
      package_data={'hashpy':['libhashpy.so','src/*.inc','Makefile','data/*', 'src/*.f']},
)


