#!/usr/bin/env python

from distutils.core import setup

setup(name='AUG contrib',
      version='1.0',
      description='Antelope Users Contributed Modules',
      author='Mark Williams',
      url='https//github.com/NVSeismoLab',
      packages=['aug','aug.contrib',
                'aug.contrib.orm'],
)


