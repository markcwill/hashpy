##HASHpype

This is a fork of the first motion focal mechanism program HASH by Hardebeck and Shearer. The subroutines (in Fortran 77, which I did not write) now contain f2py compiler directives. The Makefile has a new target, 'hashpy', which will make a library, 'libhashpy.so'. This is a python module which will import all the subs and common blocks into the python namespace.

The driver program has been ported to python, and the end goal is to replace it entirely with Python classes and methods, which can then be called from a script. Right now the HashPype class in the 'hashpype' module has all the methods needed for a HASH run. One can then create a pipe-like workflow using this class which holds all the necessary data. See the 'driver.py' script for details.

Right now the only supported input is BRTT Antelope Datascope databases. The initial port of the driver (hashdriver2.py) supports using the FPFIT-style input files. Eventually, this flat text file functionality may be ported over to the pipe class as well.

This is just a draft implementation, lots of non-bakcwards compatibillty changes ahead...

###Dependencies
* NumPy & SciPy (main dependancy, for numerical arrays and f2py)
* ObsPy (only if plotting and for additional functionality) 
* mplstereonet (only for plotting)
* antelope (if using BRTT Datascope databases)

###HASH references

* Hardebeck, Jeanne L. and Peter M. Shearer, A new method for determining first-
  motion focal mechanisms, Bulletin of the Seismological Society of America, 92,
  2264-2276, 2002.
* Hardebeck, Jeanne L. and Peter M. Shearer, Using S/P Amplitude Ratios to
  Constrain the Focal Mechanisms of Small Earthquakes, Bulletin of the
  Seismological Society of America, 93, 2434-2444, 2003.

