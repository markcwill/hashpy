##HASHpype

This is a fork of the first motion focal mechanism program HASH by Hardebeck and Shearer. The subroutines (in Fortran 77, which I did not write) now contain f2py compiiler directives. The Makefile has a new target, 'hashpy', which will make a library, 'libhashpy.so'. This is a python module which will import all the subs and common blocks into the python namespace.

The driver program has been ported to python, and the end goal is to replace it entirely with Python classes and methods, which can then be called from a script. Right now the HashPype class in the 'hashpype' module has all the methods needed for a HASH run. see the 'driver.py' script for details.

Right now the only supported input is BRTT Antelope Datascope databases. The initial port of the driver (examples/hashdriver2.py) supports using the FPFIT-style input files. Eventually, this flat text file functionality may be ported over as well.

This is just a draft implementation, lots of non-bakcwards compatibillty changes ahead...

#Dependencies
* ObsPy (NumPy, SciPy, matplotlib, etc.,)
* mplstereonet
* antelope (if using BRTT Datascope databases)

#Original HASH references

* Hardebeck, Jeanne L. and Peter M. Shearer, A new method for determining first-
  motion focal mechanisms, Bulletin of the Seismological Society of America, 92,
  2264-2276, 2002.
* Hardebeck, Jeanne L. and Peter M. Shearer, Using S/P Amplitude Ratios to
  Constrain the Focal Mechanisms of Small Earthquakes, Bulletin of the
  Seismological Society of America, 93, 2434-2444, 2003.

