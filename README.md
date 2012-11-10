##HASHpype

This is a fork of the first motion focal mechanism program HASH by Hardebek and Shearer. The subroutines (in Fortran 77, which I did not write) now contain f2py compiiler directives. The Makefile has a new target, 'hashpy', which will make a library, 'hashpy.so'. This is a python module which will import all the subs and common blocks into the python namespace.

The driver program has been ported to python, and the end goal is to replace it entirely with classes and methods, which can then be called from a script. Right now the FocalMech class hass all the methods needed for a HASH run. see the 'hashpype' script for details.

Right now the only supportly input is BRTT Antelope Datascope databases. The port of the driver supports using the FPFIT-style input files. Eventually, this flat text file functionality will be ported over.

This is just a draft implementation, lots of non-bakcwards compatibillty changes ahead...

* Hardebeck, Jeanne L. and Peter M. Shearer, A new method for determining first-
  motion focal mechanisms, Bulletin of the Seismological Society of America, 92,
  2264-2276, 2002.
* Hardebeck, Jeanne L. and Peter M. Shearer, Using S/P Amplitude Ratios to
  Constrain the Focal Mechanisms of Small Earthquakes, Bulletin of the
  Seismological Society of America, 93, 2434-2444, 2003.

