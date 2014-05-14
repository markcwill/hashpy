HASHpy
------

[![DOI](https://zenodo.org/badge/3723/markcwill/hashpy.png)](http://dx.doi.org/10.5281/zenodo.9808)

This is a fork of HASH v1.2, the first motion focal mechanism program by Hardebeck and Shearer. The subroutines (in Fortran 77, which I did not write) are compiled into a python module, 'libhashpy.so', which will import all the subs and common blocks into the python namespace. There is a base class, HashPype, that contains attributes which hold data for a HASH calculation, and methods which can be called to do the HASH calculation. This class facilitates easily writing a 'hash driver' script in python. See below for details.

Note: As stated in the in-code docs,  the current code is based on the 'hashdriver2' script, and as such, does not utilize the 'amp' routines for S/P amplitude measurements. This will most likely be added in the near future. If someone wants to take it upon themselves to add the compiler directives in the source and the methods to the HashPype class, throw me a pull request.

### Installation

The latest version of this code uses `numpy.distutils`, which can nicely compile Fortran code. The source code and makefiles to remake the libhashpy module are installed to the hashpy folder, so with a little hacking, one could just remake them in place, if one wanted to change the source code.

For basic install:

```shell
# First, install dependencies. Then cd to this top directory and:
python setup.py install

# or use pip
pip install git+https://github.com/markcwill/hashpy.git

```

### Testing and Usage

For lower-level programming, the Fortran subroutines are available in `hashpy.libhashpy`

```python
In [1]: from hashpy import libhashpy

In [2]: libhashpy?
Type:       module
Base Class: <type 'module'>
String Form:<module 'hashpy.libhashpy' from '/opt/antelope/python2.7.2-64/lib/python2.7/site-packages/hashpy/libhashpy.so'>
Namespace:  Interactive
File:       /opt/antelope/python2.7.2-64/lib/python2.7/site-packages/hashpy/libhashpy.so
Docstring:
This module 'libhashpy' is auto-generated with f2py (version:1).
Functions:
  nf,strike,dip,rake,faults,slips = focalmc(p_azi_mc,p_the_mc,p_pol,p_qual,nmc,dang,maxout,nextra,ntotal,npsta=len(p_pol))
  nsltn,str_avg,dip_avg,rak_avg,prob,rms_diff = mech_prob(norm1in,norm2in,cangle,prob_max,nf=shape(norm1in,1))
  norm1_avg,norm2_avg = mech_avg(norm1,norm2,nf=shape(norm1,1))
  rota = mech_rot(norm1,norm2,slip1,slip2)
  v3 = cross(v1,v2)
  to_car(the,phi,r,x,y,z)
  fpcoor(strike,dip,rake,fnorm,slip,idir)
  fran = ran_norm()
  mfrac,stdr = get_misf(p_azi_mc,p_the_mc,p_pol,p_qual,str_avg,dip_avg,rak_avg,npol=len(p_azi_mc))
  magap,mpgap = get_gap(p_azi_mc,p_the_mc,npol=len(p_azi_mc))
  sort(ra,n=len(ra))
  ntab = mk_table(ntab)
  tt,iflag = get_tts(ip,del,qdep)
  layertrace(p1,h1,utop1,ubot1,imth,dx1,dt1,irtr)
  flat,flon,felev = getstat_tri(stlfile,snam,scom,snet)
  stpol = check_pol(polfile,snam,evyr,evmon,evdy,evhr)
  qcor = get_cor(stlfile,snam,scom,snet)
  ntab = mk_table_add(ind,vmodel)
COMMON blocks:
  /angtable/ table(101,14,10),delttab(101),deptab(14),ndel,ndep
.

In [3]: 
```

### Input/Output

Data can be input into HASH in various formats. This is currently handled in HASHpy by registering a format in the `hashpy.io` module by adding a module/functions to the dictionary in `hashpy.io.core`. (May change in future releases). I/O functions are then called by the HashPype methods `HashPype.input()` and `HashPype.output()`. If the 'output' method is called with no format it will return a simple string with the event ID with the best strike/dip/rake

Currently Supports:
* ObsPy Event object I/O
* Antelope Datascope database I/O

```python
# Usage example:
# Typical "hash_driver2" style script
# Using the ObsPy Event format as input for origin, picks, and arrivals

from hashpy import HashPype, HashError

# Make an ObsPy Event, or get from a QuakeML file
from obspy.core.event import readEvents
event = readEvents('my_quakeml_file.xml').events[0]

# Set configuration at creation with a dict...
# ...can from file or interactively, etc
config = { "npolmin" : 10,
           "max_agap": 90,
           "vmodels" : ['/path/to/my/vmodel/file1.vz', 
                        '/new/picking/model/file2.vz',
                       ] 
           }

hp = HashPype(**config)
hp.input(event, format="OBSPY")
hp.load_velocity_models()
hp.generate_trial_data()
hp.calculate_takeoff_angles()

pass1 = hp.check_minimum_polarity()
pass2 = hp.check_maximum_gap()

if pass1 and pass2:
    hp.calculate_hash_focalmech()
    hp.calculate_quality()
    print hp.output() # default output is a simple string
else:
    raise HashError("Didn't pass user checks!")

```

### Plotting

A trial implementation of plotting exists, using `matplotlib` and the `mplstereonet` package, as the  `hashpy.plotting.focalmechplotter.FocalMechPlotter` class. It accepts an ObsPy Event containing Picks, Origin/Arrivals, FocalMechanism, etc, objects (as output from HashPype) and generates a stereonet plot. Multiple FocalMechansim solutions from HASH are accessible through the navigation toolbar 'back' and 'forward' arrows.

```python
# Get an obspy Event object as output
>>> event = hp.output(format="OBSPY")
# Pass to plotter class as constructor variable
>>> fmp = FocalMechPlotter(event)
# Plots a figure, accessible as 'fmp.fig'
```
![](http://markcwill.github.io/hashpy/images/979567_focalmech.png)

### Dependencies

#### Required
* Fortran compiler + library (tested with gfortran + libgfortran3)
* NumPy (main dependancy, for numerical arrays and f2py compiling)

#### Optional
* [ObsPy](https://github.com/obspy/obspy.git) (Only for plotting and ObsPy I/O))
* [mplstereonet](https://github.com/joferkington/mplstereonet.git) (Only for plotting)
* [Antelope](http://www.brtt.com) (Only for Antelope database I/O)

### HASH references

* Hardebeck, Jeanne L. and Peter M. Shearer, A new method for determining first-
  motion focal mechanisms, Bulletin of the Seismological Society of America, 92,
  2264-2276, 2002.
* Hardebeck, Jeanne L. and Peter M. Shearer, Using S/P Amplitude Ratios to
  Constrain the Focal Mechanisms of Small Earthquakes, Bulletin of the
  Seismological Society of America, 93, 2434-2444, 2003.

### Future
Add remaining HASH routines (like those in fmamp_subs.f) to wrapped library and class methods.

If I ever get ambitious, I would redo the structure of some of the Fortran subroutines (already started this with vel-subs2) and rewrite them in Fortran 90/95/2003. Probably couldn't ditch the common blocks without a full rewrite, but at least try and get rid of the GOTOs. Biggest improvement would be allocatable arrays, so you could avoid the includes and max-size arrays.

There will probably be small adjustments to the locations and structure of what functions are in what files, but the HashPype class and methods will be the main way to interact with HASH.

