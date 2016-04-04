HASHpy
------

[![DOI](https://zenodo.org/badge/3723/markcwill/hashpy.png)](http://dx.doi.org/10.5281/zenodo.9808)

This is a fork of HASH v1.2, the first motion focal mechanism program by Hardebeck and Shearer. The original Fortran subroutines are compiled into a python module, 'libhashpy.so', which will import all the subs and common blocks into the python namespace. There is a base class, HashPype, that contains attributes which hold data for a HASH calculation, and methods which can be called to do the HASH calculation. This class facilitates easily writing a 'hash driver' script in python. See below for details.

The code base has been forked at HASH 1.2, and re-syntaxed from Fortran 77 to Fortran 95. More work is needed to take advantage of F95 features (allocatable arrays, etc) but many GOTOs have been removed from the base subroutines. A more complete refactor would involve changing the Python wrapper code as well, and will be reserved for a '2.0' hashpy version.


### Installation

This package setup uses `numpy.distutils`, which can nicely compile Fortran code. The source code and makefiles to remake the libhashpy module are installed to the hashpy folder, so with a little hacking, one could just remake them in place, if one wanted to change the source code.

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
Type:       module
String Form:<module 'hashpy.libhashpy' from '/home/seismo/homes/markwill/Code/python/hashpy/hashpy/libhashpy.so'>
File:       /home/seismo/homes/markwill/Code/python/hashpy/hashpy/libhashpy.so
Docstring:
This module 'libhashpy' is auto-generated with f2py (version:2).
Functions:
  nf,strike,dip,rake,faults,slips = focalamp_mc(p_azi_mc,p_the_mc,sp_amp,p_pol,nmc,dang,maxout,nextra,ntotal,qextra,qtotal,npsta=len(sp_amp))
  mfrac,mavg,stdr = get_misf_amp(p_azi_mc,p_the_mc,sp_ratio,p_pol,str_avg,dip_avg,rak_avg,npol=len(p_azi_mc))
  nf,strike,dip,rake,faults,slips = focalmc(p_azi_mc,p_the_mc,p_pol,p_qual,nmc,dang,maxout,nextra,ntotal,npsta=len(p_pol))
  nsltn,str_avg,dip_avg,rak_avg,prob,rms_diff = mech_prob(norm1in,norm2in,cangle,prob_max,nf=shape(norm1in,1))
  norm1_avg,norm2_avg = mech_avg(norm1,norm2,nf=shape(norm1,1))
  rota = mech_rot(norm1,norm2,slip1,slip2)
  v3 = cross(v1,v2)
  x,y,z = to_car(the,phi,r)
  strike,dip,rake,fnorm,slip = fpcoor(strike,dip,rake,fnorm,slip,idir)
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
  mk_table_add(itab,vmodel)
COMMON blocks:
  /angtable/ table(101,14,10),delttab(101),deptab(14),ndel,ndep
.

```

### Input/Output

Data can be input into HASH in various formats. There are some I/O routine submodules in the `hashpy.io` module. I/O functions can be imported from here and used directly, or set to the class attributes `input_factory` and `output_factory`. The functions will then called by the HashPype methods `HashPype.input()` and `HashPype.output()`.

Current IO submodules:
* `hashpy.io.obspyIO` - ObsPy Event object I/O (`inputOBSPY` and `outputOBSPY`)
* `hashpy.io.antelopeIO` - Antelope Datascope database and other I/O (`input`, `output`, plus other fxns)
* `hashpy.io.fpfitIO` - Old style FPFIT file I/O (still in develoment)

The default I/O functions input keyword arguments as HashPype attribute and output a string of the "best" solution.


```python
# Usage example:
# Typical "hash_driver2" style script
# Using the ObsPy Event format as input for origin, picks, and arrivals

from hashpy import HashPype, HashError
from hashpy.io.obspyIO import inputOBSPY, outputOBSPY

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

hp = HashPype(input_factory=inputOBSPY, **config)
hp.input(event)
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

### Plotting (Experimental)

A trial implementation of plotting exists, using `matplotlib` and the `mplstereonet` package, as the  `hashpy.plotting.focalmechplotter.FocalMechPlotter` class. It accepts an ObsPy Event containing Picks, Origin/Arrivals, FocalMechanism, etc, objects (as output from HashPype) and generates a stereonet plot. Multiple FocalMechansim solutions from HASH are accessible through the navigation toolbar 'back' and 'forward' arrows.

```python
# Get an obspy Event object as output
>>> event = outputOBSPY(hp)
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
* [curds2](http://github.com/NVSeismoLab/curds2.git) (Only for Antelope database I/O)
* [Antelope](http://www.brtt.com) (Required by curds2)

### HASH references

* Hardebeck, Jeanne L. and Peter M. Shearer, A new method for determining first-
  motion focal mechanisms, Bulletin of the Seismological Society of America, 92,
  2264-2276, 2002.
* Hardebeck, Jeanne L. and Peter M. Shearer, Using S/P Amplitude Ratios to
  Constrain the Focal Mechanisms of Small Earthquakes, Bulletin of the
  Seismological Society of America, 93, 2434-2444, 2003.

### Future
Redo the structure of some of the Fortran subroutines (already started this with vel-subs2) to use Fortran 90/95/2003 new features. Probably couldn't ditch the common blocks without a full rewrite, but at least try and get rid of the GOTOs. Biggest improvement would be allocatable arrays, so you could avoid the includes and max-size arrays.

There will probably be small adjustments to the locations and structure of what functions are in what files, but the HashPype class and methods will be the main way to interact with HASH.

