HASHpy
------

This is a fork of HASH v1.2, the first motion focal mechanism program by Hardebeck and Shearer. The subroutines (in Fortran 77, which I did not write) are compiled into a python module, 'libhashpy.so', which will import all the subs and common blocks into the python namespace. There is a base class, HashPype, that contains attributes which hold data for a HASH calculation, and methods which can be called to do the HASH calculation. This class facilitates easily writing a 'hash driver' script in python. See below for details.

### Installation
The latest version of this code uses `numpy.distutils`, which can nicely compile Fortran code. The source code and makefiles to remake the libhashpy module are installed to the hashpy folder, so with a little hacking, one could just remake them in place, if one wanted to change the source code.

For basic install:

```shell
# First, install dependencies. Then cd to this top directory and:
python setup.py install
```

### Testing and Usage
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

Right now, there are no input/output methods. There is one VERY simple output which will print out the event ID with the best strike/dip/rake, but that is it. The idea is to use this as a metaclass, and build classes which inherit from HashPype and define their own I/O. See the EventHashPype class which uses `obspy.core.event` objects for I/O.

```python
# Usage example:
class MyHashRun(HashPype):
    def input(self, input_files):
        # read in from files, database, STDIN, etc

    def output(self):
        # output to a file, command line, etc...

# Then script:
hro = MyHashRun()
hro.input(my_input)
hro.load_velocity_models()
hro.generate_trial_data()
hro.calculate_takeoff_angles()                
check1 = hro.check_minimum_polarity()
check2 = hro.check_maximum_gap()
hro.calculate_hash_focalmech()
hro.add_solution_to_dict()
hro.print_solution_line()
```

### Dependencies

* Fortran compiler + library (tested with gfortran + libgfortran3)
* NumPy (main dependancy, for numerical arrays and f2py)
* [ObsPy](https://github.com/obspy/obspy.git) (Only if using EventHashPype class)
* [mplstereonet](https://github.com/joferkington/mplstereonet.git) (only for Plotter)

### HASH references

* Hardebeck, Jeanne L. and Peter M. Shearer, A new method for determining first-
  motion focal mechanisms, Bulletin of the Seismological Society of America, 92,
  2264-2276, 2002.
* Hardebeck, Jeanne L. and Peter M. Shearer, Using S/P Amplitude Ratios to
  Constrain the Focal Mechanisms of Small Earthquakes, Bulletin of the
  Seismological Society of America, 93, 2434-2444, 2003.

### Future
If I ever get ambitious, I would redo the structure of some of the Fortran subroutines (already started this with vel-subs2) and rewrite them in Fortran 90/95/2003. Probably couldn't ditch the common blocks without a full rewrite, but at least try and get rid of the GOTOs. Biggest improvement would be allocatable arrays, so you could avoid the includes and max-size arrays.

There will probably be small adjustments to the locations and structure of what functions are in what files, but the HashPype class and methods will be the main way to interact with HASH.

Currently working on generic FocalMech and Plotter classes to abstract out the specific algorithm. These will likely be moved to their own module.
