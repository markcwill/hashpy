##HASHpy

This is a fork of HASH v1.2, the first motion focal mechanism program by Hardebeck and Shearer. The subroutines (in Fortran 77, which I did not write) now contain f2py compiler directives. The Makefile has a new target, 'hashpy', which will make a library, 'libhashpy.so'. This is a python module which will import all the subs and common blocks into the python namespace.

### Installation
The current version of this code uses the standard python distutils, NOT the numpy version, which can compile Fortran code. So, for now, install using the Makefile, which first builds the Fortran library, then runs the classic 'python setup.py install'. The source code and makefiles to remake the libhashpy module are installed to the hashpy folder, so with a little hacking, one could just remake them in place, if one wanted to change the source code.

For basic install:
```
make install
```
### Testing and Usage
```python
# This will give you access to all the Fortran compiled subs from python
from hashpy.libhashpy import *
# To use the HashPype class example:
from hashpy import HashPype
```
Right now, there are no input/output methods. There is one VERY simple output which will print out the event ID with the best strike/dip/rake, but that is it. The idea is to use this as a metaclass, and build classes which inherit from HashPype and define their own I/O. See the dbhash project, which uses the hashpy module to build a DbHashPype class that uses Antelope databases for I/O.

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


This is just a draft implementation, lots of non-bakcwards compatibillty changes ahead...

###Dependencies
* Fortran compiler (tested with gfortran)
* NumPy (main dependancy, for numerical arrays and f2py)
* [ObsPy](https://github.com/obspy/obspy.git) (only if plotting and for additional functionality) 
* [mplstereonet](https://github.com/joferkington/mplstereonet.git) (only for plotting)

###HASH references

* Hardebeck, Jeanne L. and Peter M. Shearer, A new method for determining first-
  motion focal mechanisms, Bulletin of the Seismological Society of America, 92,
  2264-2276, 2002.
* Hardebeck, Jeanne L. and Peter M. Shearer, Using S/P Amplitude Ratios to
  Constrain the Focal Mechanisms of Small Earthquakes, Bulletin of the
  Seismological Society of America, 93, 2434-2444, 2003.

###Future
If I ever get ambitious, I would redo the structure of some of the Fortran subroutines (already started this with vel-subs2) and rewrite them in Fortran 90/95/2003. Probably couldn't ditch the common blocks without a full rewrite, but at least try and get rid of the GOTOs. Biggest improvement would be allocatable arrays, so you could avoid the includes and max-size arrays.

There will probably be small adjustments to the locations and structure of what functions are in what files, but the HashPype class and methods will be the main way to interact with HASH.


