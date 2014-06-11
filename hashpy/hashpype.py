# -*- coding: utf-8 -*-
"""
Module to call and run HASH subroutines in Fortran
 
    Contains: HashPype
    First motion focal mechanism class for running HASH


:copyright:
    Mark Williams (2013) 
    Nevada Seismological Laboratory

"""

# import all Fortran subroutines and common blocks from the mod
# plus my custom utils and numpy arrays

import numpy as np
import os
from pwd import getpwuid
from hashpy.io.core import Inputter, Outputter
from hashpy.libhashpy import (mk_table_add, angtable, ran_norm, get_tts, get_gap,
    focalmc, mech_prob, get_misf, focalamp_mc, get_misf_amp)

def parameter(**kwargs):
    """
    Returns variables inside a fortran 'parameter' call
    """
    return kwargs

def fortran_include(fname):
    """
    Fortran include
    
    Quickie for access to 'include' files, needed for access to
    parameters useful for interacting with f2py wrapped functions
    
    Input
    -----
    fname : str of fortran include filename
    
    Returns : dict of parameters described in any 'parameter' call
    """
    inc_params = {}
    f = open(fname)
    for line in f:
        if 'c' in line[0]:
            continue
        elif 'parameter' in line:
            # evaluate the line 'parameter(foo=999)' as python expr
            inc_params.update(eval(line))
        else:
            pass
    f.close()
    return inc_params
    
    
class HashPype(object):
    """
    Object which holds all data from a HASH instance for one event
    
    Methods are based on the 'hash_driver2.f' program from H & S
    
    The variables are named and structured as close to the original code
    as possible. Methods act on the variables within this namespace and
    most are accessible as attributes.
    
    One can make a "HASH" driver program by calling the methods of this
    class on the data held within it, see the 'driver2' method. Input and
    output formats can be registered in the 'io' package, in the 'core'
    module.
    
    Methods
    -------
      input(data, format=None, *args, **kwargs)
      load_velocity_models(model_list=None)
      generate_trial_data()
      calculate_takeoff_angles()
      view_polarity_data()
      check_minimum_polarity()
      check_maximum_gap()
      calculate_hash_focalmech()
      calculate_quality()
      output(format=None, *args, **kwargs)

      driver2(check_for_minimum_picks=True, check_for_maximum_gap_size=True)

    """
    # These MUST be the same as the fortran includes!!
    # (They are compiled into the Fortran subroutines)
    #--- param.inc ---#
    npick0, nmc0, nmax0 = None, None, None
    #--- rot.inc ---#
    dang0, ncoor        = None, None
    
    # initialize arrays
    
    # Input arrays
    sname     = None
    scomp     = None
    snet      = None
    pickpol   = None
    pickonset = None
    p_pol     = None
    p_qual    = None
    p_azi_mc  = None
    p_the_mc  = None
    index     = None
    qdep2     = None
    sp_ratio  = None
    
    # Output arrays
    f1norm  = None
    f2norm  = None
    strike2 = None
    dip2    = None
    rake2   = None
    str_avg = None
    dip_avg = None
    rak_avg = None
    var_est = None
    var_avg = None
    mfrac   = None
    mavg    = None
    stdr    = None
    prob    = None
    qual    = None
    
    # Chosen fm run setings: (a config file could look like this:)
    # dbhash.conf -----------------------------------------------------
    # Defaults
    npolmin  = 8    # Enter mininum number of polarities (e.g., 8)
    max_agap = 90   # Enter maximum azimuthal gap (e.g., 90)
    max_pgap = 60   # Enter maximum takeoff angle gap (e.g., 60)
    dang     = 5    # Enter grid angle for focal mech search, in degrees 5(min {0})
    nmc      = 30   # Enter number of trials (e.g., 30)
    maxout   = 500  # Enter maxout for focal mech. output (e.g., 500)
    badfrac  = 0.1  # Enter fraction of picks presumed bad (e.g., 0.10)
    delmax   = 500  # Enter maximum allowed source-station distance, in km (e.g., 120)
    cangle   = 45   # Enter angle for computing mechanisms probability, in degrees (e.g., 45)
    prob_max = 0.1  # Enter probability threshold for multiples (e.g., 0.1)
    ratmin   = 3.0  # Enter minimum allowed signal-to-noise ratio
    qbadfrac = 0.3  # Enter the assumed noise in amplitude ratios, log10 (e.g. 0.3 for factor of 2)
    #----------------------------------------------------------------
    
    # Mark's spec'd object variables
    vmodel_dir = None # string of directory containing velocity models
    vmodels = []      # list of string names of velocity model files
    vtables = None    # actual travel time tables in common block
    arid    = None    # unique number for each pick...
    author  = None    # logged in user
    
    # Variables HASH keeps internally, for ref and passing to fxns
    ntab = 0     # number of tables loaded
    npol = 0     # number of observations
    nppl = 0     # number of p polarities
    nspr = 0     # number of s/p ratios
    dist = None  # distance from source im km
    qazi = None  # azimuth from event to station
    flat = None  # pick station lat
    flon = None  # pick station lon
    felv = None  # pick station elv
    nf2 = None   # number of fm's returned
    nmult = None # number of fm solutions returned
    magap = None # calulated max azi gap
    mpgap = None # calculated max plunge gap (I think)
    
    tstamp  = None
    qlat    = None
    qlon    = None
    qdep    = None
    qmag    = None
    icusp   = None
    seh     = None
    sez     = None
    rms     = None
    terr    = None
    
    # polarity reversals, [-1,1] stub for now
    spol = 1
    
    def __init__(self, **kwargs):
        """
        Make an empty HASH run object.
        
        Initialize all the arrays needed for a HASH run. They are based
        on the maximum size of the arrays passed to the FORTRAN
        subroutines. 
        
        ANYTHING initialized here can be changed by passing
        a keyword arg at the Constructor call. Useful for using an
        input parameter file for HASH.
        
        Example
        -------
        >>> param_dict_from_some_input = { 'npolmin': 20, 'badfrac': 0.25 }
        >>> h = HashPype(**param_dict_from_some_input)
        """
        
        # INCLUDES ------------------------------------------------
        # contain array sizes, get em
        directory = os.path.dirname(__file__)
        param_inc_file = os.path.join(directory,'src','param.inc')
        rot_inc_file   = os.path.join(directory,'src','rot.inc')
        
        param_inc = fortran_include(param_inc_file) # npick0, nmc0, nmax0 
        rot_inc = fortran_include(rot_inc_file) # dang0, ncoor
        
        # Save include vars for other fucntions to access
        self.__dict__.update(param_inc)
        self.__dict__.update(rot_inc)
        
        # For convenience, throw in fxn namespace
        npick0 = self.npick0
        nmc0   = self.nmc0
        nmax0  = self.nmax0
        dang0  = self.dang0
        ncoor  = self.ncoor
        
        # ARRAYS ---------------------------------------------------
        # initialize arrays for HASH
        self.dang2=max(self.dang0, self.dang)
        
        # Input arrays
        self.sname     = np.empty(npick0, 'a6', 'F')
        self.scomp     = np.empty(npick0, 'a3', 'F')
        self.snet      = np.empty(npick0, 'a2', 'F')
        self.pickpol   = np.empty(npick0, 'a1', 'F')
        self.pickonset = np.empty(npick0, 'a1', 'F')
        self.p_pol     = np.empty(npick0, int, 'F')
        self.p_qual    = np.empty(npick0, int, 'F')
        self.spol      = np.empty(npick0, int, 'F')
        self.sp_ratio  = np.empty(npick0, float, 'F')
        self.p_azi_mc  = np.empty((npick0,nmc0), float, 'F')
        self.p_the_mc  = np.empty((npick0,nmc0), float, 'F')
        self.index     = np.empty(nmc0, int, 'F')
        self.qdep2     = np.empty(nmc0, float, 'F')
        
        # Output arrays
        self.f1norm  = np.empty((3,nmax0), float, 'F')
        self.f2norm  = np.empty((3,nmax0), float, 'F')
        self.strike2 = np.empty(nmax0, float, 'F')
        self.dip2    = np.empty(nmax0, float, 'F')
        self.rake2   = np.empty(nmax0, float, 'F')
        self.str_avg = np.empty(5, float, 'F')
        self.dip_avg = np.empty(5, float, 'F')
        self.rak_avg = np.empty(5, float, 'F')
        self.var_est = np.empty((2,5), float, 'F')
        self.var_avg = np.empty(5, float, 'F')
        self.mfrac   = np.empty(5, float, 'F')
        self.mavg    = np.empty(5, float, 'F')
        self.stdr    = np.empty(5, float, 'F')
        self.prob    = np.empty(5, float, 'F')
        self.qual    = np.empty(5, 'a', 'F')
        
        # Internals (added for Python classiness)
        self.dist    = np.empty(npick0, float)
        self.qazi    = np.empty(npick0, float)
        self.flat    = np.empty(npick0, float)
        self.flon    = np.empty(npick0, float)
        self.felv    = np.empty(npick0, float)
        self.esaz    = np.empty(npick0, float)
        self.arid    = np.empty(npick0, int) * 0
        
        self.author = getpwuid(os.getuid()).pw_name
        
        if kwargs:
            self.__dict__.update(kwargs)
    
    def __repr__(self):
        """
        String saying what you are
        """
        return '{0}(icusp={1})'.format(self.__class__.__name__, self.icusp)
    
    def input(self, data, format=None, *args, **kwargs):
        """
        Input data using a formatting standard from the 'io' module
        
        An input function is passed the current HashPype, the object for
        input (could be object instance, open file handle, string of filename, etc)
        and any other args and kwargs required. New inputs can be added to the
        `hashpy.io` module and registered in `hashpy.io.core`

        :param data: The picks/data to be input into a HashPype run, format
            depends on the input function
        :param str format: A  format type registered in `hashpy.io` module
        
        :param args:    Additional arguments are passed to the input function
        :param kwargs:  Additioanl keywords are passed to the input function
        
        """
        try:
            inputter = Inputter(format=format)
        except ImportError as ierr:
            raise ImportError("Couldn't import module for the format {0}: {1}".format(format, ierr.message))
        inputter(self, data, *args, **kwargs)

    def output(self, format=None, *args, **kwargs):
        """
        Output data using a formatting standard from the 'io' module
        
        The current module will output a simple one-line string of the
        best solution stored in the run.
        
        :params str format: A format type registered in `hashpy.io` module

        """
        try:
            outputter = Outputter(format=format)
        except:
            raise IOError("Can't find a module for the format {0}".format(format))
        return outputter(self, *args, **kwargs)
        
    def load_velocity_models(self, model_list=None):
        """
        Load velocity model data
        
        If a None is specified (the default), will use the list in 
        `HashPype.vmodels`
        
        :type model_list:  list of str
        :param model_list: Optional velocity model filenames
        
        """
        # Future -- allow adding on fly, check and append to existing
        #
        # THIS LOADS INTO hashpy.angtable.table (nx,nd,nindex) !!!
        
        # take care of 
        if model_list:
            models = model_list
        else:
            models = self.vmodels
        
        for n,v in enumerate(models):
            self.ntab = mk_table_add(n+1,v)
            self.vtable = angtable.table
        
    def generate_trial_data(self):
        """
        Make data for running trials
        
        (MUST have loaded data and vel mods already)
        
        Algorithm by H & S (From HASH driver script)
        """
        # choose a new source location and velocity model for each trial
        self.qdep2[0] = self.qdep
        self.index[0] = 1
        for nm in range(1,self.nmc):
            val = ran_norm()
            self.qdep2[nm] = abs(self.qdep + self.sez * val) # randomly perturbed source depth
            self.index[nm] = (nm % self.ntab) + 1  # index used to choose velocity model
            
    def calculate_takeoff_angles(self):
        """
        Use HASH fortran subroutine to calulate takeoff angles for each trial
        """
        # loop over k picks b/c I broke it out -- NOTE: what does iflag do?
        #
        # find azimuth and takeoff angle for each trial
        for k in range(self.npol):
            for nm in range(self.nmc):
                self.p_azi_mc[k,nm] = self.qazi[k]
                self.p_the_mc[k,nm], iflag = get_tts(self.index[nm],self.dist[k],self.qdep2[nm])
        self.magap, self.mpgap = get_gap(self.p_azi_mc[:self.npol,0],self.p_the_mc[:self.npol,0],self.npol)
    
    def view_polarity_data(self):
        """
        Print out a list of polarity data for interactive runs
        """
        for k in range(self.npol):
            print('{0}   {1} {2} {3} {4}'.format(k,self.sname[k],self.p_azi_mc[k,0],self.p_the_mc[k,0],self.p_pol[k]))
    
    def check_minimum_polarity(self):
        """
        Polarity check
        """
        if self.npol >= self.npolmin:
            return True
        else:
            return False
    
    def check_maximum_gap(self):
        """
        Gap check
        """
        if ((self.magap > self.max_agap) or (self.mpgap > self.max_pgap)):
            return False
        else:
            return True
    
    def calculate_hash_focalmech(self, use_amplitudes=False):
        """
        Run the actual focal mech calculation, and find the probable mech

        use_amplitudes : bool of whether to use the 'focalamp_mc' routine
        """
        if use_amplitudes:
            # determine maximum acceptable number misfit polarities
            nmismax = max(int(self.nppl * self.badfrac),2)        # nint
            nextra  = max(int(self.nppl * self.badfrac * 0.5),2)  # nint
            qmismax = max(int(self.nspr * self.qbadfrac),2)        # nint
            qextra  = max(int(self.nspr * self.qbadfrac * 0.5),2)  # nint
            # find the set of acceptable focal mechanisms for all trials
            self.nf2, self.strike2, self.dip2, self.rake2, self.f1norm, self.f2norm = focalamp_mc(self.p_azi_mc, self.p_the_mc, self.sp_ratio[:self.npol], self.p_pol[:self.npol], self.nmc, self.dang2, self.nmax0, nextra, nmismax, self.qextra, self.qmismax, self.npol)
        else:
            # determine maximum acceptable number misfit polarities
            nmismax = max(int(self.npol * self.badfrac),2)        # nint
            nextra  = max(int(self.npol * self.badfrac * 0.5),2)  # nint
            
            # find the set of acceptable focal mechanisms for all trials
            self.nf2, self.strike2, self.dip2, self.rake2, self.f1norm, self.f2norm = focalmc(self.p_azi_mc, self.p_the_mc, self.p_pol[:self.npol], self.p_qual[:self.npol], self.nmc, self.dang2, self.nmax0, nextra, nmismax, self.npol)
        
        self.nout2 = min(self.nmax0, self.nf2)  # number mechs returned from sub
        self.nout1 = min(self.maxout, self.nf2) # number mechs to return
        
        # find the probable mechanism from the set of acceptable solutions
        self.nmult, self.str_avg, self.dip_avg, self.rak_avg, self.prob, self.var_est = mech_prob(self.f1norm[:,:self.nout2], self.f2norm[:,:self.nout2], self.cangle, self.prob_max, self.nout2) # nout2
    
    def calculate_quality(self, use_amplitudes=False):
        """
        Do the quality calulations
        """
        for imult in range(self.nmult):
            self.var_avg[imult] = (self.var_est[0,imult] + self.var_est[1,imult]) / 2.
            
            if use_amplitudes:
                self.mfrac[imult], self.mavg[imult], self.stdr[imult] = get_misf_amp(self.p_azi_mc[:self.npol,0], self.p_the_mc[:self.npol,0], self.sp_ratio[:self.npol], self.p_pol[:self.npol], self.str_avg[imult], self.dip_avg[imult], self.rak_avg[imult], self.npol) # npol
            else:
                # find misfit for prefered solution
                self.mfrac[imult], self.stdr[imult] = get_misf(self.p_azi_mc[:self.npol,0], self.p_the_mc[:self.npol,0], self.p_pol[:self.npol], self.p_qual[:self.npol], self.str_avg[imult], self.dip_avg[imult], self.rak_avg[imult], self.npol) # npol
            
            # HASH default solution quality rating
            if ((self.prob[imult] > 0.8) and (self.var_avg[imult] < 25) and (self.mfrac[imult] <= 0.15) and (self.stdr[imult] >= 0.5)):
                self.qual[imult]='A'
            elif ((self.prob[imult] > 0.6) and (self.var_avg[imult] <= 35) and (self.mfrac[imult] <= 0.2) and (self.stdr[imult] >= 0.4)):
                self.qual[imult]='B'
            elif ((self.prob[imult] > 0.5) and (self.var_avg[imult] <= 45) and (self.mfrac[imult] <= 0.3) and (self.stdr[imult] >= 0.3)):
                self.qual[imult]='C'
            else:
                self.qual[imult]='D'
    
    @property
    def _best_quality_index(self):
        """
        Returns index of highest quality solution
        ( imult of nmult in HASH-speak )
        
        Just use the A-D quality for now, could use
        (self.prob self.var_avg, self.mfrac, self.stdr)
        to make your own quality assessment...
        """
        # todo: make more sophisticated "Best" function
        return self.qual[:self.nmult].argsort()[0]


    def driver2(self, check_for_minimum_picks=True, check_for_maximum_gap_size=True):
        """
        Convenience method for hash_driver2
        
        This approximately acts like the "hash_driver2.f" program in the original HASH code.
        One must first input data from one's chosen sources. This method will generate tables,
        trial data, and takeoff angles from velocity models, don't use it if you make your own.
        
        """
        # Generate preliminary data for run
        self.load_velocity_models() # file list in 'self.vmodels'
        if not self.ntab:
            raise RuntimeWarning("No velocity tables loaded, continuing would be futile!")
        self.generate_trial_data()
        self.calculate_takeoff_angles()
        
        pass1 = self.check_minimum_polarity()
        pass2 = self.check_maximum_gap()
        
        # If it passes checks, run HASH
        if check_for_minimum_picks and not pass1:
            raise ValueError("Didn't pass check: # picks = {0} | Minimum = {1}".format(self.npol, self.npolmin))
        if check_for_maximum_gap_size and not pass2:
            raise ValueError("Didn't pass check: agap/pgap = {0}/{1} | Max allowed = {2}/{3}".format(self.magap, self.mpgap, self.max_agap, self.max_pgap))

        self.calculate_hash_focalmech()
        self.calculate_quality()
    
    def driver3(self, check_for_minimum_picks=True, check_for_maximum_gap_size=True):
        """
        Convenience method for hash_driver3
        
        This approximately acts like the "hash_driver3.f" program in the original HASH code.
        One must first input data from one's chosen sources. This method will generate tables,
        trial data, and takeoff angles from velocity models, don't use it if you make your own.
        
        NOT TESTED!!
        """
        # Generate preliminary data for run
        self.load_velocity_models() # file list in 'self.vmodels'
        if not self.ntab:
            raise RuntimeWarning("No velocity tables loaded, continuing would be futile!")
        self.generate_trial_data()
        self.calculate_takeoff_angles()
        
        pass1 = self.check_minimum_polarity()
        pass2 = self.check_maximum_gap()
        
        # If it passes checks, run HASH
        if check_for_minimum_picks and not pass1:
            raise ValueError("Didn't pass check: # picks = {0} | Minimum = {1}".format(self.npol, self.npolmin))
        if check_for_maximum_gap_size and not pass2:
            raise ValueError("Didn't pass check: agap/pgap = {0}/{1} | Max allowed = {2}/{3}".format(self.magap, self.mpgap, self.max_agap, self.max_pgap))

        self.calculate_hash_focalmech(use_amplitudes=True)
        self.calculate_quality(use_amplitudes=True)


class HashError(StandardError):
    """Throw this if something happens while running"""
    pass

