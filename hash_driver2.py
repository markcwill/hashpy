#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# hash_driver2.py
# -by Mark Williams 2012
# This is a port of the file 'hash_driver2.f' rewritten in python.
#
# This is just for compatibility and testing. See the 'driver.py'
# script for creating a Hash run pipeline object.
#
# ----------------------------------------------------------------
# Sample main driver program for using the focal mechanism subroutines.  
# The polarity input format is the "phase" format of the SCEDC, similar to FPFIT.
# Takeoff angle uncertainty is expresses as a suite of 1D velocity models.
#
# INPUT FILE: "fpfile" FPFIT-like format - ** YOU MAY NEED TO CHANGE THE INPUT FORMAT **
#
# event line:
#     columns  format   value
#     -------------------------
#     1-4        i4     origin time, year
#     5-12       4i2    origin time, month, day, hour, minute
#     13-17      f5.2   origin time, seconds
#     18-19      i2     latitude, degrees
#     20         a1     latitude, 'S'=south
#     21-25      f5.2   latitude, minutes
#     26-28      i3     longitude, degrees
#     29         a1     longitude, 'E'=east
#     30-34      f5.2   longitude, minutes
#     35-39      f5.2   depth, km
#     89-93      f5.2   horizontal location uncertainty, km
#     95-99      f5.2   vertical location uncertainty, km
#     140-143    f4.2   magnitude
#     150-165    a16    event ID
#
# polarity line:
#     columns  format   value
#     -------------------------
#     1-4        a4     station name
#     6-7        a2     station network
#     10-12      a3     station component
#     14         a1     P onset, I or E
#     16         a1     P polarity : U, u, +, D, d, or -
#
from sys import argv
from libhashpy import *
import numpy as np
from hash_utils import fortran_include, get_sta_coords, test_stereo


npick0, nmc0, nmax0 = fortran_include('src/param.inc')
dang0, ncoor        = fortran_include('src/rot.inc')

# initialize arrays

# input arrays
sname     = np.empty(npick0, 'a4', 'F')
scomp     = np.empty(npick0, 'a3', 'F')
snet      = np.empty(npick0, 'a2', 'F')
pickpol   = np.empty(npick0, 'a1', 'F')
pickonset = np.empty(npick0, 'a1', 'F')
p_pol     = np.empty(npick0, int, 'F')
p_qual    = np.empty(npick0, int, 'F')
spol      = np.empty(npick0, int, 'F')
p_azi_mc  = np.empty((npick0,nmc0), float, 'F')
p_the_mc  = np.empty((npick0,nmc0), float, 'F')
index     = np.empty(nmc0, int, 'F')
qdep2     = np.empty(nmc0, float, 'F')

# Output arrays
f1norm  = np.empty((3,nmax0), float, 'F')
f2norm  = np.empty((3,nmax0), float, 'F')
strike2 = np.empty(nmax0, float, 'F')
dip2    = np.empty(nmax0, float, 'F')
rake2   = np.empty(nmax0, float, 'F')
str_avg = np.empty(5, float, 'F')
dip_avg = np.empty(5, float, 'F')
rak_avg = np.empty(5, float, 'F')
var_est = np.empty((2,5), float, 'F')
var_avg = np.empty(5, float, 'F')
mfrac   = np.empty(5, float, 'F')
stdr    = np.empty(5, float, 'F')
prob    = np.empty(5, float, 'F')
qual    = np.empty(5, 'a', 'F')

degrad = 180. / np.pi
rad = 1. / degrad
  
summary = '''
c variables for storing earthquake input information  
      character*16 icusp            ! event ID
      real qlat,qlon,qdep           ! location
      real qmag                     ! magnitude
      integer iyr,imon,idy,ihr,imn  ! origin time, year, month, day, hour, minute
      real qsec                     ! origin time, seconds
      real seh, sez                 ! location error, horizontal & vertical 
      real rms                      ! residual travel time error 
      real terr                     ! origin time error 
      character*1 evtype            ! event type
      character*1 magtype           ! magnitude type
      character*1 locqual           ! location quality
      character*1 cns,cew           ! north/south and east/west codes

c variables for polarity information, input to HASH subroutines
      character*4 sname(npick0)                        ! station name
      character*3 scomp(npick0)                        ! station component
      character*2 snet(npick0)                         ! network code
      character*1 pickpol(npick0),pickonset(npick0)    ! polarity pick : U, u, +, D, d, or - ; onset, I or E
      integer p_pol(npick0),p_qual(npick0),spol        ! polarity pick (+/-1),  quality (0/1), and reversal (+/-1)
      real p_azi_mc(npick0,nmc0),p_the_mc(npick0,nmc0) ! azimuth and takeoff angle for each trial
      integer index(nmc0)                              ! index into velocity models, for each trial
      real qdep2(nmc0)                                 ! new depth, for each trail
      integer nmc                                      ! number of trails with different azimuth and take-off angles
      integer npol                                     ! number of polarity readings

c variables for set of acceptable mechanisms, output of HASH subroutines
      integer nout2                                    ! number of acceptable mechanisms returned
      integer nmult                                    ! number of solutions (1 if no mulitples)
      real str_avg(5),dip_avg(5),rak_avg(5)            ! solution(s)
      real f1norm(3,nmax0),f2norm(3,nmax0)             ! normal vectors to the fault planes
      real strike2(nmax0),dip2(nmax0),rake2(nmax0)     ! strike, dip and rake
      real var_est(2,5),var_avg(5)                     ! variance of each plane, average
      real mfrac(5),stdr(5)                            ! fraction misfit polarities, station distribution
      real prob(5)                                     ! probability true mechanism is "close" to preferred solution
      character*1 qual(5),mflag                        ! solution quality rating, multiple flag
c
c control parameters
      integer npolmin,max_agap,max_pgap                ! minimum polarities, maximum azimuthal & "plungal" gap
      real delmax                                      ! maximum station distance
      real dang,dang2                                  ! grid angle to sample focal sphere
      integer maxout                                   ! max number of acceptable mechanisms output
      real badfrac                                     ! assumed rate of polarity error (fraction)
      real cangle                                      ! definition of "close" == 45 degrees
'''

interactive = True
# Parse the input file, instead of interactively...
if not interactive:
    input_file = open(argv[1]).read()
    inp_vars = input_file.split()
    (stfile,plfile,fpfile,outfile1,outfile2,npolmin,max_agap,max_pgap,
      dang,nmc,maxout,badfrac,delmax,cangle,prob_max) = inp_vars[:15]
    num_vel_mods = int(inp_vars[15])
    vmfile = inp_vars[16:16+num_vel_mods]
    for i in range(num_vel_mods):
        ntab = mk_table_add(i+1,vmfile[i])  
else:
    stfile   = raw_input('Enter station list file: ')
    plfile   = raw_input('Enter station polarity reversal file: ')
    fpfile   = raw_input('Enter name of input file (FPFIT-like format): ')
    outfile1 = raw_input('Enter output file name for focal mechanisms: ')
    outfile2 = raw_input('Enter output file name for acceptable planes: ')
    npolmin  = raw_input('Enter mininum number of polarities (e.g., 8) :')
    max_agap = raw_input('Enter maximum azimuthal gap (e.g., 90): ')
    max_pgap = raw_input('Enter maximum takeoff angle gap (e.g., 60): ')
    dang     = raw_input('Enter grid angle for focal mech search, in degrees (min {0}): '.format(dang0))
    nmc      = raw_input('Enter number of trials (e.g., 30): ')
    maxout   = raw_input('Enter maxout for focal mech. output (e.g., 500): ')
    badfrac  = raw_input('Enter fraction of picks presumed bad (e.g., 0.10): ')
    delmax   = raw_input('Enter maximum allowed source-station distance, in km (e.g., 120): ')
    cangle   = raw_input('Enter angle for computing mechanisms probability, in degrees (e.g., 45): ')
    prob_max = raw_input('Enter probability threshold for multiples (e.g., 0.1): ')
    ntab = 10
    ntab = mk_table(ntab)  


sta_coords = get_sta_coords(stfile)
fh12 = open(fpfile)
fh13 = open(outfile1, 'w')
fh11 = open(outfile2, 'w')
npolmin = int(npolmin)
max_agap = int(max_agap)
max_pgap = int(max_pgap)
dang = float(dang)
dang2=max(dang0,dang) # don't do finer than dang0
nmc = int(nmc)
maxout = int(maxout)
badfrac = float(badfrac)
delmax = float(delmax)
cangle = float(cangle)
prob_max = float(prob_max)


# make tables of takeoff angles for various velocity models


# start of event loop - not necessary for one event?
#20    continue

while True:
    # read in earthquake location, etc SCEDC format - ** YOU MAY NEED TO CHANGE THE INPUT FORMAT **
    line = fh12.readline()
    if not line:
        break
    
    iyr    = int(line[0:4])
    imon   = int(line[4:6])
    idy    = int(line[6:8])
    ihr    = int(line[8:10])
    imn    = int(line[10:12])
    qsec   = float(line[12:17])
    ilatd  = int(line[17:19])
    cns    = line[19]
    qlatm  = float(line[20:25])
    ilond  = int(line[25:28])
    cew    = line[28]
    qlonm  = float(line[29:34])
    qdep   = float(line[34:39]) # then 49x
    seh    = float(line[88:93]) # then 1x
    sez    = float(line[94:99]) # then 40x
    qmag   = float(line[139:143]) # then 6x
    icusp  = int(line[149:165])
    
    qlat = ilatd + (qlatm / 60.0)
    if ('S' in cns):
        qlat *= -1
    qlon = -(ilond + (qlonm / 60.0))
    if ('E' in cew):
        qlon *= -1
    aspect = np.cos(qlat / degrad) # convert using python later.
    
    # set parameters not given in input file
    if not sez:
        sez = 1.
    terr = -9
    rms = -9
    nppick = -9
    nspick = -9
    evtype = 'L'
    magtype = 'X'
    locqual = 'X'
    
    # choose a new source location and velocity model for each trial
    qdep2[0] = qdep
    index[0] = 1
    for nm in range(1,nmc):
        val = ran_norm()
        qdep2[nm] = qdep + sez * val # randomly perturbed source depth
        index[nm] = (nm % ntab) + 1  # index used to choose velocity model
        
    # read in polarities - SCEDC format - ** YOU MAY NEED TO CHANGE THE INPUT FORMAT **
    k = 0
    #30    continue
    while True:
        ph = fh12.readline()
        sname[k]     = ph[0:4]
        snet[k]      = ph[5:7]
        scomp[k]     = ph[9:12]
        pickonset[k] = ph[13]
        pickpol[k]   = ph[15]
        #read (12,35,end=40) sname(k),snet(k),scomp(k),pickonset(k),pickpol(k)
        #format (a4,1x,a2,2x,a3,1x,a1,1x,a1)
        
        if (sname[k] == '    '):
            #goto 40 ! end of data for this event
            break
        
        #flat,flon,felv = getstat_tri(stfile,sname[k],scomp[k],snet[k])
        # SCSN station information - ** YOU MAY NEED TO USE YOUR OWN SUBROUTINE **
        #
        #-- the above sucks, I needed to use my own subroutine --#
        # 
        if sname[k] in sta_coords:
            flat,flon,felv = sta_coords[sname[k]]
        else:
            continue
        
        #print '{0} {1} {2} {3} {4} {5} {6} {7}'.format(k, sname[k], snet[k], scomp[k], pickonset[k], pickpol[k], flat, flon)
        # obspy or another python mod (antelope) could do this on WGS84
        dx = (flon - qlon) * 111.2 * aspect
        dy = (flat - qlat) * 111.2
        dist = np.sqrt(dx**2 + dy**2)
        qazi = 90. - np.arctan2(dy,dx) * degrad
        if (qazi < 0.):
            qazi = qazi + 360.
        if (dist > delmax):
            continue
        if (pickpol[k] in 'Uu+'):
            p_pol[k] = 1
        elif (pickpol[k] in 'Dd-'):
            p_pol[k] = -1
        else:
            continue
        
        if (pickonset[k] in 'Ii'):
            p_qual[k] = 0
        else:
            p_qual[k] = 1
        
        spol = check_pol(plfile,sname[k],iyr,imon,idy,ihr)  # SCSN station polarity reversal information - ** YOU MAY NEED TO USE YOUR OWN SUBROUTINE **
        p_pol[k] = p_pol[k] * spol
        
        # find azimuth and takeoff angle for each trial
        for nm in range(nmc):
            p_azi_mc[k,nm] = qazi
            p_the_mc[k,nm], iflag = get_tts(index[nm],dist,qdep2[nm])
        
        k += 1
        continue
    npol = k - 1
    
    # view polarity data
    for k in range(npol):
        print '{0}   {1} {2} {3} {4}'.format(k,sname[k],p_azi_mc[k,0],p_the_mc[k,0],p_pol[k])
    
    # stop if there aren't enough polarities
    print 'cid = {0}  npol = {1}'.format(icusp,npol)
    if (npol < npolmin):
        str_avg[0] = 999
        dip_avg[0] = 99
        rak_avg[0] = 999
        var_est[0,0] = 99
        var_est[1,0] = 99
        mfrac[0] = 0.99
        qual[0] = 'F'
        prob[0] = 0.0
        nout1 = 0
        nout2 = 0
        nmult = 0
        #goto 400
        
    if str_avg[0] == 999:
        pass
    else:   
        # determine maximum azimuthal and takeoff gap in polarity observations and stop if either gap is too big
        magap,mpgap = get_gap(p_azi_mc[:npol,0],p_the_mc[:npol,0],npol)
        if ((magap > max_agap) or (mpgap > max_pgap)):
            str_avg[0] = 999
            dip_avg[0] = 99
            rak_avg[0] = 999
            var_est[0,0] = 99
            var_est[1,0] = 99
            mfrac[0] = 0.99
            qual[0] = 'E'
            prob[0] = 0.0
            nout1 = 0
            nout2 = 0
            nmult = 0
            #goto 400
    
    if str_avg[0] == 999:
        pass
    else:
        # determine maximum acceptable number misfit polarities
        nmismax = max(int(npol * badfrac),2)        # nint
        nextra  = max(int(npol * badfrac * 0.5),2)  # nint
        
        # find the set of acceptable focal mechanisms for all trials
        nf2,strike2,dip2,rake2,f1norm,f2norm = focalmc(p_azi_mc,p_the_mc,p_pol[:npol],p_qual[:npol],nmc,dang2,nmax0,nextra,nmismax,npol)
        nout2 = min(nmax0,nf2)  # number mechs returned from sub
        nout1 = min(maxout,nf2) # number mechs to return
        
        # find the probable mechanism from the set of acceptable solutions
        nmult,str_avg,dip_avg,rak_avg,prob,var_est = mech_prob(f1norm[:,:nout2],f2norm[:,:nout2],cangle,prob_max,nout2) # nout2
        
        for imult in range(nmult):
            var_avg[imult] = (var_est[0,imult] + var_est[1,imult]) / 2.
            print 'cid = {0} {1}  mech = {2} {3} {4}'.format(icusp,imult,str_avg[imult],dip_avg[imult],rak_avg[imult])
            # find misfit for prefered solution
            mfrac[imult],stdr[imult] =  get_misf(p_azi_mc[:npol,0],p_the_mc[:npol,0],p_pol[:npol],p_qual[:npol],str_avg[imult],dip_avg[imult],rak_avg[imult],npol) # npol
            
            # solution quality rating  ** YOU MAY WISH TO DEVELOP YOUR OWN QUALITY RATING SYSTEM **
            if ((prob[imult] > 0.8) and (var_avg[imult] < 25) and (mfrac[imult] <= 0.15) and (stdr[imult] >= 0.5)):
                qual[imult]='A'
            elif ((prob[imult] > 0.6) and (var_avg[imult] <= 35) and (mfrac[imult] <= 0.2) and (stdr[imult] >= 0.4)):
                qual[imult]='B'
            elif ((prob[imult] > 0.5) and (var_avg[imult] <= 45) and (mfrac[imult] <= 0.3) and (stdr[imult] >= 0.3)):
                qual[imult]='C'
            else:
                qual[imult]='D'
    
    if (nmult > 1):
        mflag='*'
    else:
        mflag=' '
    
    format13 = ' '.join(['{{{0}}}'.format(s) for s in range(32)]) + '\n'
    for i in range(nmult):
        # output prefered mechanism  ** YOU MAY WISH TO CHANGE THE OUTPUT FORMAT **
        fh13.write(format13.format(icusp,iyr,imon,idy,ihr,imn,qsec,evtype,
                  qmag,magtype,qlat,qlon,qdep,locqual,rms,seh,sez,terr,
                  nppick+nspick,nppick,nspick,int(round(str_avg[i])),
                  int(round(dip_avg[i])),int(round(rak_avg[i])),
                  int(round(var_est[0,i])),int(round(var_est[1,i])),
                  npol,int(round(mfrac[i]*100.)),qual[i],
                  int(round(100*prob[i])),int(round(100*stdr[i])),mflag))
    
    # output set of acceptable mechanisms  ** YOU MAY WISH TO CHANGE THE OUTPUT FORMAT **
    format11 = ' '.join(['{{{0}}}'.format(s) for s in range(24)]) + '\n'
    fh11.write(format11.format(iyr,imon,idy,ihr,imn,qsec,qmag,qlat,qlon,
              qdep,sez,seh,npol,nout2,icusp,str_avg[0],dip_avg[0],
              rak_avg[0],var_est[0,0],var_est[1,0], mfrac[0],qual[0],
              prob[0],stdr[0]))
    format11 = ' '.join(['{{{0}}}'.format(s) for s in range(9)]) + '\n'
    for ic in range(nout1):
        fh11.write(format11.format(strike2[ic],dip2[ic],rake2[ic],
                   f1norm[0,ic], f1norm[1,ic],f1norm[2,ic],f2norm[0,ic],
                   f2norm[1,ic],f2norm[2,ic]))
    
    #test_stereo(p_azi_mc[:npol,0],p_the_mc[:npol,0],p_pol[:npol],sdr=[str_avg[0],dip_avg[0],rak_avg[0]])

    go_again = raw_input("Continue or stop (and graph)? [Y/n]: ")
    
    if go_again == '' or go_again in 'Yesyes':
        continue # next earthquake line
    else:
        test_stereo(p_azi_mc[:npol,0],p_the_mc[:npol,0],p_pol[:npol],sdr=[str_avg[0],dip_avg[0],rak_avg[0]])
        fh11.close()
        fh12.close()
        fh13.close()
        break
fh11.close()
fh12.close()
fh13.close()
#test_stereo(p_azi_mc[:npol,0],p_the_mc[:npol,0],p_pol[:npol],sdr=[str_avg[0],dip_avg[0],rak_avg[0]])



