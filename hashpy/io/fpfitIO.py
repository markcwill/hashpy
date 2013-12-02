# -*- coding: utf-8 -*-

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
import numpy as np
import datetime

degrad = 180. / np.pi

def get_sta_coords(f_sta):
	'''
    Stub which reads in the HASH example station file
    
    :param f_sta: file handle with station info

    :returns: dict of list of floats -> { 'station_code' : [lat, lon, elv] }
    
    '''
	coords = {}
	for s in f_sta:
		sta = s[0:4]
		lat = float(s[42:50])
		lon = float(s[51:61])
		elv = float(s[62:67])/1000.
		if sta not in coords:
			coords[sta] = [lat,lon,elv]
	return coords
    
def check_polarity_file(polarity_file, sta_name, year, month, day, hour):
    """
    Stub for checking polarity file FPFIT style, there's a fortran sub for this,
    ideally, rewrite, but this will work for now...
    
    :param str polarity_file: name of polarity file in rando SC esoteric anti-format

    """
    from hashpy import libhashpy
    return libhashpy.check_pol(polarity_file, sta_name, year, month, day, hour)

def input(hp, files={} ):
    """
    Input function for FPFIT-style HASH

    :param dict files: dict of str filenames

    Input dictionary has keys: ('station', 'polarity', 'input')
    for station file, polarity reversal file, and main event/phase
    input, respectively.

    """
    with open(files['station']) as fs:
        sta_coords = get_sta_coords(fs)
    
    with open(files['input']) as fi:
    
        while True:
            # read in earthquake location, etc SCEDC format
            line = fi.readline()
            if not line:
                break
            
            iyr    = int(line[0:4])
            imon   = int(line[4:6])
            idy    = int(line[6:8])
            ihr    = int(line[8:10])
            imn    = int(line[10:12])
            sec    = line[12:17]
            ilatd  = int(line[17:19])
            cns    = line[19]
            qlatm  = float(line[20:25])
            ilond  = int(line[25:28])
            cew    = line[28]
            qlonm  = float(line[29:34])
            
            hp.qdep   = float(line[34:39]) # then 49x
            hp.seh    = float(line[88:93]) # then 1x
            hp.sez    = float(line[94:99]) # then 40x
            hp.qmag   = float(line[139:143]) # then 6x
            hp.icusp  = int(line[149:165])
            
            qsec = float(sec)
            secs = sec.split('.')
            isec = int(secs[0])
            imsec = int(secs[1]+'0000')
            dt = datetime.datetime(iyr,imon,idy,ihr,imn,isec,imsec)
            hp.tstamp = (dt - datetime.datetime(1970, 1, 1)).total_seconds()

            hp.qlat = ilatd + (qlatm / 60.0)
            if ('S' in cns):
                hp.qlat *= -1
            hp.qlon = -(ilond + (qlonm / 60.0))
            if ('E' in cew):
                hp.qlon *= -1
            aspect = np.cos(hp.qlat / degrad)
            
            # set parameters not given in input file
            if not hp.sez:
                hp.sez = 1.
            terr = -9
            rms = -9
            nppick = -9
            nspick = -9
            evtype = 'L'
            magtype = 'X'
            locqual = 'X'
            
            # read in polarities - SCEDC format - ** YOU MAY NEED TO CHANGE THE INPUT FORMAT **
            k = 0
            while True:
                ph = fi.readline()
                hp.sname[k]     = ph[0:4]
                hp.snet[k]      = ph[5:7]
                hp.scomp[k]     = ph[9:12]
                hp.pickonset[k] = ph[13]
                hp.pickpol[k]   = ph[15]
                hp.arid[k]      = k

                if (hp.sname[k] == '    '):
                    break
                
                #flat,flon,felv = getstat_tri(stfile,sname[k],scomp[k],snet[k])
                # SCSN station information - ** YOU MAY NEED TO USE YOUR OWN SUBROUTINE **
                #
                #-- the above sucks, I needed to use my own subroutine --#
                # 
                if hp.sname[k] in sta_coords:
                    hp.flat[k], hp.flon[k], hp.felv[k] = sta_coords[hp.sname[k]]
                else:
                    continue
                
                dx = (hp.flon[k] - hp.qlon) * 111.2 * aspect
                dy = (hp.flat[k] - hp.qlat) * 111.2
                hp.dist[k] = np.sqrt(dx**2 + dy**2)
                hp.qazi[k] = 90. - np.arctan2(dy,dx) * degrad
                
                if (hp.qazi[k] < 0.):
                    hp.qazi[k] += 360.
                
                if (hp.dist[k] > hp.delmax):
                    continue
                
                if (hp.pickpol[k] in 'Uu+'):
                    hp.p_pol[k] = 1
                elif (hp.pickpol[k] in 'Dd-'):
                    hp.p_pol[k] = -1
                else:
                    continue
                
                if (hp.pickonset[k] in 'Ii'):
                    hp.p_qual[k] = 0
                else:
                    hp.p_qual[k] = 1
               
                # Try to check a polarity file using sub
                try:
                    spol = check_polarity_file(files['polarity'], hp.sname[k], iyr, imon, idy, ihr)
                except:
                    spol = 1
                finally:
                    hp.p_pol[k] *= spol
                
                k += 1
                continue
            hp.npol = k - 1
        
