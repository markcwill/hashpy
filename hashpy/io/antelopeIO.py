# -*- coding: utf-8 -*-
#
# dbhashpype.py
#
# by Mark Williams 2012.313
#
# Class to run HASHpy using Antelope database I/O
#
import os.path
import numpy as np
from antelope.datascope import Dbptr, dbtmp, dblookup, dbprocess

# for HASH compatiblity, change later.
degrad = 180. / np.pi
rad = 1. / degrad

class RowPointerDict(dict):
    
    _dbptr = None

    def __init__(self, dbptr=None, record=None):
        self._dbptr = Dbptr(dbptr)
        if record is not None:
            self._dbptr.record = record
        if self._dbptr.record < 0:
            self._dbptr.record = 0

    def __getitem__(self, key):
        return self._dbptr.getv(key)[0]

    def __setitem__(self, key, value):
        self._dbptr.putv(key, value)

    def __len__(self):
        return self._dbptr.nrecs()



def input(hp, dbname, evid=None, orid=None):
    '''Input HASH data from Antelope database
    
    This will accept a database name OR Antelope Dbptr, and either
    an ORID, or an EVID, in which case the 'prefor' ORID is looked
    up and used.
    
    Inputs
    ------
    dbname	:	str or antelope.datascope.Dbptr
    orid	:	int of ORID
    evid	:	int of EVID
    '''
    
    db = Dbptr(dbname)

    if orid is None:
        dbv = dbprocess(db,['dbopen event', 'dbsubset evid == '+str(evid)])
        orid = RowPointerDict(dbv)['prefor']
    
    db = dbprocess(db,[ 'dbopen origin', 'dbsubset orid == '+str(orid),
                    'dbjoin origerr', 'dbjoin assoc',  'dbjoin arrival',
                    'dbjoin affiliation', 'dbjoin site',
                    'dbsubset iphase =~ /.*[Pp].*/',
                    'dbsubset (ondate <= time)',
                    'dbsubset (time <= offdate) || (offdate == -1)']
                    )
    
    ph = RowPointerDict(db, record=0)

    hp.nrecs = len(ph)
    if len(ph) > 0:
        raise ValueError("No picks for this ORID: {0}".format(orid) )
        
    hp.tstamp = ph['origin.time']
    hp.qlat   = ph['origin.lat']
    hp.qlon   = ph['origin.lon']
    hp.qdep   = ph['origin.depth']
    hp.qmag   = ph['origin.ml']
    hp.icusp  = ph['origin.orid']
    hp.seh	  = ph['origerr.smajax']
    hp.sez	  = ph['origerr.sdepth']
    
    aspect = np.cos(hp.qlat / degrad) # convert using python later.
    
    # The index 'k' is deliberately non-Pythonic to deal with the fortran
    # subroutines which need to be called and the structure of the original HASH code.
    # May be able to update with a rewrite... YMMV
    k = 0
    for n in range(len(ph)):
        
        ph = RowPointerDict(db, record=n)

        hp.sname[k]	  = ph['sta']
        hp.snet[k]	  = ph['net']
        hp.scomp[k]	  = ph['chan']
        hp.pickonset[k] = 'I'
        hp.pickpol[k]   = ph['fm']
        hp.arid[k]	  = ph['arid']
        
        flat, flon, felv = ph['site.lat'],ph['site.lon'],ph['site.elev']
        hp.esaz[k] = ph['esaz']

        
        # dist @ azi, get from db OR obspy or another python mod (antelope) could do this on WGS84
        dx = (flon - hp.qlon) * 111.2 * aspect
        dy = (flat - hp.qlat) * 111.2
        dist = np.sqrt(dx**2 + dy**2)
        qazi = 90. - np.arctan2(dy,dx) * degrad
        
        if (qazi < 0.):
            qazi = qazi + 360.
        if (dist > hp.delmax):
            continue
        if (hp.pickpol[k] in 'CcUu'):
            hp.p_pol[k] = 1
        elif (hp.pickpol[k] in 'RrDd'):
            hp.p_pol[k] = -1
        else:
            continue
        
        # save them for other functions -MCW
        hp.dist[k] = dist
        hp.qazi[k] = qazi
        hp.flat[k] = flat
        hp.flon[k] = flon
        hp.felv[k] = felv
        
        if (hp.pickonset[k] in 'Ii'):
            hp.p_qual[k] = 0
        else:
            hp.p_qual[k] = 1
        
        # polarity check in original code... doesn't work here
        #hp.p_pol[k] = hp.p_pol[k] * hp.spol
        k += 1
    #npol = k - 1
    hp.npol = k # k is zero indexed in THIS loop
    db.close()


def output(hp, dbout=None, solution=0, schema="css3.0"):
    '''Write the preferred HASH solution to Datascope database.
    
    This writes the strike, dip, rakes to 'fplane', arids used for a
    given mech in 'predmech' and the takeoffs in 'predarr'.
    
    Input
    -----
    dbout	: str or antelope.datascope.Dbptr to database
    solution : <STUB> int of desired solution.
    
    '''
    from hashpy.doublecouple import DoubleCouple

    x = solution
    dc = DoubleCouple([hp.str_avg[x], hp.dip_avg[x], hp.rak_avg[x]])
    str1, dip1, rak1 = dc.plane1 
    str2, dip2, rak2 = dc.plane2
    axes = dc.axis
    
    if dbout is not None:
        db = Dbptr(dbout, perm='r+')
    else:
        db = dbtmp(schema)

    mechid = db.nextid('mechid')
    
    dbfpln = dblookup(db,table='fplane')
    
    dbfpln.record = dbfpln.addnull()
    dbfpln.putv(
            'orid' , hp.icusp,
            'str1' , round(str1,1) ,
            'dip1' , round(dip1,1) ,
            'rake1', round(rak1,1),
            'algorithm', "HASH",
            'mechid', mechid,
            'str2' , round(str2,1) ,
            'dip2' , round(dip2,1) ,
            'rake2', round(rak2,1),
            'taxazm', round(axes['T']['azimuth'])
            'taxplg', round(axes['T']['dip'])
            'paxazm', round(axes['P']['azimuth'])
            'paxplg', round(axes['P']['dip'])
            )
    
    dbpmec = dblookup(db,table='predmech')
    dbparr = dblookup(db,table='predarr')
    for k in range(hp.npol):
        if hp.p_pol[k] > 0:
            fm = 'U'
        else:
            fm = 'D'
        dbpmec.record = dbpmec.addnull()
        dbpmec.putv('arid', int(hp.arid[k]) ,
                    'orid', hp.icusp,
                    'mechid', mechid,
                    'fm', fm,
                    )
        dbparr.record = dbparr.addnull()
        dbparr.putv('arid', int(hp.arid[k]),
                    'orid', hp.icusp, 
                    'esaz', hp.qazi[k], 
                    'dip' , hp.p_the_mc[k,0],
                    )
    return db


#
# Utility for AntelopeIO, for sucking in a pf file
#
def load_pf(hp, pffile='dbhash.pf'):
    '''Update runtime settings from a pf file
    
    This inputs from pf and converts the original HASH command line params to
    proper python type. One can also specify names of velocity model files in
    the pf file. See an actual file for example.
    
    Right now these settings are inherited from the HashPype class,
     and are not instance attributes.
     
     Input
     -----
     pffile : string of full path to pf file
     
    '''
    # Version change 5.2->5.3 broke Antelope API, quick fix for now
    # TODO: 1) use Mark's version agnostic 'pfgetter' function
    #       2) man up and use JSON for your configs
    try:
        from antelope.stock import pfget
    except ImportError:
        from antelope.stock import pfread as pfget
    
    pf_settings = pfget(pffile)
    
    # Little hack to do type conversions 
    for key in pf_settings:
        pfi = pf_settings[key]
        if key in ['badfrac','prob_max']:
            pfi = float(pfi)
        elif key in ['npolmin','max_agap','max_pgap','dang','nmc','maxout', 'delmax','cangle']:
            pfi = int(pfi)
        else:
            pass
        hp.__setattr__(key, pfi)
    
    if 'vmodel_dir' in pf_settings and 'vmodels' in pf_settings:
        hp.vmodels = [os.path.join(hp.vmodel_dir, table) for table in hp.vmodels]


