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
    dbname  :   str or antelope.datascope.Dbptr
    orid    :   int of ORID
    evid    :   int of EVID
    '''
    
    db = Dbptr(dbname)

    if orid is None:
        dbv = dbprocess(db,['dbopen event', 'dbsubset evid == '+str(evid)])
        orid = RowPointerDict(dbv)['prefor']
    
    db = dbprocess(db,['dbopen origin', 'dbsubset orid == '+str(orid),
                    'dbjoin origerr', 'dbjoin assoc',  'dbjoin arrival',
                    'dbjoin affiliation', 'dbjoin site',
                    'dbsubset iphase =~ /.*[Pp].*/',
                    'dbsubset (ondate <= time)',
                    'dbsubset (time <= offdate) || (offdate == -1)']
                    )
    
    ph = RowPointerDict(db, record=0)

    hp.nrecs = len(ph)
    if len(ph) <= 0:
        raise ValueError("No picks for this ORID: {0}".format(orid) )
        
    hp.tstamp = ph['origin.time']
    hp.qlat = ph['origin.lat']
    hp.qlon = ph['origin.lon']
    hp.qdep = ph['origin.depth']
    hp.qmag = ph['origin.ml']
    hp.icusp = ph['origin.orid']
    hp.seh = ph['origerr.smajax']
    hp.sez = ph['origerr.sdepth']
    
    aspect = np.cos(hp.qlat / degrad) # convert using python later.
    
    # The index 'k' is deliberately non-Pythonic to deal with the fortran
    # subroutines which need to be called and the structure of the original HASH code.
    # May be able to update with a rewrite... YMMV
    k = 0
    for n in range(len(ph)):
        # Extract pick data from the db 
        ph = RowPointerDict(db, record=n)

        hp.sname[k] = ph['sta']
        hp.snet[k] = ph['net']
        hp.scomp[k] = ph['chan']
        hp.pickonset[k] = ph['qual'].strip('.')
        hp.pickpol[k] = ph['fm']
        hp.arid[k] = ph['arid']
        
        flat, flon, felv = ph['site.lat'],ph['site.lon'],ph['site.elev']
        hp.esaz[k] = ph['esaz']

        # Distance and Azimuth filtering
        dx = (flon - hp.qlon) * 111.2 * aspect
        dy = (flat - hp.qlat) * 111.2
        dist = np.sqrt(dx**2 + dy**2)
        qazi = 90. - np.arctan2(dy,dx) * degrad
        
        if (qazi < 0.):
            qazi = qazi + 360.
        
        if (dist > hp.delmax):
            continue

        # Try to get an up/down polarity
        if not hp.pickpol[k].lower():
            continue
        if (hp.pickpol[k].lower() in 'cu'):
            hp.p_pol[k] = 1
        elif (hp.pickpol[k].lower() in 'dr'):
            hp.p_pol[k] = -1
        else:
            continue

        # Save them for other functions
        hp.dist[k] = dist
        hp.qazi[k] = qazi
        hp.flat[k] = flat
        hp.flon[k] = flon
        hp.felv[k] = felv
        
        # Try to get the onset, impulsive if none
        if (hp.pickonset[k].lower() == 'i'):
            hp.p_qual[k] = 0
        elif (hp.pickonset[k].lower() == 'e'):
            hp.p_qual[k] = 1
        elif (hp.pickonset[k].lower() == 'w'):
            hp.p_qual[k] = 1
        else:
            hp.p_qual[k] = 0
        
        # polarity check in original code... doesn't work here
        #hp.p_pol[k] = hp.p_pol[k] * hp.spol
        k += 1
    hp.npol = k # k is zero indexed in THIS loop
    db.close()


def output(hp, dbout=None, solution=0, schema="css3.0"):
    '''Write the preferred HASH solution to Datascope database.
    
    This writes the strike, dip, rakes to 'fplane', arids used for a
    given mech in 'predmech' and the takeoffs in 'predarr'.
    
    Input
    -----
    dbout   : str or antelope.datascope.Dbptr to database
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
            'orid', hp.icusp,
            'str1', round(str1,1) ,
            'dip1', round(dip1,1) ,
            'rake1', round(rak1,1),
            'algorithm', "HASH",
            'mechid', mechid,
            'auth', 'hashpy:'+ hp.author,
            'str2', round(str2,1) ,
            'dip2', round(dip2,1) ,
            'rake2', round(rak2,1),
            'taxazm', round(axes['T']['azimuth'],1),
            'taxplg', round(axes['T']['dip'],1),
            'paxazm', round(axes['P']['azimuth'],1),
            'paxplg', round(axes['P']['dip'],1),
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
# Utility for AntelopeIO, for adding HASH parameters from a pf file
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
    #       2) man up and use JSON/ini/native python for your configs
    try:
        from antelope.stock import pfget
    except ImportError:
        from antelope.stock import pfread as pfget
    
    pf_settings = pfget(pffile)
    pf_keys = list(pf_settings.keys()) # b/c v5.3 broke backwards dict compat

    # Little hack to do type conversions 
    for key in pf_keys:
        pfi = pf_settings[key]
        if key in ['badfrac','prob_max']:
            pfi = float(pfi)
        elif key in ['npolmin','max_agap','max_pgap','dang','nmc','maxout', 'delmax','cangle']:
            pfi = int(pfi)
        else:
            pass
        hp.__setattr__(key, pfi)
    
    if 'vmodel_dir' in pf_keys and 'vmodels' in pf_keys:
        hp.vmodels = [os.path.join(hp.vmodel_dir, table) for table in hp.vmodels]

##############################################################################
# Antelope <-> ObsPy utilities
##############################################################################
#
# These are not necessary to run Antelope with HASH, for interactive/GUI
# scripts and housekeeping utils. 
#
def readANTELOPE(database, station=None, channel=None, starttime=None, endtime=None):
    """
    Reads a portion of a Antelope wfdisc table to a Stream.
    
    Attempts to return one Trace per line of the 'wfdisc' view passed.  
    Additionally, will filter and cut with respect to any of the fields
    in the primary key IF specified. (sta chan time::endtime)
    
    NOTE: Currently MUST have both times (start/end) or neither.
    the returned Traces will have a new attribute, 'db'

    :type database: string or antelope.datascope.Dbptr
    :param database: Antelope database name or pointer
    :type station: string
    :param station: Station expression to subset
    :type channel: string
    :param channel: Channel expression to subset
    :type starttime: :class: `~obspy.core.utcdatetime.UTCDateTime`
    :param starttime: Desired start time
    :type endtime: :class: `~obspy.core.utcdatetime.UTCDateTime`
    :param endtime: Desired end time
        
    :rtype: :class: `~obspy.core.stream.Stream'
    :return: Stream with one Trace for each row of the database view
    
    .. rubric:: Example
    
    >>> st = readANTELOPE('/opt/antelope/example/db', station='TOLO', channel='LH.',
                        starttime=UTCDateTime(2008,6,13), endtime=UTCDateTime(2008,6,14))
    >>> print(st)
    6 Trace(s) in Stream:
    XA.TOL0..LHE | 2008-06-12T23:59:59.640000Z - 2008-06-13T00:04:11.640000Z | 1.0 Hz, 253 samples
    XA.TOL0..LHE | 2008-06-13T00:04:12.640000Z - 2008-06-13T23:59:59.640000Z | 1.0 Hz, 86148 samples
    XA.TOL0..LHN | 2008-06-12T23:59:59.640000Z - 2008-06-13T00:04:11.640000Z | 1.0 Hz, 253 samples
    XA.TOL0..LHN | 2008-06-13T00:04:12.640000Z - 2008-06-13T23:59:59.640000Z | 1.0 Hz, 86148 samples
    XA.TOL0..LHZ | 2008-06-12T23:59:59.640000Z - 2008-06-13T00:04:21.640000Z | 1.0 Hz, 263 samples
    XA.TOL0..LHZ | 2008-06-13T00:04:22.640000Z - 2008-06-13T23:59:59.640000Z | 1.0 Hz, 86138 samples
    
    Also adds a Dbrecord as an attribute of the Trace
    
    >>> st[0].db
    Dbrecord('View43' -> TOL0 LHE 1213229044.64::1213315451.64)
 
    """
    from obspy.core import read, Stream, UTCDateTime

    if isinstance(database,Dbptr):
        db = Dbptr(database)
        db = db.lookup(table='wfdisc')
    else:
        raise TypeError("Must input a string or pointer to a valid database")
        
    if station is not None:
        db = db.subset('sta=~/{0}/'.format(station))
    if channel is not None:
        db = db.subset('chan=~/{0}/'.format(channel))
    if starttime is not None and endtime is not None:
        ts = starttime.timestamp
        te = endtime.timestamp
        db = db.subset('endtime > {0} && time < {1}'.format(ts,te) )
    else:
        ts = starttime
        te = endtime
    assert db.nrecs() is not 0, "No records for given time period"
    
    st = Stream()
    for db.record in range(db.nrecs() ):
        fname = db.filename() 
        dbr = RowPointerDict(db)
        t0 = UTCDateTime(dbr['time'])
        t1 = UTCDateTime(dbr['endtime'])
        if dbr['time'] < ts:
            t0 = starttime
        if dbr['endtime'] > te:
            t1 = endtime
        if os.path.exists(fname):
            _st = read(fname, starttime=t0, endtime=t1)      # add format?
            _st = _st.select(station=dbr['sta'], channel=dbr['chan']) #not location aware
            #_st[0].db = dbr
            if dbr['calib'] < 0:
                _st[0].data *= -1
            st += _st
    # Close what we opened, BUT garbage collection may take care of this:
    # if you have an open pointer but pass db name as a string, global
    # use of your pointer won't work if this is uncommented:
    #
    #if isinstance(database,str):
    #   db.close()
    return st

def dbloc_source_db(db, pointer=True):
    """
    Checks if you are in a dbloc2 'trial' db and returns the source
    one if you are, otherwise returns the same Dbptr. This is for running
    interactive scripts lauched from dbloc2 and writing to a non-volitile
    original db.
    
    INPUT: Dbptr of current temp database in dbloc2
    OUTPUT: Dbptr to database that dbloc2 is using.
    """
    try:
        from antelope.stock import pfget
    except ImportError:
        from antelope.stock import pfread as pfget
    
    db = Dbptr(db, perm='r+') 
    dbname = db.query('dbDATABASE_NAME')
    pf_settings = pfget('dbloc2')
    pfdef = pf_settings['Define']
    tempdb = pfdef['Temporary_db']
    workdir = pfdef['Work_dir']
    dblocdb = os.path.join(workdir,tempdb)
    if dbname.endswith(tempdb):
        # path of trial db from dbloc2
        dbcwd = os.path.dirname(dbname)
        # relative name of 1st db in 'trial' database decriptor file
        dbpath0 = db.query('dbDBPATH').split(':')[0].translate(None,'{}')
        # full absolute path database name to source
        dbname = os.path.abspath(os.path.join(dbcwd, dbpath0))
        db.close()
        db = Dbptr(dbname, perm='r+')
    if pointer:
        return db
    else:
        db.close()
        return dbname

def eventfocalmech2db(event=None, database=None):
    """
    Write the preferred HASH solution to Datascope database.
    
    Writes to 'fplane', 'predmech' and 'predarr' tables
    """
    focm = event.preferred_focal_mechanism()
    o = focm.triggering_origin_id.getReferredObject()
    
    plane1 = focm.nodal_planes.nodal_plane_1
    plane2 = focm.nodal_planes.nodal_plane_2
    T = focm.principal_axes.t_axis
    P = focm.principal_axes.p_axis
    orid = int(o.creation_info.version)
    
    db = Dbptr(database, perm='r+')
    try:
        # Use the original db if in a dbloc2 'tmp/trial' db
        #db = dbloc_source_db(db)
        # save solution as a new mechid
        mechid = db.nextid('mechid')
        # in fplane...
        dbfpln = dblookup(db,table='fplane')
        dbfpln.record = dbfpln.addnull()
        dbfpln.putv('orid', orid,
            'str1', round(plane1.strike,1) ,
            'dip1', round(plane1.dip,1) ,
            'rake1',round(plane1.rake,1),
            'str2', round(plane2.strike,1) ,
            'dip2', round(plane2.dip,1) ,
            'rake2',round(plane2.rake,1),
            'taxazm',round(T.azimuth,1),
            'taxplg',round(T.plunge,1),
            'paxazm',round(P.azimuth,1),
            'paxplg',round(P.plunge,1),
            'algorithm', focm.method_id.resource_id,
            'auth', focm.creation_info.author,
            'mechid', mechid,
            )
        dbpmec = dblookup(db,table='predmech')
        dbparr = dblookup(db,table='predarr')
        for av in o.arrivals:
            pk = av.pick_id.getReferredObject()
            if pk.polarity is 'positive':
                fm = 'U'
            elif pk.polarity is 'negative':
                fm = 'D'
            else:
                continue
            
            arid = int(av.creation_info.version)
            
            # ..and predmech
            dbpmec.record = dbpmec.addnull()
            dbpmec.putv('arid', arid,
                        'orid', orid,
                        'mechid', mechid,
                        'fm', fm,
                        )
            # if there are entries for this arrival already, write over it...
            dbparr.record = dbparr.find('arid=={0} && orid=={1}'.format(arid, orid))
            if dbparr.record < 0:
                dbparr.record = dbparr.addnull()
            dbparr.putv('arid', arid,
                        'orid', orid, 
                        'esaz', av.azimuth, 
                        'dip' , av.takeoff_angle,
                        )
    except Exception as e:
        raise e
    finally:
        db.close()

def get_first_motions(dbname, orid=None):
    """
    Port of Gabe/Mark dbprocess for getting info to pass to an FM calulator
    
    Right now, gets origin, arrival info, and joins wfdisc for filenames to the waveform
    """
    db = Dbptr(dbname)
    db = dbprocess(db ,['dbopen origin', 'dbsubset orid=={0}'.format(orid),
            'dbjoin origerr', 'dbjoin assoc',  'dbjoin arrival', 'dbsubset iphase =~ /.*[Pp].*/',
            'dbsubset fm =~ /.*[UuCcDdRr.].*/',
            'dbjoin wfdisc', 'dbsubset chan==wfdisc.chan', 'dbsort arrival.time'])
            #'dbjoin -o affiliation', 'dbjoin -o site',
            #
            #
            #'dbsubset (ondate <= time)',
            #'dbsubset (time <= offdate) || (offdate == -1)']
            #)
    return db

