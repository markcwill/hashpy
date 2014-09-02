# -*- coding: utf-8 -*-
"""
hashpy.io.antelopeIO

by Mark Williams 2012.313

Class to run HASHpy using Antelope database I/O

"""
import os
import math

import curds2.dbapi2 as dbapi2
from curds2.cursors import InteractiveCursor

from antelope import stock

# for HASH compatiblity, change later.
degrad = 180. / math.pi
rad = 1. / degrad


# Ported from nsl.commons
def get_pf(pfname):
    """Return a dict from a pf file"""
    if hasattr(stock, 'pfread'):
        return stock.pfread(pfname).pf2dict()
    elif hasattr(stock, 'pfget'):
        return stock.pfget(pfname)
    else:
        raise AttributeError("No pf function available")


def input(hp, dbname, evid=None, orid=None):
    """
    Input HASH data from Antelope database
    
    This will accept a database name OR Antelope Dbptr, and either
    an ORID, or an EVID, in which case the 'prefor' ORID is looked
    up and used.
    
    Inputs
    ------
    dbname  :   str or antelope.datascope.Dbptr
    orid    :   int of ORID
    evid    :   int of EVID
    """
    
    with dbapi2.connect(dbname) as conn:
        conn.cursor_factory = InteractiveCursor
        curs = conn.cursor()

        if orid is None:
            n = curs.execute.process(['dbopen event', 'dbsubset evid == '+str(evid)])
            orid = curs.fetchone()['prefor']
        
        n = curs.execute.process([
            'dbopen origin', 'dbsubset orid == '+str(orid),
            'dbjoin origerr', 'dbjoin assoc',  'dbjoin arrival',
            'dbjoin affiliation', 'dbjoin site',
            'dbsubset iphase =~ /.*[Pp].*/',
            'dbsubset (ondate <= time)',
            'dbsubset (time <= offdate) || (offdate == -1)',
            ])
        
        if n <= 0:
            raise ValueError("No picks for this ORID: {0}".format(orid) )
        
        ph = curs.fetchone()
            
        hp.tstamp = ph['origin.time']
        hp.qlat = ph['origin.lat']
        hp.qlon = ph['origin.lon']
        hp.qdep = ph['origin.depth']
        hp.qmag = ph['origin.ml']
        hp.icusp = ph['origin.orid']
        hp.seh = ph['origerr.smajax']
        hp.sez = ph['origerr.sdepth']
        hp.nrecs = n
        
        aspect = math.cos(hp.qlat / degrad) # convert using python later.
        
        # The index 'k' is deliberately non-Pythonic to deal with the fortran
        # subroutines which need to be called and the structure of the original HASH code.
        # May be able to update with a rewrite... YMMV
        k = 0
        for ph in curs:
            # Extract pick data from the db 
            hp.sname[k] = ph['sta']
            hp.snet[k] = ph['net']
            hp.scomp[k] = ph['chan']
            hp.pickonset[k] = ph['qual'] #.strip('.')
            hp.pickpol[k] = ph['fm']
            hp.arid[k] = ph['arid']
            
            flat, flon, felv = ph['site.lat'],ph['site.lon'],ph['site.elev']
            hp.esaz[k] = ph['esaz']

            # Distance and Azimuth filtering
            dx = (flon - hp.qlon) * 111.2 * aspect
            dy = (flat - hp.qlat) * 111.2
            dist = math.sqrt(dx**2 + dy**2)
            qazi = 90. - math.atan2(dy,dx) * degrad
            
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
        args = tuple(dbout,) 
        kwargs = dict(perm='r+')
    else:
        args = tuple(':memory:',)
        kwargs = dict()
    
    conn = dbapi2.connect(*args, **kwargs)
    conn.cursor_factory = InteractiveCursor

    mechid = curs.execute.nextid('mechid')
    
    n = curs.execute.lookup(table='fplane')
    
    recnum = curs.execute.addnull()
    curs.scroll(recnum, 'absolute')
    row = curs.fetchone()
    row.update({
            'orid': hp.icusp,
            'str1': round(str1,1) ,
            'dip1': round(dip1,1) ,
            'rake1': round(rak1,1),
            'algorithm': "HASH",
            'mechid': mechid,
            'auth': 'hashpy:'+ hp.author,
            'str2': round(str2,1) ,
            'dip2': round(dip2,1) ,
            'rake2': round(rak2,1),
            'taxazm': round(axes['T']['azimuth'],1),
            'taxplg': round(axes['T']['dip'],1),
            'paxazm': round(axes['P']['azimuth'],1),
            'paxplg': round(axes['P']['dip'],1),
            })
    
    curs_m = conn.cursor()
    curs_a = conn.cursor()
    nm = curs_m.execute.lookup(table='predmech')
    na = curs_a.execute.lookup(table='predarr')
    for k in range(hp.npol):
        if hp.p_pol[k] > 0:
            fm = 'U'
        else:
            fm = 'D'
        recnum = curs_m.execute.addnull()
        curs_m.scroll(recnum, 'absolute')
        row = curs.fetchone()
        row.update({
            'arid': int(hp.arid[k]) ,
            'orid': hp.icusp,
            'mechid': mechid,
            'fm': fm,
            })
        recnum = curs_a.execute.addnull()
        curs_a.scroll(recnum, 'absolute')
        row = curs.fetchone()
        row.update({
            'arid': int(hp.arid[k]),
            'orid': hp.icusp, 
            'esaz': hp.qazi[k], 
            'dip' : hp.p_the_mc[k,0],
            })
    return conn


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
    # 1) use Mark's version agnostic 'pfgetter' function
    # 2) TODO: man up and use JSON/YAML for configs, get typing for free
    pf_settings = get_pf(pffile)    

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


##############################################################################
# Antelope <-> ObsPy utilities
##############################################################################
#
# These are not necessary to run Antelope with HASH, for interactive/GUI
# scripts and housekeeping utils. 
#
def dbloc_source_db(db, pointer=True):
    """
    Checks if you are in a dbloc2 'trial' db and returns the source
    one if you are, otherwise returns the same Dbptr. This is for running
    interactive scripts lauched from dbloc2 and writing to a non-volitile
    original db.
    
    INPUT: Dbptr of current temp database in dbloc2
    OUTPUT: Dbptr to database that dbloc2 is using.
    """
    pf_settings = get_pf('dbloc2')

    conn = dbapi2.connect(db, perm='r+')
    curs = conn.cursor()
    dbname = curs.execute.query('dbDATABASE_NAME')
    pfdef = pf_settings['Define']
    tempdb = pfdef['Temporary_db']
    workdir = pfdef['Work_dir']
    dblocdb = os.path.join(workdir,tempdb)
    if dbname.endswith(tempdb):
        # path of trial db from dbloc2
        dbcwd = os.path.dirname(dbname)
        # relative name of 1st db in 'trial' database decriptor file
        dbpath0 = curs.execute.query('dbDBPATH').split(':')[0].translate(None,'{}')
        # full absolute path database name to source
        dbname = os.path.abspath(os.path.join(dbcwd, dbpath0))
        conn.close()
        conn = dbapi2.connect(dbname, perm='r+')
    if pointer:
        return conn
    else:
        conn.close()
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
    
    conn = dbapi2.connect(database, perm='r+')
    try:
        # Use the original db if in a dbloc2 'tmp/trial' db
        #db = dbloc_source_db(db)
        # save solution as a new mechid
        mechid = curs.execute.nextid('mechid')
        if mechid < 0:
            raise IOError("Error writing to database, check permissions")
        # in fplane...
        n = curs.execute.lookup(table='fplane')
        
        recnum = curs.execute.addnull()
        curs.scroll(recnum, 'absolute')
        row = curs.fetchone()
        row.update({
            'orid': orid,
            'str1': round(plane1.strike,1) ,
            'dip1': round(plane1.dip,1) ,
            'rake1':round(plane1.rake,1),
            'str2': round(plane2.strike,1) ,
            'dip2': round(plane2.dip,1) ,
            'rake2': round(plane2.rake,1),
            'taxazm': round(T.azimuth,1),
            'taxplg': round(T.plunge,1),
            'paxazm': round(P.azimuth,1),
            'paxplg': round(P.plunge,1),
            'algorithm': "HASHpy",
            'auth': focm.creation_info.author,
            'mechid': mechid,
            })
        curs_m = conn.cursor()
        curs_a = conn.cursor()
        nm = curs_m.execute.lookup(table='predmech')
        na = curs_a.execute.lookup(table='predarr')
        for av in o.arrivals:
            pk = av.pick_id.getReferredObject()
            if pk.polarity is 'positive':
                fm = 'U'
            elif pk.polarity is 'negative':
                fm = 'D'
            else:
                continue
            
            arid = int(av.creation_info.version)
            
            recnum = curs_m.execute.addnull()
            curs_m.scroll(recnum, 'absolute')
            row = curs.fetchone()
            row.update({
                'arid': arid,
                'orid': orid,
                'mechid': mechid,
                'fm': fm,
                })
            # if there are entries for this arrival already, write over it...
            recnum = dbparr.find('arid=={0} && orid=={1}'.format(arid, orid))
            if dbparr.record < 0:
                recnum = curs_a.execute.addnull()
            curs_a.scroll(recnum, 'absolute')
            row = curs.fetchone()
            row.update({
                'arid': arid,
                'orid': orid, 
                'esaz': av.azimuth, 
                'dip' : av.takeoff_angle,
                })
    except Exception as e:
        raise e
    finally:
        conn.close()


def get_first_motions(dbname, orid=None):
    """
    Port of Gabe/Mark dbprocess for getting info to pass to an FM calulator
    
    Right now, gets origin, arrival info, and joins wfdisc for filenames to the waveform
    """
    conn = dbapi2.connect(dbname)
    curs = conn.cursor()
    n = curs.execute.process([
        'dbopen origin', 'dbsubset orid=={0}'.format(orid),
        'dbjoin origerr', 'dbjoin assoc',  'dbjoin arrival', 'dbsubset iphase =~ /.*[Pp].*/',
        'dbsubset fm =~ /.*[UuCcDdRr.].*/',
        'dbjoin wfdisc', 'dbsubset chan==wfdisc.chan', 'dbsort arrival.time'
        ])
        #'dbjoin -o affiliation', 'dbjoin -o site',
        #
        #
        #'dbsubset (ondate <= time)',
        #'dbsubset (time <= offdate) || (offdate == -1)']
        #)
    return curs

