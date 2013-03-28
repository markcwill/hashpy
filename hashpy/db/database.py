#
# database.py 
# by Mark
# 2013-2-13
#
# functions to map Antelope css3.0 to obspy Event class
# which can then write out QuakeML
#
# Requirements: 'antelope' Antelope python API
#               'obspy'    ObsPy (version with event, quakeml support)
#               'aug'      Antelope Users Group contributed python module
#
import settings # NSL specific module
from antelope.datascope import Dbptr, dbopen, dbprocess
from aug.contrib.orm import AttribDbptr, open_db_or_string
from numpy import array
from obspy.core.utcdatetime import UTCDateTime
from obspy.core.quakeml import Pickler
from obspy.core.util import tostring, gps2DistAzimuth
from obspy.core.event import (Catalog, Event, Origin, CreationInfo, Magnitude,
    EventDescription, OriginUncertainty, OriginQuality, CompositeTime,
    ConfidenceEllipsoid, StationMagnitude, Comment, WaveformStreamID, Pick,
    QuantityError, Arrival, FocalMechanism, MomentTensor, NodalPlanes,
    PrincipalAxes, Axis, NodalPlane, SourceTimeFunction, Tensor, DataUsed,
    ResourceIdentifier, StationMagnitudeContribution)


#############################################################################
# Tools for creating obspy event objects
#############################################################################

agency   = settings.AGENCY_CODE 
placedb  = settings.PLACE_DB
site_url = settings.URL
site_tag = settings.TAG

def quakeml_rid(obj, url=site_url, tag=site_tag):
    """
    Return a resource identifier for quakeml (for NSL)
    
    *** BASED ON NSL QuakeML CONVENTIONS! ***

    This takes a generic site tag and returns a ResourceIdentifier
    with a resource_id based on the name of the object and the 'version'
    in its creation_info

    Could also use creation_info.agency_id...
    """
    # Make publicID using links to NSL site, if possible.
    # else, just make up a generic one, like 'quakeml:nsl/Origin/908347'
    if isinstance(obj, Event):
        evid = obj.creation_info.version
        l = ['quakeml:'+ url + evid]
    else:
        prefix = 'quakeml:'+ tag
        name   = obj.__class__.__name__
        id_num = obj.creation_info.version
        l = [prefix, name, id_num]
    # In case of multiple magnitudes, make Mag unique with type
    if isinstance(obj, Magnitude):
        l.insert(2, obj.magnitude_type)
        
    ridstr = '/'.join(l)
    return ResourceIdentifier(ridstr)

def azimuth2compass(azimuth):
    """
    Return 1 of 8 compass directions from an azimuth in degrees from N
    """
    needle = None
    if azimuth < 22.5:
        needle = 'N'
    elif azimuth < 67.5:
        needle = 'NE'
    elif azimuth < 112.5:
        needle = 'E'
    elif azimuth < 157.5:
        needle = 'SE'
    elif azimuth < 202.5:
        needle = 'S'
    elif azimuth < 247.5:
        needle = 'SW'
    elif azimuth < 292.5:
        needle = 'W'
    elif azimuth < 337.5:
        needle = 'NW'
    else:
        needle = 'N'
    return needle


class DbConnection(object):
    '''
    Connction to an Antelope Datascope database
    
    Methods return AttribDbptrs which contain records from views
    needed to convert to ObsPy classes.
    '''
    ptr = None
    
    def __init__(self, database):
        """Open a db"""
        self.ptr, _opened = open_db_or_string(database)
    
    def close(self):
        self.ptr.close()
    
    def get_event_origins(self, evid=None, orid=None):
        """
        Gets the origin-origerr joined line for a given event
        
        Finds the preferred orid, joins origin and origerr, and
        subsets the view for that orid. Returns object pointer
        to all the fields in that view for that record.
        
        Input:  int of evid, orid
        Output: aug.contrib.orm.dbpointers.AttribDbptr of origins
        """
        d = Dbptr(self.ptr)
        if orid is not None:
            substr = 'dbsubset origin.orid=={0}'.format(orid)
        else:
            substr = 'dbsubset origin.evid=={0}'.format(evid)
        
        d = dbprocess(d, ['dbopen origin', 'dbjoin -o origerr', substr, 'dbsort lddate'])
        return AttribDbptr(d)
    
    def get_origin_arrivals(self, orid=None):
        """
        Gets the assoc-arrival joined lines for a given origin
    
        Joins assoc to arrival table, then subsets on orid
        """
        d = Dbptr(self.ptr)
        d = dbprocess(d, ['dbopen assoc', 'dbsubset orid=={0}'.format(orid), 'dbjoin arrival'])
        # took out 'dbjoin snetsta' and 'dbjoin schanloc'
        # schanloc wouldn't join right, and not all dbsnapshots have snetsta for some weird reason...
        return AttribDbptr(d)
    
    def get_first_motions(self, orid=None):
        """
        Port of Gabe/Mark dbprocess for getting info to pass to an FM calulator
        
        (Could build this from the other two functions? OR take out the origin-origerr)
        """
        d = Dbptr(self.ptr)
        d = dbprocess(d,['dbopen origin', 'dbsubset orid=={0}'.format(orid),
                'dbjoin origerr', 'dbjoin assoc',  'dbjoin arrival', 'dbsubset iphase =~ /.*[Pp].*/',
                'dbsubset fm =~ /.*[UuCcDdRr.].*/',
                'dbjoin wfdisc', 'dbsubset chan==wfdisc.chan', 'dbsort arrival.time'])
                #'dbjoin -o affiliation', 'dbjoin -o site',
                #
                #
                #'dbsubset (ondate <= time)',
                #'dbsubset (time <= offdate) || (offdate == -1)']
                #)
        return AttribDbptr(d)


class Converter(object):
    '''
    Converter functions for Antelope --> ObsPy
    
    Methods/Functions use Object Relational Mapper objects as inputs
    
    Inputs usually reference databases
    Outputs are ObsPy event classes instances or lists of instances
    '''
    
    def _dbmapper(self, keymap, source={}, dest={}):
        """Maps from AttribDbptr object to an AttribDict (or any dict)
    
        Has support for Antelope NULLs (inserts a None)
        Catches if field doesn't exist (inserts a None)
    
        Input
        -----
        keymap => dict where keys are keys of dest,
                  values are valid for source.get(value)
    
        source => AttribDbptr (or anything with a 'get' fxn)
    
        dest   => Any instance of a dict of values with new keys
    
        Output
        ------
        returns dest
        """
        for _k in keymap:
            try:
                dest[_k] = source.get(keymap[_k])
            except TypeError:
                dest[_k] = None
        return dest
    
    def _create_dict(self, db, field):
        """Make a dict of {field:value} only if field is not NULL"""
        value = db.get(field)
        if value:
            return { field : value }
        else:
            return None
    
    def ptr2arrivals(self, adbp):
        """
        Return lists of obspy Arrivals and Picks from a DbrecordPtrList (or AttribDbptr)
        to records in an assoc-arrival join
        
        Inputs
        ------
        adbp : AttribDbptr to an assoc-arrival join
        
        Returns
        -------
        picks    : list of obspy Pick objects referenced by arrivals
        arrivals :  list of obspy Arrival objects with refs to picks
        
        """
        picks = []
        arrivals = []
        for db in adbp:
            p = Pick()
            p.time = UTCDateTime(db['time'])
            p.waveform_id = WaveformStreamID(station_code=db.get('sta'), channel_code=db.get('chan'))
            p.horizontal_slowness = db.get('slow')
            p.horizontal_slowness_errors = self._create_dict(db, 'delslo')
            p.backazimuth = db.get('azimuth')
            p.backazimuth_errors = self._create_dict(db, 'delaz')
            
            on_qual = db['qual'].lower()
            if 'i' in on_qual:
                p.onset = "impulsive"
            elif 'e' in on_qual:
                p.onset = "emergent"
            elif 'w' in on_qual:
                p.onset = "questionable"
            else:
                pass # or set to None
            
            p.phase_hint = db.get('iphase')
            
            pol = db['fm'].lower()
            if 'c' in pol or 'u' in pol:
                p.polarity = "positive"
            elif 'd' in pol or 'r' in pol:
                p.polarity = "negative"
            else:
                p.polarity = "undecidable"
            
            p.evaluation_mode = "automatic"
            if 'orbassoc' not in db['auth']:
                p.evaluation_mode = "manual"
            
            p.evaluation_status = "preliminary"
            if p.evaluation_mode is "manual":
                p.evaluation_status = "reviewed"
            
            p.creation_info = CreationInfo(version=db['arid'], creation_time=UTCDateTime(db['arrival.lddate']), agency_id=agency, author=db['auth'])
            p.resource_id = quakeml_rid(p)
            picks.append(p)
    
            # Now do the arrival
            a = Arrival()
            a.pick_id = ResourceIdentifier(p.resource_id.resource_id, referred_object=p)
            a.phase = db.get('phase')
            a.azimuth = db.get('esaz')
            a.distance = db.get('delta')
            a.takeoff_angle = db.get('ema')
            a.takeoff_angle_errors = self._create_dict(db, 'emares')
            a.time_residual = db.get('timeres')
            a.horizontal_slowness_residual = db.get('slores')
            a.time_weight = db.get('wgt')
            # TEMP HACK TO SET id
            a.earth_model_id = ResourceIdentifier('quakeml:' + site_tag +'/VelocityModel/'+ db['vmodel'])
            a.creation_info = CreationInfo(version=db['arid'], creation_time=UTCDateTime(db['assoc.lddate']))
            a.resource_id = quakeml_rid(a)
            arrivals.append(a)
        return picks, arrivals
    
    def record2maglist(self, db, mtypes=['ml','mb']):
        """Return list of obspy event Magnitudes from a DbrecordPtr to an origin table"""
        mags = []
        for mtype in mtypes:
            if db.get(mtype):
                mag_n = Magnitude(
                    mag=db[mtype],
                    magnitude_type=mtype, 
                    creation_info=CreationInfo(creation_time=UTCDateTime(db.lddate), agency_id=agency.upper(), version=str(db.orid)),
                    )
                mag_n.resource_id = quakeml_rid(mag_n)
                mags.append(mag_n)
        return mags
    
    def record2origin(self, db):
        """Return obspy event Origin from DbrecordPtr of origin/origerr view"""
        
        # Maps for mapper funtion: obspy <- antelope
        origin_map = {
            'latitude'  : 'origin.lat',
            'longitude' : 'origin.lon',
            'depth'     : 'origin.depth',
        }
    
        origin_quality_map = {
            'associated_phase_count': 'origin.nass',
            'used_phase_count'      : 'origin.ndef',
            'standard_error'        : 'origerr.sdobs',
        }
    
        confidence_ellipsoid_map = {
            'semi_minor_axis_length'         : 'origerr.sminax',
            'semi_major_axis_length'         : 'origerr.smajax',
            'semi_intermediate_axis_length'  : 'origerr.sdepth',
            'major_axis_azimuth'             : 'origerr.strike',
        }
        
        ellipse = ConfidenceEllipsoid(major_axis_plunge=0, major_axis_rotation=90)
        ellipse = self._dbmapper(confidence_ellipsoid_map, db, ellipse)
        # Now construct Origin object
        origin               = self._dbmapper(origin_map, db, Origin() )
        origin.time          = UTCDateTime(db['origin.time'])
        origin.quality       = self._dbmapper(origin_quality_map, db, OriginQuality() )
        origin.creation_info = CreationInfo(creation_time=UTCDateTime(db.lddate), agency_id=agency.upper(), version=str(db.orid))
        origin.origin_uncertainty = OriginUncertainty(confidence_ellipsoid=ellipse) 
        origin.resource_id   = quakeml_rid(origin)
    
        if 'orbassoc' in db.auth:
            origin.evaluation_mode = "automatic"
            origin.evaluation_status = "preliminary"
        else:
            origin.evaluation_mode = "manual"
            origin.evaluation_status = "reviewed"
        return origin
        
    def coords2evDescription(self, database, latitude, longitude):
        """
        Get the nearest place to a lat/lon from a db with a 'places' table
        """
        d = dbopen(database).lookup(table='places')
        adbp = AttribDbptr(d)
        stats = array([gps2DistAzimuth(latitude, longitude, r.lat, r.lon) for r in adbp])
        ind = stats.argmin(0)[0]
        minstats = stats[ind]
        minrec = adbp[int(ind)]
        dist, azi, backazi = minstats
        compass = azimuth2compass(backazi)
        place_info = {'distance': dist/1000., 'direction': compass, 'city': minrec.place, 'state': minrec.state}
        d.close()
        nearest_city_string = "{distance:0.1f} km {direction} of {city}, {state}".format(**place_info)
        return EventDescription(nearest_city_string, "nearest cities")


def db2event(database, evid=None, orid=None, phase_data=False):
    """
    Creates an ObsPy Event object from an Antelope Datascope database
    
    Currently, passing just an evid is handled by getting all origins
    sorted by creation date, and taking the most recent.
    
    Inputs
    ------
    database   : str of database name OR antelope.Dbptr to open database
    evid       : int of CSS3.0 Event ID
    orid       : int of CSS3.0 Origin ID
    phase_data : bool of whether to include associated picks (False)
    
    Returns
    -------
    obspy.core.event.Event object instance
    """
    event = Event()
    dbc = DbConnection(database)
    adbp = dbc.get_event_origins(evid, orid)
    if len(adbp) > 0:
        # for now, take the most recent origin (sorted by mod time)
        db = adbp[-1]
        # build Origin and list of Magnitude objects
        origin  = Converter().record2origin(db)
        maglist = Converter().record2maglist(db)
        ev_desc = Converter().coords2evDescription(placedb, db.lat, db.lon)
        if phase_data:
            adbp1 = dbc.get_origin_arrivals(db.orid)
            picks, arrivals = Converter().ptr2arrivals(adbp1)
            origin.arrivals = arrivals
            event.picks = picks
        event.origins = [origin]
        event.preferred_origin_id = origin.resource_id.resource_id
        event.event_descriptions = [ev_desc]
        # If magnitude is null or negative, maglist will be empty.
        if maglist:
            event.magnitudes = maglist
            event.preferred_magnitude_id = maglist[0].resource_id.resource_id
        event.creation_info = origin.creation_info.copy()
        event.creation_info.version = db.evid
        event.resource_id = quakeml_rid(event)
    else:
        # Punt for now if no valid records in db
        pass
    dbc.close()
    return event
