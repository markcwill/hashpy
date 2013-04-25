#
# HashPype class which uses ObsPy Event as I/O for focal mech data
# 
from hashpy.hashpype import HashPype, HashError
from hashpy.doublecouple import DoubleCouple
from obspy.core.utcdatetime import UTCDateTime
from obspy.core.event import (Catalog, Event, Origin, CreationInfo, Magnitude,
    EventDescription, OriginUncertainty, OriginQuality, CompositeTime,
    ConfidenceEllipsoid, StationMagnitude, Comment, WaveformStreamID, Pick,
    QuantityError, Arrival, FocalMechanism, MomentTensor, NodalPlanes,
    PrincipalAxes, Axis, NodalPlane, SourceTimeFunction, Tensor, DataUsed,
    ResourceIdentifier, StationMagnitudeContribution)


class EventHashPype(HashPype):
    """
    Hashpype with Obspy Event I/O
    
    This is a class which inherits HashPype (methods to set up and run HASH)
    and defines input and output to HASH using ObsPy Event instances.
    
    Methods
    -------
    load_event(event):
     - Loads Event data into internal HashPype parameters
    
    output(event=None, only_fm_picks=False):
     - Adds HASH data (FocalMechanism, pick takeoff angles, etc) to
       new or existing Event
       
    run():
     - Convenience function that runs several HashPype methods in a row
       to set up and calculate a focal mechanism
    """
    p_index = None   # index of which picks were used in HASH
                     # len(p_index) == self.npol
    
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
    
    def load_event(self, event):
        """
        Load Event into HASH
        
        Takes origin, arrival, pick information from an Event and loads
        up the HASH arrays needed for the focal mechanism calculation.
        
        Input
        -----
        event : obspy.core.event.Event instance
        
         ** event.preferred_origin_id should be set to the origin in
            event.origins you want to use
        """
        # Takes an obspy event and loads FM data into HASH
        _o = event.preferred_origin()
        _m = event.preferred_magnitude()
        
        self.tstamp = _o.time.timestamp
        self.qlat   = _o.latitude
        self.qlon   = _o.longitude
        self.qdep   = _o.depth
        self.icusp  = _o.creation_info.version
        self.seh    = _o.origin_uncertainty.confidence_ellipsoid.semi_major_axis_length
        self.sez    = _o.origin_uncertainty.confidence_ellipsoid.semi_intermediate_axis_length
        if _m:
            self.qmag   = _m.mag
        
        # The index 'k' is deliberately non-Pythonic to deal with the fortran
        # subroutines which need to be called and the structure of the original HASH code.
        # May be able to update with a rewrite... YMMV
        self.p_index = []
        k = 0
        for _i, arrv in enumerate(_o.arrivals):
            # load up params
            pick = arrv.pick_id.getReferredObject()
            self.sname[k]     = pick.waveform_id.station_code
            self.snet[k]      = pick.waveform_id.network_code
            self.scomp[k]     = pick.waveform_id.channel_code
            self.arid[k]      = pick.creation_info.version
            
            self.qazi[k] = arrv.azimuth
            self.dist[k] = arrv.distance * 111.2
            
            if (self.qazi[k] < 0.):
                self.qazi[k] += 360.
            
            if (self.dist[k] > self.delmax):
                continue
                
            if arrv.phase not in 'Pp':
                continue
            
            if (pick.polarity is 'positive'):
                self.p_pol[k] = 1
            elif (pick.polarity is 'negative'):
                self.p_pol[k] = -1
            else:
                continue
            
            if  (pick.onset is 'impulsive'):
                self.p_qual[k] = 0
            elif (pick.onset is 'emergent'):
                self.p_qual[k] = 1
            elif (pick.onset is 'questionable'):
                self.p_qual[k] = 1
            else:
                self.p_qual[k] = 0
                
            # polarity check in original code... doesn't work here
            #self.p_pol[k] = self.p_pol[k] * self.spol
            self.p_index.append(_i) # indicies of [arrivals] which passed
            k += 1
        self.npol = k # k is zero indexed in THIS loop
    
    def output(self, event=None, only_fm_picks=False):
        """
        Make an Event which includes the current focal mechanism information from HASH
        
        Use the 'only_fm_picks' flag to only include the picks HASH used for the FocalMechanism.
        This flag will replace the 'picks' and 'arrivals' lists of existing events with new ones.
        
        Inputs
        -------
        event : obspy.core.event.Event
        
        only_fm_picks : bool of whether to overwrite the picks/arrivals lists
        
        
        Returns
        -------
        obspy.core.event.Event
        
        Event will be new if no event was input, FocalMech added to existing event
        """
        # Returns new (or updates existing) Event with HASH solution
        n = self.npol
        if event is None:
            event = Event(focal_mechanisms=[], picks=[], origins=[])
            origin = Origin(arrivals=[])
            origin.time = UTCDateTime(self.tstamp)
            origin.latitude = self.qlat
            origin.longitude = self.qlon
            origin.depth = self.qdep
            origin.creation_info = CreationInfo(version=self.icusp)
            origin.resource_id = ResourceIdentifier('smi:hash/Origin/{0}'.format(self.icusp))
            for _i in range(n):
                p = Pick()
                p.creation_info = CreationInfo(version=self.arid[_i])
                p.resource_id = ResourceIdentifier('smi:nsl/Pick/{0}'.format(p.creation_info.version))
                p.waveform_id = WaveformStreamID(network_code=self.snet[_i], station_code=self.sname[_i], channel_code=self.scomp[_i])
                if self.p_pol[_i] > 0:
                    p.polarity = 'positive'
                else:
                    p.polarity = 'negative'
                a = Arrival()
                a.creation_info = CreationInfo(version=self.arid[_i])
                a.resource_id = ResourceIdentifier('smi:nsl/Arrival/{0}'.format(p.creation_info.version))
                a.azimuth = self.p_azi_mc[_i,0]
                a.takeoff_angle = self.p_the_mc[_i,0]
                a.pick_id = p.resource_id
                origin.arrivals.append(a)
                event.picks.append(p)
            event.origins.append(origin)
            event.preferred_origin_id = origin.resource_id.resource_id
        else: # just update the changes
            origin = event.preferred_origin()
            picks = []
            arrivals = []
            for _i in range(n):
                ind = self.p_index[_i]
                a = origin.arrivals[ind]
                p = a.pick_id.getReferredObject()
                a.takeoff_angle = self.p_the_mc[_i,0]
                picks.append(p)
                arrivals.append(a)
            if only_fm_picks:
                origin.arrivals = arrivals
                event.picks = picks
        # Use me double couple calculator and populate planes/axes etc
        s = self._best_quality_index
        # Right now, only FocalMechansim for 's' is used, for speed.
        # todo: loop through possible solutions, add them all, make 's' the preferred one
        #       so could use a future GUI to look at them all???
        dc = DoubleCouple([self.str_avg[s], self.dip_avg[s], self.rak_avg[s]])
        ax = dc.axis
        focal_mech = FocalMechanism()
        focal_mech.creation_info = CreationInfo(creation_time=UTCDateTime(), author=self.author)
        focal_mech.triggering_origin_id = origin.resource_id
        focal_mech.resource_id = ResourceIdentifier('smi:hash/FocalMechanism/{0}'.format(self.icusp))
        focal_mech.method_id = ResourceIdentifier('HASH')
        focal_mech.nodal_planes = NodalPlanes()
        focal_mech.nodal_planes.nodal_plane_1 = NodalPlane(*dc.plane1)
        focal_mech.nodal_planes.nodal_plane_2 = NodalPlane(*dc.plane2)
        focal_mech.principal_axes = PrincipalAxes()
        focal_mech.principal_axes.t_axis = Axis(azimuth=ax['T']['azimuth'], plunge=ax['T']['dip'])
        focal_mech.principal_axes.p_axis = Axis(azimuth=ax['P']['azimuth'], plunge=ax['P']['dip'])
        focal_mech.station_polarity_count = n
        focal_mech.azimuthal_gap = self.magap
        focal_mech.misfit = self.mfrac[s]
        focal_mech.station_distribution_ratio = self.stdr[s]
        focal_mech.comments.append(Comment(self.qual[s], resource_id=ResourceIdentifier('comment/quality')))
        #----------------------------------------
        event.focal_mechanisms.append(focal_mech)
        event.preferred_focal_mechanism_id = focal_mech.resource_id.resource_id
        return event
    
    def __str__(self):
        """
        One line statement for print function
        """
        x = self._best_quality_index
        s,d,r = self.str_avg[x], self.dip_avg[x], self.rak_avg[x]
        return 'Solution:{orid} |  STRIKE: {st:0.1f}  DIP: {dp:0.1f}  RAKE: {rk:0.1f}  | Quality:{q}'.format(orid=self.icusp,
            st=float(s), dp=float(d), rk=float(r), q=self.qual[x])
    
    def run(self, min_check=True, gap_check=True):
        """
        Script that runs all the methods for a run
        """
        self.load_velocity_models()
        self.generate_trial_data()
        self.calculate_takeoff_angles()
        check1 = self.check_minimum_polarity()
        check2 = self.check_maximum_gap()
        if min_check and not check1:
            raise HashError("Low number of picks = {0} / {1}".format(self.npol, self.npolmin))
        if gap_check and not check2:
            raise HashError("Exceeded maximum gap = azi:{0}/{1} plunge:{2}/{3}".format(self.magap, self.max_agap, self.mpgap, self.max_pgap))
        self.calculate_hash_focalmech()
        self.calculate_quality()
