#
#
# Replace * with functions
from hashpy.db import DbHashPype
from db.utils import DoubleCouple
from obspy.core.event import *

class HashRun(DbHashPype):
    '''
    STUB - Hashpype with Obspy Event I/O
    '''
    def input(self, event):
        '''stub'''
        # Takes an obspy event and loads FM data into HASH
        pass
    
    def output(self, event=None):
        '''stub'''
        # Returns new (or updates existing) Event with HASH solution
        if event is None:
            event = Event(focal_mechanisms=[], picks=[], origins=[])
        origin = Origin(arrivals=[])
        origin.latitude = self.qlat
        origin.longitude = self.qlon
        origin.creation_info = CreationInfo(version=self.icusp)
        origin.resource_id = ResourceIdentifier('smi:hash/Origin/{0}'.format(self.icusp))
        focal_mech = FocalMechanism()
        focal_mech.resource_id = ResourceIdentifier('smi:hash/FocalMechanism/{0}'.format(self.icusp))
        n = self.npol
        for _i in range(n):
            p = Pick()
            p.creation_info = CreationInfo(version=self.arid[_i])
            p.resource_id = ResourceIdentifier('smi:nsl/Pick/{0}'.format(p.creation_info.version))
            p.waveform_id = WaveformStreamID(station_code=self.sname[_i], channel_code=self.scomp[_i])
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
        # Use double couple calculator and populate planes/axes etc
        dc = DoubleCouple([self.str_avg[0], self.dip_avg[0], self.rak_avg[0]])
        ax = dc.axis
        focal_mech.nodal_planes = NodalPlanes()
        focal_mech.nodal_planes.nodal_plane_1 = NodalPlane(*dc.plane1)
        focal_mech.nodal_planes.nodal_plane_2 = NodalPlane(*dc.plane2)
        focal_mech.principal_axes = PrincipalAxes()
        focal_mech.principal_axes.t_axis = Axis(azimuth=ax.T.azi, plunge=ax.T.dip)
        focal_mech.principal_axes.p_axis = Axis(azimuth=ax.P.azi, plunge=ax.P.dip)
        # and various error info to be added here
        # todo: Find arrivals of existing shit brute force?
        event.focal_mechanisms.append(focal_mech)
        event.origins.append(origin)
        event.preferred_origin_id = origin.resource_id.resource_id
        event.preferred_focal_mechanism_id = focal_mech.resource_id.resource_id
        return event
        
class FocalMechEvent(Event):
    '''Stub'''
    pass

hro = HashRun()
hro.load_pf()
hro.get_phases_from_db('reno', orid=979567)
hro.load_velocity_models()
hro.generate_trial_data()
hro.calculate_takeoff_angles()
hro.calculate_hash_focalmech()
