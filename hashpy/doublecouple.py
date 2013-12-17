# -*- coding: utf-8 -*-
"""
Calulator for getting double couple information

:copyright:
    Mark Williams (2013) 
    Nevada Seismological Laboratory

"""

import numpy as np

class NodalPlane(list):
    """
    List to hold strike, dip, rake of a nodal plane
    
    Overview
    --------
    Basically, a list consisting of:
    [strike, dip, rake]
    with each element accessible by name as well as index.
    
    Construct with sequence, list, or named keyword args, see
    constructor doc for details.
    
   
    :type  strike: int or float 
    :param strike: degrees of strike
    :type  dip:    int or float 
    :param dip:    degrees of dip
    :type  rake:   int or float 
    :param rake:   degrees of rake

    """
    
    @property
    def strike(self):
        return self[0]
    @strike.setter
    def strike(self, value):
        self[0] = value

    @property
    def dip(self):
        return self[1]
    @dip.setter
    def dip(self, value):
        self[1] = value

    @property
    def rake(self):
        return self[2]
    @rake.setter
    def rake(self, value):
        self[2] = value

    def __init__(self, *args, **kwargs):
        """
        Build as a list or use keywords
        
        :param args:   One list or three numbers (strike,dip,rake)
        :param kwargs: Optionally specify 'strike','dip', 'rake' by name

        .. rubric:: Constructor Forms

        NodalPlane(strk, dp, rk)
        NodalPlane([strk,dp,rk])
        NodalPlane(strike=strk, dip=dp, rake=rk)
        
        .. rubric:: Examples
        
        >>> l = [123, 45, 67]
        >>> p = NodalPlane(l)
        >>> p = NodalPlane(145, 45, 67)
        >>> p = NodalPlane(strike=145, dip=45, rake=67)
        >>> p.dip = 30
        """
        super(NodalPlane,self).__init__([None,None,None])
        
        if args:
            if isinstance(args[0], list) or isinstance(args[0], tuple) and len(args[0]) == 3 :
                self[:] = [float(n) for n in args[0]]
            elif len(args) == 3:
                self[:] = [float(n) for n in args]
            else:
                pass
        for key, value in kwargs.items():
            if hasattr(self, key):
                self.__setattr__(key, float(value))


class DoubleCouple(object):
    """
    Calulate nodal planes and P and T axes of a double couple focal mech
    
    The attributes are set up to calulate everything on the fly from the
    initial plane (strike, dip, rake), so you can change something (like
    a rake in your primary plane), and calling for a 'P' axis, e.g. will
    give you the new answer...
    
    :type  plane1:   :class:`~hashpy.doublecoule.NodalPlane`
    :param plane1:   Primary plane containing strike/dip/rake
    :type  plane2:   :class:`~hashpy.doublecoule.NodalPlane`
    :param plane2:   Auxiliary plane caculated from primary
    :type  axis:     dict of key/dict pairs
    :param axis:     Keys ('P' and 'T') contain dict with keys 
        'azimuth' and 'dip' for that axis
    
    .. rubric:: Example

    >>> dc = DoubleCouple([270, 90, 180])
    >>> dc.plane1
    >>> dc.plane2
    >>> dc.axis['P']

    """
    
    _plane = None
    
    @property
    def plane1(self):
        '''Return Preferred plane'''
        return NodalPlane(*self._plane)
    
    @property
    def plane2(self):
        '''Return Auxiliary plane'''
        auxplane = self.aux_plane(*self._plane)
        return NodalPlane(*auxplane)
    
    @property
    def axis(self):
        '''return direction and dip for P and T axes'''
        dipP, dipT, aziP, aziT = self.nodal2pt(*self.plane1+self.plane2)
        return {'P': {'azimuth': aziP, 'dip': dipP}, 'T': {'azimuth': aziT, 'dip': dipT}}
    
    def __init__(self, nodal_plane=None):
        self._plane = nodal_plane
    
    @staticmethod
    def zero_360(str1):
        '''Put an angle between 0 and 360 degrees
    
        Genevieve Patau
        '''
        if str1 >= 360:
            str1 -= 360
        elif str1 < 0:
            str1 += 360
        else:
            pass
        return str1
    
    @classmethod
    def nodal2pt(cls, str1,da1,sa1,str2,da2,sa2):
        '''Compute azimuth and plungement of P-T axis 
        (from nodal plane strikes, dips and rakes.)
        
        Mark's python port from Gabe's perl port from:
        FORTRAN routines of Anne Deschamps ::
        
        Inputs
        ------
        *args == (str1,da1,sa1,str2,da2,sa2)
        For each plane:
        str : strike angle in degrees
        da  : dip angle in degrees
        sa  : rake (slip angle) in degrees
        
        Returns
        -------
        Dips and azimuths of P and T axes
        (dip_p, dip_t, azi_p, azi_t)
        
        (Original fxn used azimuth of dip plane, not strike)
        str1 = dd1 - 90;
        str2 = dd2 - 90;
        '''
        # Constants, mostly unnecessary, fix later:
        # e.g. M_PI = np.pi
        EPSIL   = .0001
        M_PI    = 3.14159265358979323846
        M_PI_2  = 1.57079632679489661923
        M_PI_4  = 0.78539816339744830962
        M_SQRT2 = 1.41421356237309504880
        TWO_PI  = 6.28318530717958647692
        D2R  = M_PI / 180.0
        
        #my ($pp, $dp, $pt, $dt, $xp, $yp);
        
        im = 0
        pure_strike_slip = 0
        
        if abs(np.sin(sa1 * D2R)) > EPSIL:
            im = sa1 / abs(sa1)
        elif abs(np.sin(sa2 * D2R)) > EPSIL:
            im = sa2 / abs(sa2)
        else:
            pure_strike_slip = 1
        
        if pure_strike_slip:
            if np.cos(sa1 * D2R) < 0:
                pp = cls.zero_360(str1 + 45)
                pt = cls.zero_360(str1 - 45)
            else:
                pp = cls.zero_360(str1 - 45);
                pt = cls.zero_360(str1 + 45);
            dp = 0
            dt = 0
        else:
            cd1 = np.cos(da1 * D2R) *  M_SQRT2
            sd1 = np.sin(da1 * D2R) *  M_SQRT2
            cd2 = np.cos(da2 * D2R) *  M_SQRT2
            sd2 = np.sin(da2 * D2R) *  M_SQRT2
            cp1 = -(np.cos(str1 * D2R)) * sd1
            sp1 = np.sin(str1 * D2R) * sd1
            cp2 = -(np.cos(str2 * D2R)) * sd2
            sp2 = np.sin(str2 * D2R) * sd2
            
            amz = -(cd1 + cd2)
            amx = -(sp1 + sp2)
            amy =   cp1 + cp2
            dx  = np.arctan2(np.sqrt(amx * amx + amy * amy), amz) - M_PI_2
            px  = np.arctan2(amy, -amx)

            if px < 0:
                px += TWO_PI
            
            amz   = cd1 - cd2
            amx   = sp1 - sp2
            amy   = - cp1 + cp2
            dy = np.arctan2(np.sqrt(amx * amx + amy * amy), -abs(amz)) - M_PI_2
            py = np.arctan2(amy, -amx)

            if amz > 0:
                py -= M_PI
            if py < 0:
                py += TWO_PI

            if im == 1:
                dp = dy
                pp = py
                dt = dx
                pt = px
            else:
                dp = dx
                pp = px
                dt = dy
                pt = py

            pp *= 180 / M_PI
            dp *= 180 / M_PI
            pt *= 180 / M_PI
            dt *= 180 / M_PI
        
        # I added this line b/c the names are confusing - MCW
        dip_p, dip_t, azi_p, azi_t = dp, dt, pp, pt
        return dip_p, dip_t, azi_p, azi_t

    # These two taken from obspy.imaging.beachball
    @staticmethod
    def get_strike_dip(n, e, u):
        """
        Finds strike and dip of plane given normal vector having components n, e,
        and u.

        Adapted from MATLAB script
        `bb.m <http://www.ceri.memphis.edu/people/olboyd/Software/Software.html>`_
        written by Andy Michael and Oliver Boyd.
        """
        r2d = 180 / np.pi
        if u < 0:
            n = -n
            e = -e
            u = -u

        strike = np.arctan2(e, n) * r2d
        strike = strike - 90
        while strike >= 360:
                strike = strike - 360
        while strike < 0:
                strike = strike + 360
        x = np.sqrt(np.power(n, 2) + np.power(e, 2))
        dip = np.arctan2(x, u) * r2d
        return (strike, dip)

    @classmethod
    def aux_plane(cls, s1, d1, r1):
        """
        Get Strike and dip of second plane.

        Adapted from MATLAB script
        `bb.m <http://www.ceri.memphis.edu/people/olboyd/Software/Software.html>`_
        written by Andy Michael and Oliver Boyd.
        """
        r2d = 180 / np.pi

        z = (s1 + 90) / r2d
        z2 = d1 / r2d
        z3 = r1 / r2d
        # slick vector in plane 1
        sl1 = -np.cos(z3) * np.cos(z) - np.sin(z3) * np.sin(z) * np.cos(z2)
        sl2 = np.cos(z3) * np.sin(z) - np.sin(z3) * np.cos(z) * np.cos(z2)
        sl3 = np.sin(z3) * np.sin(z2)
        (strike, dip) = cls.get_strike_dip(sl2, sl1, sl3)

        n1 = np.sin(z) * np.sin(z2)  # normal vector to plane 1
        n2 = np.cos(z) * np.sin(z2)
        h1 = -sl2  # strike vector of plane 2
        h2 = sl1
        # note h3=0 always so we leave it out
        # n3 = np.cos(z2)

        z = h1 * n1 + h2 * n2
        z = z / np.sqrt(h1 * h1 + h2 * h2)
        z = np.arccos(z)
        rake = 0
        if sl3 > 0:
            rake = z * r2d
        if sl3 <= 0:
            rake = -z * r2d
        return (strike, dip, rake)


