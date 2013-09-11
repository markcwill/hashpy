# -*- coding: utf-8 -*-
"""
 focalmechplotter.py

 by Mark Williams 2012.313

 Focal Mechanism Plotter using ObsPy Event as I/O
 (See docstring for class 'FocalMechPlotter')

"""

from matplotlib import pyplot as plt
from matplotlib.gridspec import GridSpec
import mplstereonet

class FocalMechPlotter(object):
    """
    Class to create stereonet plot of a HASH first motion
    
    This is a simple, slow, matplotlib script to plot a first motion
    focal mechanism using ObsPy object input. Currently supports:

    Event with a FocalMechanism (marked as preferred)
                 Origin/Arrivals (marked as preferred)
                 Picks (referred to by each Arrival)
    
    Notes
    -----
    Requires the mplstereonet matplotlib extension.

    """
    fig = None      # handle to figure
    ax = None       # list of axes in the figure
    h = None        # list of handles of FM picks 
    gs = None       # GridSpec instance of plot figure
    h_text = None   # handle for text labels
    event = None    # Event object
    ind = None      # indicies of which picks in Origin.arrivals are plotted
    _arrv = None    # current arrival
    _curh = None    # current handle
    _axis = None    # dict of axes for quick ref
    save = None     # function handle to a 'save' function
    source = None   # temp string of where to save
    
    picks_list = ('undecidable', 'positive', 'negative')
    wf_color = { 'positive' : 'red', 'negative' : 'blue', 'undecidable' : 'black' }
    
    @property       # preferred origin
    def _orig(self):
        '''Stub'''
        return self.event.preferred_origin()
    @property       # preferred focal mechanism
    def _focm(self):
        '''Stub'''
        return self.event.preferred_focal_mechanism()
    @property       # current pick
    def _pick(self):
        '''Pick pointed to by current arrival'''
        return self._arrv.pick_id.getReferredObject()
        
    def plot_on_stereonet(self, axis=None):
        '''Plot first motions from an Event instance
        
        Stored in self.event
        '''
        if axis:
            ax = axis
        else:
            ax = self._axis['stereonet']
        
        ax.clear() 
        ax.set_title('Origin: {0}'.format(self._orig.creation_info.version))
        tlab  = ax.set_azimuth_ticklabels([])
        
        h = []
        index = []
       
        # Planes
        np1 = self._focm.nodal_planes.nodal_plane_1
        np2 = self._focm.nodal_planes.nodal_plane_2
        strike1,dip1,rake1 = np1.strike, np1.dip, np1.rake
        strike2,dip2,rake2 = np2.strike, np2.dip, np2.rake
        
        # Plotting --------------------------------------#
        # plot best fit average plane, aux plane
        h_rk = ax.plane(strike1, dip1, color='black', linewidth=3)
        h_rk = ax.rake( strike1, dip1, -rake1, 'k^', markersize=8)
        h_rk = ax.plane(strike2, dip2, color='black', linewidth=3)
        # plot station takeoffs
        plot_specs = { 'picker'          : 5,
                       'markersize'      : 8,
                       'markeredgewidth' : 2,
                       }
        for ind, a in enumerate(self._orig.arrivals):
            #--- HASH takeoffs are 0-180 from vertical UP!!
            #--- Stereonet angles 0-90 inward (dip)
            #--- Classic FM's are toa from center???
            p = a.pick_id.getReferredObject()
            azi = a.azimuth - 90.
            toa = abs(a.takeoff_angle - 90)
            if p.polarity is 'positive':
                #plot_specs.update({'markeredgecolor' : 'black', 'markerfacecolor' : 'red'   })
                h += ax.rake(azi, toa, 90, 'o', markeredgecolor='black', markerfacecolor='red', **plot_specs)
            if p.polarity is 'negative':
                #plot_specs.update({'markeredgecolor' : 'blue', 'markerfacecolor' : 'white' })
                h += ax.rake(azi, toa, 90, 'o', markeredgecolor='blue', markerfacecolor='white', **plot_specs)
            index.append(ind)
            if True:
                h_text = ax.rake(azi, toa+5, 90, marker='$   {0}$'.format(p.waveform_id.station_code), color='black',markersize=20)
        
        for comm in self._focm.comments:
            if 'quality' in comm.resource_id.resource_id:
                qual = comm.text
            else:
                qual = None
        plane_str = "STRIKE1:{0: > 7.1f}\nDIP1:{1: > 7.1f}\nRAKE1:{2: > 7.1f}\n\nSTRIKE2:{3: > 7.1f}\nDIP2:{4: > 7.1f}\nRAKE2:{5: > 7.1f}"
        h_text = self.fig.text(0.25, 0.88, plane_str.format(strike1, dip1, rake1, strike2, dip2, rake2), ha='right', va='top', family='monospace')
        h_text = self.fig.text(0.25, 0.11, 'Quality:  {0}\n# of picks: {1}'.format(qual, len(self._orig.arrivals)), ha='right', va='top', family='monospace')  
        
        if not axis:
            self.h = h
            self.ind = index
            
    def __init__(self, event=None, save=None, source=None):
        """
        Create a plot for focal mechanisms

        Input
        -----
        event : obspy.core.event.Event containing 
            FocalMechanism
            Picks
            Origin/Arrivals

        """
        self.event = event
        self.ax = []
        self._axis = {}
        self.save = save
        self.source = source
        
        # Draw figure and set up 
        self.fig = plt.figure(facecolor='#D9D9EE')
        self.fig.canvas.set_window_title('Focal Mechanism')
        
        ### INTERACTIVE - MAIN MODE ###
        self.gs = GridSpec(8,5) # 8x5 grid of axis space
        
        # Stereonet axis
        ax = self.fig.add_subplot(self.gs[:,1:], projection='stereonet') # gs[:-2,:-1]
        self.ax.append(ax)
        self._axis['stereonet'] = ax
        
        # Plot first motions
        self.plot_on_stereonet()
        
        # Save and draw
        plt.show()
    
