# -*- coding: utf-8 -*-
"""
 focalmechplotter.py

 by Mark Williams 2012.313

 Focal Mechanism Plotter using ObsPy Event as I/O
 (See docstring for class 'FocalMechPlotter')

"""
import matplotlib
matplotlib.use('TkAgg')
from matplotlib import pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.backend_bases import NavigationToolbar2
from matplotlib.backends.backend_tkagg import NavigationToolbar2TkAgg
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
    focm = None     # Current FM plotted
    ind = None      # indicies of which picks in Origin.arrivals are plotted
    _arrv = None    # current arrival
    _curh = None    # current handle
    _axis = None    # dict of axes for quick ref
    save = None     # function handle to a 'save' function
    source = None   # temp string of where to save
    
    picks_list = ('undecidable', 'positive', 'negative')
    wf_color = { 'positive' : 'red', 'negative' : 'blue', 'undecidable' : 'black' }
    
    @property
    def _fm_index(self):
        return self.event.focal_mechanisms.index(self.focm)
    @property
    def _num_fms(self):
        return len(self.event.focal_mechanisms)

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
        
    def plot_on_stereonet(self, axis=None, fm=None):
        '''Plot first motions from an Event instance
        
        Stored in self.event
        '''

        if axis:
            ax = axis
        else:
            ax = self._axis['stereonet']
        
        if fm is None:
            self.focm = self._focm
        else:
            self.focm = self.event.focal_mechanisms[fm]
        
        ax.clear() 
        ax.set_title('Origin: {0}'.format(self._orig.creation_info.version))
        tlab  = ax.set_azimuth_ticklabels([])
        
        h = []
        index = []
       
        # Planes
        np1 = self.focm.nodal_planes.nodal_plane_1
        np2 = self.focm.nodal_planes.nodal_plane_2
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
            p = a.pick_id.getReferredObject()
            # Calculate strike azi from direct (dip-pointing) azi 
            azi = a.azimuth - 90.
            #--- HASH takeoffs are 0-180 from vertical UP!!
            #--- obspy QuakeML are 0-180 from vertical DOWN
            #--- Stereonet angles 0-90 inward (dip)
            #--- Classic FM's are toa from center???
            if 0. <= a.takeoff_angle < 90.:
                toa = 90. - a.takeoff_angle  # complement for downward angles
            elif 90. <= a.takeoff_angle <= 180.:
                toa = 270. - a.takeoff_angle  # project upward angles
            else:
                raise ValueError("Takeoff angle ({0}) must be in [0, 180]".format(a.azimuth))
            
            if p.polarity is 'positive':
                #plot_specs.update({'markeredgecolor' : 'black', 'markerfacecolor' : 'red'   })
                h += ax.rake(azi, toa, 90, 'o', markeredgecolor='black', markerfacecolor='red', **plot_specs)
            if p.polarity is 'negative':
                #plot_specs.update({'markeredgecolor' : 'blue', 'markerfacecolor' : 'white' })
                h += ax.rake(azi, toa, 90, 'o', markeredgecolor='blue', markerfacecolor='white', **plot_specs)
            index.append(ind)
            if True:
                h_text = ax.rake(azi, toa+5, 90, marker='$   {0}$'.format(p.waveform_id.station_code), color='black',markersize=20)
        
        for comm in self.focm.comments:
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
        
        self._set_window_title()

    def draw_stereonet_axis(self, gridspec_slice=None):
        # Stereonet axis
        if gridspec_slice is None:
            gridspec_slice = self.gs_slice
        ax = self.fig.add_subplot(gridspec_slice, projection='stereonet') # gs[:-2,:-1]
        self.ax.append(ax)
        self._axis['stereonet'] = ax

    def _set_window_title(self):
        self.fig.canvas.set_window_title('Focal Mechanism: {0} of {1}'.format(self._fm_index+1, self._num_fms))
    
    def plot(self, solution=None):
        self.fig.clear()
        self.draw_stereonet_axis()
        self.plot_on_stereonet(fm=solution)
        plt.draw()

    def __init__(self, event=None, save=None):
        """
        Create a plot for focal mechanisms

        Input
        -----
        event : obspy.core.event.Event containing 
            FocalMechanism
            Picks
            Origin/Arrivals

        save : function handle, the instance is passed as 1st arg
               "save(self)", called on save button press

        """
        self.event = event
        self.ax = []
        self._axis = {}
        self.save = save
        self.gs = GridSpec(8,5) # 8x5 grid of axis space
        self.gs_slice = self.gs[:,1:]
        
        # Hijack buttons
        # -------------------------------------------------------------------#
        # Make forward plot next solution
        forward = NavigationToolbar2.forward
        def next_fm(navtb, *args, **kwargs):
            m = self._fm_index
            n = self._num_fms
            if m < n-1:
                self.plot(solution=m+1)
            forward(navtb, *args, **kwargs)
        NavigationToolbar2.forward = next_fm
        
        # Make back plot previous solution
        back = NavigationToolbar2.back
        def prev_fm(navtb, *args, **kwargs):
            m = self._fm_index
            if m > 0:
                self.plot(solution=m-1)
            back(navtb, *args, **kwargs)
        NavigationToolbar2.back = prev_fm
        
        # Make home plot preferred solution
        home = NavigationToolbar2.home
        def pref_fm(navtb, *args, **kwargs):
            self.plot()
            home(navtb, *args, **kwargs)
        NavigationToolbar2.home = pref_fm
        
        # Take a 'save' function passed on creation and try to call it
        if self.save:
            savefig = NavigationToolbar2TkAgg.save_figure
            def save_fm(navtb, *args, **kwargs):
                try:
                    navtb.set_message("Saving...")
                    self.save(self)
                    navtb.set_message("DONE")
                except Exception as e:
                    navtb.set_message("FAILED:{0}".format(e.message))
                    
            NavigationToolbar2TkAgg.save_figure = save_fm
        # -------------------------------------------------------------------#

        # Draw figure
        self.fig = plt.figure(facecolor='#D9D9EE')
        self.plot()
        
        
        plt.show()
    
