# -*- coding: utf-8 -*-
"""
 plotter.py

 by Mark Williams 2012.313

 Interactive plotter using ObsPy Event as I/O
 (See docstring for class 'FocalMechPlotter')

"""

from numpy import arange
from matplotlib import pyplot as plt
from matplotlib.gridspec import GridSpec
import mplstereonet

class HashPlotter(object):
    """
    Class to create interactive stereonet plot of a HASH first motion
    
    This is a simple, slow, matplotlib script to plot a first motion
    focal mechanism using ObsPy object input. Currently supports:

    Event with a FocalMechanism
                 Origin/Arrivals
                 Picks
    
    <optional>
    Stream containing waveform data for first motion picks
    
    Notes
    -----
    Requires the mplstereonet matplotlib extension.

    Future plan would be to do this in straight Tkinter or
    preferable Qt toolkit.

    """
    fig = None      # handle to figure
    ax = None       # list of axes in the figure
    h = None        # list of handles of FM picks 
    gs = None       # GridSpec instance of plot figure
    l = None        # handle to current waveform line plot
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
        
    def onpick(self, event):
        '''Plot waveform of Z channel when a station is clicked
        
        When a first motion is clicked on the stereonet, get the
        waveform from the database and plot it. This will be the data
        associated with the first motion P pick.
        
        See the 'get_waveform_from_arid' in hashpy.db.utils 
        '''
        N = len(event.ind)
        if not N: return True
        
        if event.artist in self.h:
            self._pickevent = event
        else:
            return True
        
        # Clear figure, set up to plot
        fig = self.fig 
        ax = self.ax[1]
        ax.clear()
        
        # Find the right pick in our data that was clicked on
        #dataind = event.ind[0]              # index of events IN that handle
        dataind = self.h.index(event.artist) # index of which handle
        self._curh = self.h[dataind]         # set current handle
        k = self.ind[dataind]                # index of which arrival that handle plot is
        self._arrv = self._orig.arrivals[k]  # set current arrival
        pick = self._pick
        
        # Data fetch (in here for now)
        if self.data:
            #st = get_waveform_from_arid(self.database, pick.creation_info.version, window=0.5)
            st = self.data.select(station=pick.waveform_id.station_code, channel=pick.waveform_id.channel_code)
            half_window = 0.5/2.
            st = st.slice(pick.time - half_window, pick.time + half_window)
            if len(st) > 0:
                # Plot em if you got em
                self.l, = ax.plot(st[0].data, color=self.wf_color[pick.polarity], lw=2)
                # to stick a vertical line at the pick... try ax.axvline(x=xpick) or axvspan for rect
                xpick = len(st[0].data)/2
                #yt = abs(st[0].data).max()
                #yb = -yt
                v = ax.axvline(x=xpick, ls='--', lw=2, color='black')
                #v, = ax.plot([xpick,xpick],[yb,yt],'--k', lw=2)
            else:
                self.l = None
        else:
            self.l = None
        if self.l is None:
            h_text = ax.text(0.5, 0.5,'No data for this station', 
                horizontalalignment='center',
                verticalalignment='center',
                transform = ax.transAxes)
        ax.set_xlabel("{0} -- {1} -> {2}".format(pick.waveform_id.station_code, pick.creation_info.version, pick.polarity))
        ax.set_xticklabels([])
        ax.set_yticklabels([])
        plt.draw()
        return True
    
    def enter_axes(self, event):
        '''When mouse enters button area'''
        if event.inaxes in self.ax[2:]:
            event.inaxes.patch.set_facecolor('gray')
            event.canvas.draw()
        return True
    
    def leave_axes(self, event):
        '''When mouse leaves button area'''
        if event.inaxes in self.ax[2:]:
            event.inaxes.patch.set_facecolor('white')
            event.canvas.draw()
        return True
    
    def click_on_button(self, event):
        '''Change the first motion designation of a pick
        
        When waveform axes is clicked on, change the polarity to its
        opposite, e.g., from 'UP' to 'DOWN'. The waveform will change
        color correspondingly.
        '''
        if event.inaxes is self._axis['change']:
            # Cycle through 'up, down, exclude' first motion picks
            if self._pick.polarity is not None:
                ind = self.picks_list.index(self._pick.polarity)
                if ind < 2:
                    ind += 1
                else:
                    ind = 0
                self._pick.polarity = self.picks_list[ind]
                if self.l:
                    self.l.set_color(self.wf_color[self.picks_list[ind]])
                    self.ax[1].set_ylabel("{0}".format(self.picks_list[ind]))
                event.canvas.draw()
        elif event.inaxes is self._axis['save']:
            # Save planes/P,T axes/takeoffs/picks to db tables
            event.inaxes.patch.set_facecolor('red')
            event.canvas.draw()
            # should be a referred save function to support external namespaces
            if self.save:
                self.save(self.event, self.source)
            else:
                self.fig.savefig('focal_mech_'+ self._orig.creation_info.version +'.png')
            # done saving, rc=True is success
            event.inaxes.patch.set_facecolor('white')
            event.canvas.draw()
        elif event.inaxes is self._axis['quit']:
            # Close window / go back to previous program
            plt.close(self.fig)
        else:
            pass
        return True
    
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
            plane_str = "STRIKE:{0: > 7.1f}\nDIP:{1: > 7.1f}\nRAKE:{2: > 7.1f}"
            h_text = self.fig.text(0.25, 0.88, plane_str.format(strike1, dip1, rake1), ha='right', va='top', family='monospace')
            h_text = self.fig.text(0.25, 0.33, 'Quality: {0}'.format(qual), ha='right', va='top', family='monospace')  
        
        if not axis:
            self.h = h
            self.ind = index
            
    def _add_button(self, gridslice, textlabel):
        '''Add a button to the thing'''
        ax = self.fig.add_subplot(gridslice) # waveform
        ax.set_xticks([])
        ax.set_yticks([])
        ax.text(0.5, 0.5, textlabel,
            horizontalalignment='center',
            verticalalignment='center',
            transform = ax.transAxes)
        return ax
    
    def __init__(self, event=None, datastream=None, save=None, source=None):
        """
        Create a plot for focal mechanisms

        Input
        -----
        event : obspy.core.event.Event containing 
            FocalMechanism
            Picks
            Origin/Arrivals

        datastream : Stream with Traces - sta/chan names matching Event Picks.
        
        """
        self.event = event
        self.data = datastream
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
        ax = self.fig.add_subplot(self.gs[:-2,:-1], projection='stereonet') # net
        self.ax.append(ax)
        self._axis['stereonet'] = ax
        
        # Waveform display
        ax = self.fig.add_subplot(self.gs[-2:,:])
        ax.set_xticklabels([])
        ax.set_yticklabels([])
        self.h_text = ax.text(0.5, 0.5,'Waveform data - click station',
                    horizontalalignment='center',
                    verticalalignment='center',
                    transform = ax.transAxes)
        self._axis['waveform'] = ax
        self.ax.append(self._axis['waveform'])
        
        # Buttons that do things
        self._axis['quit'] = self._add_button(self.gs[2,-1], 'Return')
        self.ax.append(self._axis['quit'])
        
        self._axis['save'] = self._add_button(self.gs[0:2,-1], 'SAVE') #\nfocal mech\nto db')
        self.ax.append(self._axis['save'])
        
        self._axis['change'] = self._add_button(self.gs[-3,-1], 'Change\ndirection')
        self.ax.append(self._axis['change'])
        
        # Plot first motions
        self.plot_on_stereonet()
        
        # Set mouse events to method functions.
        self.fig.canvas.mpl_connect('pick_event', self.onpick)
        self.fig.canvas.mpl_connect('axes_enter_event', self.enter_axes)
        self.fig.canvas.mpl_connect('axes_leave_event', self.leave_axes)
        self.fig.canvas.mpl_connect('button_press_event', self.click_on_button)
        # Save and draw
        plt.show()
    
