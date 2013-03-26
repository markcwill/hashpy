# -*- coding: utf-8 -*-
#
# plotter.py
#
# by Mark Williams 2012.313
#
# Interactive plotter using ObsPy Event as I/O
#

from numpy import arange
from matplotlib import pyplot as plt
from matplotlib.gridspec import GridSpec
import mplstereonet
# temporary, should pass data in as Stream, call db in script...
from hashpy.db.utils import readANTELOPE

class FocalMechPlotter(object):
    '''Class to create interactive stereonet plot of a HASH first motion'''
    fig = None      # handle to figure
    ax = None       # list of axes
    h = None        # handles to replace up/down 
    gs = None       # GridSpec instance of plot figure
    l = None        # handle to current waveform line plot
    h_text = None   # handle for text labels
    event = None    # Event object
    ind = None      # indicies of which picks in Origin.arrivals are plotted
    o = None        # current preferred origin of event (change to _o)
    fm = None       # current preferred focal mech of event
    p = None        # current Pick of current Arrival
    a = None        # current Arrival plotted
    database = None # placeholder for database to grab waveforms
    window = None   # window width for waveform plot (0.75 sec, 2/3)
    _arrv = None     # current arrival
    
    picks_list = ('undecidable', 'positive', 'negative')
    fm_list = ('..','c.', 'd.')
    wf_color = { 'positive' : 'red', 'negative' : 'blue', 'undecidable' : 'black' }
    
    @property
    def _orig(self):
        '''Stub'''
        return self.event.preferred_origin()
    @property
    def _focm(self):
        '''Stub'''
        return self.event.preferred_focal_mechanism()
    @property
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
            pass
        else:
            return True
        # Clear figure, set up to plot
        fig = self.fig 
        ax = self.ax[1]
        ax.clear()
        # Find the right pick in our data that was clicked on
        dataind = event.ind[0]
        k = self.ind[dataind]
        arrival = self.o.arrivals[k]
        self._arrv = arrival
        pick = arrival.pick_id.getReferredObject()
        # Set up and Plot (data fetch in here for now)
        ax.set_xlabel("{0} -- {1} -> {2}".format(pick.waveform_id.station_code, pick.creation_info.version, pick.polarity))
        ax.set_xticklabels([])
        ax.set_yticklabels([])
        #st = get_waveform_from_arid(self.database, pick.creation_info.version, window=0.5)
        st = readANTELOPE(self.database, station=pick.waveform_id.station_code, channel=pick.waveform_id.channel_code, time=pick.time-self.window/2., endtime=pick.time+self.window)
        l, = ax.plot(st[0].data, color=self.wf_color[pick.polarity], lw=2)
        xpick = len(st[0].data)/2
        yt = abs(st[0].data).max()
        yb = -yt
        v, = ax.plot([xpick,xpick],[yb,yt],'--k', lw=2)
        self.l = l
        self.p = pick
        self.a = arrival
        plt.show()
        return True
    
    def enter_axes(self, event):
        '''When mouse enters waveform window'''
        if event.inaxes in self.ax[2:]:
            event.inaxes.patch.set_facecolor('gray')
            event.canvas.draw()
    
    def leave_axes(self, event):
        '''When mouse leaves waveform window'''
        if event.inaxes in self.ax[2:]:
            event.inaxes.patch.set_facecolor('white')
            event.canvas.draw()
    
    def click_on_button(self, event):
        '''Change the first motion designation of a pick
        
        When waveform axes is clicked on, change the polarity to its
        opposite, e.g., from 'UP' to 'DOWN'. The waveform will change
        color correspondingly.
        '''
        if event.inaxes is self.ax[-1]:
            # Cycle through 'up, down, exclude' first motion picks
            if self.p.polarity is not None:
                ind = self.picks_list.index(self.p.polarity)
                if ind < 2:
                    ind += 1
                else:
                    ind = 0
                self.p.polarity = self.picks_list[ind]
                if self.l:
                    self.l.set_color(self.wf_color[self.picks_list[ind]])
                    self.ax[1].set_ylabel("{0}".format(self.picks_list[ind]))
                    #change_arrival_fm(database, self.p.creation_info.version, self.fm_list[ind])
                event.canvas.draw()
        elif event.inaxes is self.ax[-2]:
            # Save planes/P,T axes/takeoffs/picks to db tables
            event.inaxes.patch.set_facecolor('red')
            event.canvas.draw()
            #focalmech2db(self.mech)
            event.inaxes.patch.set_facecolor('white')
            event.canvas.draw()
        elif event.inaxes is self.ax[-3]:
            # Close window / go back to previous program
            plt.close(self.fig)
        else:
            pass
    
    def plot_on_stereonet(self):
        '''Plot first motions from an Event instance
        
        Stored in self.event
        '''
        ax = self.ax[0]
        # pull out variables from mechanism
        #--- HASH takeoffs are 0-180 from vertical UP!!
        #--- Stereonet angles 0-90 inward (dip)
        #--- Classic FM's are toa from center???
        #takeoffs = abs(self.mech.picks['takeoff'] - 90)
        #polarities = self.mech.picks['polarity']
        #azimuths = self.mech.picks['azimuth']
        
        self.h = []
        self.ind = []
       
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
            p = a.pick_id.getReferredObject()
            azi = a.azimuth - 90.
            toa = abs(a.takeoff_angle - 90)
            if p.polarity == 'positive':
                plot_specs.update({'markeredgecolor' : 'black',
                                   'markerfacecolor' : 'red'   })
            elif p.polarity == 'negative':
                plot_specs.update({'markeredgecolor' : 'blue',
                                   'markerfacecolor' : 'white' })
            else:
                continue
            _h = ax.rake(azi, toa, 90, 'o', **plot_specs)
            self.h.append(_h)
            self.ind.append(ind)
            if True:
                h_text = ax.rake(azi, toa+5, 90, marker='$   {0}$'.format(p.waveform_id.station_code), color='black',markersize=20)
    
    def __init__(self, event=None, datastream=None, window=0.5, database=None):
        '''Create a plot for focal mechanisms'''
        self.database = database # take out later
        self.event = event
        self.window = window
        self.data = datastream
        
        # Draw figure and set up 
        fig = plt.figure(facecolor='#D9D9EE')
        fig.canvas.set_window_title('Focal Mechanism - dbhash')
        gs = GridSpec(8,5) # 8x5 grid of axis space
        
        # Stereonet axis
        ax = fig.add_subplot(gs[:-2,:-1], projection='stereonet') # net
        ax.clear() 
        ax.set_title('Origin: {0}'.format(self._orig.creation_info.version))
        tlab  = ax.set_azimuth_ticklabels([])
        
        # add to plotting object
        self.fig = fig
        self.ax = [ax]
        self.gs = gs
        
        # Set up other axes (buttons and wf display area)
        self.ax.append(fig.add_subplot(self.gs[-2:,:])) # waveform
        self.ax[-1].set_xticklabels([])
        self.ax[-1].set_yticklabels([])
        self.ax[-1].text(0.5, 0.5,'Waveform data - click station',
                    horizontalalignment='center',
                    verticalalignment='center',
                    transform = self.ax[-1].transAxes)
                    
        self.ax.append(fig.add_subplot(self.gs[2,-1]))  # button
        self.ax[-1].set_xticks([])
        self.ax[-1].set_yticks([])
        self.ax[-1].text(0.5, 0.5,'Return',
                    horizontalalignment='center',
                    verticalalignment='center',
                    transform = self.ax[-1].transAxes)
                    
        self.ax.append(fig.add_subplot(self.gs[0:2,-1]))  # button
        self.ax[-1].set_xticks([])
        self.ax[-1].set_yticks([])
        self.ax[-1].text(0.5, 0.5,'SAVE\nfocal mech\nto db',
                    horizontalalignment='center',
                    verticalalignment='center',
                    transform = self.ax[-1].transAxes)
                    
        self.ax.append(fig.add_subplot(self.gs[-3,-1]))  # button
        self.ax[-1].set_xticks([])
        self.ax[-1].set_yticks([])
        self.ax[-1].text(0.5, 0.5,'Change\ndirection',
                    horizontalalignment='center',
                    verticalalignment='center',
                    transform = self.ax[-1].transAxes)
                    
        
        # Plot first motions
        self.plot_on_stereonet()
        
        # Set mouse events to method functions.
        fig.canvas.mpl_connect('pick_event', self.onpick)
        fig.canvas.mpl_connect('axes_enter_event', self.enter_axes)
        fig.canvas.mpl_connect('axes_leave_event', self.leave_axes)
        fig.canvas.mpl_connect('button_press_event', self.click_on_button)
        # Save and draw
        plt.show()
    
