#! /usr/bin/env python3
# -*- coding: utf-8 -*-

'''
This is simple tool for finding cloudlets from a loaded set of spots
Structure:
-> App (QtCore.QApplication)
--> mainWindow (QtWidgets.QMainWindow)
---> window (QtWidgets.QWidget)
----> plot(FigureCanvasQTAgg)
'''

# --- import proper modules ---
from PySide2 import QtCore, QtWidgets, QtGui
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.patches import Ellipse
from matplotlib.widgets import RectangleSelector
import numpy as np
import matplotlib
import sys
# --- custom modules ---
from os.path import realpath
scr_directory = realpath(__file__)
tmp = scr_directory.split('/')
scr_directory = ""
for i in range(len(tmp)-1):
    scr_directory += tmp[i] + '/'
sys.path.append(scr_directory + "/../data/")
from maserVlbi import maserVlbi
from cloudletClass import cloudletClass

matplotlib.use('Qt5Agg')
plt.rc('font', family='serif', style='normal', variant='normal', weight='normal', stretch='normal', size=10)
# -----------------------------

class cloudletFinder(QtWidgets.QApplication):
    def __init__(self, filename):
        '''
        Initializes the App instance
        It likely does a bit too much at the moment:
        -> loads data
        -> declares UI
        -> organizes UI
        -> connects UIelements to slots
        -> fills spot list
        '''
        super().__init__()
        self.data = maserVlbi(filename, verbose=True)

        self.__declareUIElements()
        self.__placeUIelements()
        self.__connectToSlots()
        self.__setPlotTitle(self.data.project_code)
        self.plot.setBeamProps(self.data.beam_raAxis, self.data.beam_decAxis, self.data.beam_posang)
        self.plot.makeChannelPlot(self.data.spots.dRA, self.data.spots.dDEC, self.data.spots.channels)
        self.cloudletsTab = []
        self.__fillSpotsList()
        self.__fillCloudletList()
        self.mainWindow.setGeometry(300, 300, 1366, 720)
        self.mainWindow.show()

    
    def __declareUIElements(self):
        '''
        Simply declares UI elements
        '''
        self.mainWindow = QtWidgets.QMainWindow()
        self.window = QtWidgets.QWidget()
        self.layout = QtWidgets.QGridLayout(self.window)
        self.mainWindow.setCentralWidget(self.window)
        # ----
        self.frame_spots = QtWidgets.QGroupBox("Spots")
        self.frame_spots_vbox = QtWidgets.QVBoxLayout(self.frame_spots)
        self.frame_cloudlets = QtWidgets.QGroupBox("Cloudlets")
        self.frame_cloudlets_vbox = QtWidgets.QVBoxLayout(self.frame_cloudlets)
        # ---
        self.frame_right_hand = QtWidgets.QGroupBox("Operations")
        self.frame_right_hand_vbox = QtWidgets.QVBoxLayout(self.frame_right_hand)
        # ---
        self.plot = plotCanvas()
        self.plotLayout = QtWidgets.QVBoxLayout()
        self.plotToolbar = NavigationToolbar(self.plot, self.mainWindow)
        self.plotLayout.addWidget(self.plotToolbar)
        self.plotLayout.addWidget(self.plot)
        # -- buttons --
        self.exitButton = QtWidgets.QPushButton(self.window)
        self.spotList = QtWidgets.QListWidget()
        self.cloudletList = QtWidgets.QListWidget()
        self.addToCloudlets = QtWidgets.QPushButton(self.window)
        self.removeFromCloudlets = QtWidgets.QPushButton(self.window)
        self.saveCloudlets = QtWidgets.QPushButton(self.window)
        self.clearCloudletInfo = QtWidgets.QPushButton(self.window)
        # -- checkboxes --
        self.showBeam = QtWidgets.QCheckBox(self.window)
        self.showMarkRectangle = QtWidgets.QCheckBox(self.window)
        self.showSpots = QtWidgets.QCheckBox(self.window)
        self.showCloudlets = QtWidgets.QCheckBox(self.window)
        self.showChannels = QtWidgets.QCheckBox(self.window)
        self.showMarker = QtWidgets.QCheckBox(self.window)
        # --
        self.addToCloudlets.setText("Add to cloudlets")
        self.removeFromCloudlets.setText("Remove from cloudlets")
        self.showBeam.setText("Show beam")
        self.showMarkRectangle.setText("Show mark range")
        self.showSpots.setText("Show spots")
        self.showCloudlets.setText("Show cloudlets")
        self.showChannels.setText("Show channels")
        self.showMarker.setText("Show spot marker")
        self.saveCloudlets.setText("Save cloudlets")
        self.clearCloudletInfo.setText("Clear cloudlet info")
        # -
        self.showMarkRectangle.setChecked(True)
        self.showSpots.setChecked(True)
        self.showMarker.setChecked(True)
        # --------------
        self.exitButton.setText("Exit")

    def __placeUIelements(self):
        '''
        Simply places UI elements correctly
        '''
        self.frame_spots_vbox.addWidget(self.spotList)
        self.frame_cloudlets_vbox.addWidget(self.cloudletList)
        # --
        self.frame_right_hand_vbox.addWidget(self.addToCloudlets)
        self.frame_right_hand_vbox.addWidget(self.removeFromCloudlets)
        self.frame_right_hand_vbox.addWidget(self.saveCloudlets)
        self.frame_right_hand_vbox.addWidget(self.clearCloudletInfo)
        self.frame_right_hand_vbox.addWidget(self.showBeam)
        self.frame_right_hand_vbox.addWidget(self.showMarkRectangle)
        self.frame_right_hand_vbox.addWidget(self.showSpots)
        self.frame_right_hand_vbox.addWidget(self.showCloudlets)
        self.frame_right_hand_vbox.addWidget(self.showCloudlets)
        self.frame_right_hand_vbox.addWidget(self.showChannels)
        self.frame_right_hand_vbox.addWidget(self.showMarker)
        self.frame_right_hand_vbox.addWidget(self.exitButton)
        # --
        self.layout.addWidget(self.frame_spots, 0,0, 5,1)
        self.layout.addWidget(self.frame_cloudlets, 5,0, 5,1)
        self.layout.addLayout(self.plotLayout, 0,1, 10,3)
        self.layout.addWidget(self.frame_right_hand, 0,4, 10,1)

        for i in range(self.layout.rowCount()):
            self.layout.setRowStretch(i,1)
        for i in range(self.layout.columnCount()):
            self.layout.setColumnStretch(i,1)

    def __connectToSlots(self):
        '''
        Here every connection to button / checkbox / list / plot is gathered
        Nothing special about this method tbh.
        '''
        self.exitButton.clicked.connect(QtWidgets.QApplication.exit)
        self.showBeam.clicked.connect(self.__beam_visible)
        self.kk = self.plot.fig.canvas.mpl_connect('button_press_event', self.plot.onClick)
        self.showMarkRectangle.clicked.connect(self.__rectangle_visible)
        self.showSpots.clicked.connect(self.__spots_visible)
        self.showChannels.clicked.connect(self.__channels_visible)
        self.showCloudlets.clicked.connect(self.__cloudletVisibilitySwitchSlot)
        self.addToCloudlets.clicked.connect(self.__addToCloudletList)
        self.removeFromCloudlets.clicked.connect(self.__removeFromCloudletList)
        self.saveCloudlets.clicked.connect(self.__saveCloudletsToFile)
        self.clearCloudletInfo.clicked.connect(self.__clearAllCloudletInfo)
        #QtCore.QObject.connect(self.spotList, QtCore.SIGNAL('currentItemChanged()'), self.__plotMarkerOnClick)
        self.spotList.currentItemChanged.connect(self.__plotMarkerOnClick)
        self.showMarker.clicked.connect(self.__markerVisibility)

    def plotSpots(self):
        '''
        Spots are plotted ONLY after invoking this method
        '''
        b = self.plot.addSpotPlot(self.data.spots, self.data.spectrum)
        self.plot.plotCoudlets(self.data)
        self.plot.draw()
        self.spotList.setCurrentRow(0)
    

    def __setPlotTitle(self, title):
        '''
        Dunno what is this method doing here tbh.
        '''
        self.plot.axSpots.set_title(title)

    def getSpotsFromRange(self, ranges):
        '''
        Returns attributes from spots, that are placed within ranges specified in the
        arg.
        range should be: [x_min, x_max, y_min, y_max]
        Returns 8 x numpy.ndarray
        '''
        xmin, xmax, ymin, ymax = ranges
        dRAinRange = []
        dRA_errinRange = []
        dDECnRange = []
        dDEC_errinRange = []
        fluxinRange = []
        flux_errinRange = []
        channelsinRange = []
        velocityinRange = []
        spots = self.data.spots
        for i in range(spots.dRA.size):
            if spots.dRA[i] > xmin and spots.dRA[i] < xmax and spots.dDEC[i] > ymin and spots.dDEC[i] < ymax:
                dRAinRange.append(spots.dRA[i])
                dRA_errinRange.append(spots.dRA_err[i])
                dDECnRange.append(spots.dDEC[i])
                dDEC_errinRange.append(spots.dDEC_err[i])
                fluxinRange.append(spots.flux[i])
                flux_errinRange.append(spots.flux_err[i])
                channelsinRange.append(spots.channels[i])
                velocityinRange.append(spots.velocity[i])
        return np.asarray(dRAinRange), np.asarray(dRA_errinRange), np.asarray(dDECnRange), np.asarray(dDEC_errinRange), np.asarray(fluxinRange), np.asarray(flux_errinRange), np.asarray(channelsinRange), np.asarray(velocityinRange)

    def addToList(self, cloudlet):
        '''
        Adds entry to cloudletList
        '''
        self.cloudletList.addItem(f'({round(cloudlet.dRA, 2)},{round(cloudlet.dDEC, 2)}): {round(cloudlet.maxFlux, 2)} Jy/beam')

    def __fillSpotsList(self):
        '''
        Fills spotList with entries about spots themselves
        '''
        for ra, dec, vel, flux in zip(self.data.spots.dRA, self.data.spots.dDEC, self.data.spots.velocity, self.data.spots.flux):
            self.spotList.addItem(f'({round(ra, 2)},{round(dec, 2)}): {round(vel, 2)} km/s, {round(flux, 2)} Jy / beam ')
        #self.spotList.setCurrentRow(0)

    def __fillCloudletList(self):
        '''
        Fills cloudletList with entries about cloudlets
        '''
        if len(self.data.cloudlets) > 0:
            for cloudlet in self.data.cloudlets:
                self.addToList(cloudlet)
            self.cloudletList.setCurrentRow(0)

    '''
    ==========================
    == BELOW WE STORE SLOTS ==
    ==========================
    '''

    '''
    BEAM AND RECT. SELECTOR VISIBILITY
    '''
    def __beam_visible(self):
        '''
        Sets beam visible or not - since beam and rect. zoom cannot be visible at 
        the same time due to left mouse click conflict, we handle rect. zoom 
        visibility to some extend as well
        '''
        if not self.showBeam.isChecked():
            self.plot.beam_ellipse.set_visible(False)
            self.plot.fig.canvas.mpl_disconnect(self.kk)
        else:
            self.plot.beam_ellipse.set_visible(True)
            self.kk = self.plot.fig.canvas.mpl_connect('button_press_event', self.plot.onClick)
            if self.plot.selector.active:
                self.plot.selector.set_active(False)
                self.showMarkRectangle.setChecked(False)
                self.plot.selector.set_visible(False)
        self.plot.draw()
    def __rectangle_visible(self):
        '''
        Similar to the __beam_visible, but this handles rectangle visibility and
        being active
        '''
        if not self.showMarkRectangle.isChecked():
            self.plot.selector.set_active(False)
            self.plot.selector.set_visible(False)
        else:
            if self.showBeam.isChecked():
                self.showBeam.setChecked(False)
                self.__beam_visible()
            self.plot.selector.set_active(True)
            self.plot.selector.set_visible(True)
        self.plot.fig.canvas.draw_idle()
    
    '''
    SPOTS AND CHANNEL VISIBILITY
    '''
    def __spots_visible(self):
        self.plot.spotScatterPlot.set_visible(self.showSpots.isChecked())
        self.plot.fig.canvas.draw_idle()
    def __channels_visible(self):
        self.plot.setChannelLabelsVisible(self.showChannels.isChecked())

    '''
    CLOUDLET ADDING AND REMOVING
    '''
    def __addToCloudletList(self):
        '''
        Adds cloudlet to the list - takes ranges from selector, searches for spots 
        within ranges and creates classCloudlet object, that holds gathered data.
        In the end - adds this data to the main "data" object (class: maserVlbi)
        If there is no spots in marked range, this slot should do nothing (return)
        '''
        if not self.plot.selector.active:
            print("---> selector not active!")
            return
        # -
        ranges = self.plot.getMarkedRange()
        spotsInRange = self.getSpotsFromRange(ranges)
        if len(spotsInRange[0]) < 1:
            return
        # -
        cloudlet = cloudletClass()
        cloudlet.setAttributes(*spotsInRange)
        cloudlet.calcProps()
        # -
        self.addToList(cloudlet) # <--- ADDS ONLY TO THE INTERFACE
        self.data.appendCloudlet(cloudlet) # <--- THIS IS TO HOLD ACTUAL DATA
        self.plot.setNewCloudletPlot(self.data, self.showCloudlets.isChecked()) # <--- UPDATES GRAPH
    def __removeFromCloudletList(self):
        '''
        Simply removes cloudlet from the list
        -> from the interface
        -> and from the "data" object
        '''
        try:
            index = self.cloudletList.currentRow()
            self.data.removeCloudlet(index)
            self.cloudletList.takeItem(index)
            self.plot.setNewCloudletPlot(self.data, self.showCloudlets.isChecked())
        except:
            return    
    def __clearAllCloudletInfo(self):
        '''
        Wipes all cloudlet information:
        -> from the plot
        -> from the list
        -> from dle << data >> class
        -> from the .hdf5 file
        '''
        self.cloudletList.clear()
        self.data.clearCloudletInfo()
        self.plot.setNewCloudletPlot(self.data, self.showCloudlets.isChecked())
        print('---> Removed info about cloudlets')

    def __saveCloudletsToFile(self):
        self.data.saveCloudlets()

    def __cloudletVisibilitySwitchSlot(self):
        self.plot.cloudletVisible(self.showCloudlets.isChecked())

    '''
    SPOT MARKER
    '''
    def __plotMarkerOnClick(self):
        '''
        This slot sets circle marker on the spot, chosen from spotList
        Very simple
        '''
        x,y,flux,vel = self.data.spots.getPropsFromIndex(self.spotList.currentRow())
        self.plot.plotMarker(x,y,flux)
    def __markerVisibility(self):
        '''
        This slot sets marker visibility
        '''
        self.plot.markerPlot.set_visible(self.showMarker.isChecked())
        self.plot.fig.canvas.draw_idle()

class plotCanvas(FigureCanvas):
    def __init__(self):
        self.fig = plt.figure()
        super(plotCanvas, self).__init__(self.fig)
        self.gs = gridspec.GridSpec(1,2, width_ratios=[30,1], figure=self.fig, wspace=0.0)
        self.__declareNecessaryAxes()
        self.__declareCustomWidgets()
    def __declareNecessaryAxes(self):
        self.axSpots = self.fig.add_subplot(self.gs[0,0])
        self.axSpots.invert_xaxis()
        self.axCbar = self.fig.add_subplot(self.gs[0,1])
        self.__makeFancyTicks(self.axSpots)
        self.__makeFancyTicks(self.axCbar)
        self.axSpots.set_xlabel("$\Delta$RA (mas)")
        self.axSpots.set_ylabel("$\Delta$DEC (mas)")

    def addSpotPlot(self, spots, spectrum):
        '''
        Plots spots
        Designed to be invoked only once, upon startup
        '''
        colors = spots.getJetColors(spectrum.velocity.min(), spectrum.velocity.max())
        self.spotScatterPlot = self.axSpots.scatter(spots.dRA, spots.dDEC, s=np.log(spots.flux*1000.0)**2.0 * 10, edgecolor = colors, facecolor='none')
        self.markerPlot = self.axSpots.scatter(spots.dRA[0], spots.dDEC[0], s=np.log(spots.flux[0]*1000.0)**2.0 * 10, edgecolor = 'none', facecolor='grey')
        p = plt.cm.ScalarMappable(cmap=plt.cm.jet)
        p.set_array(spectrum.velocity)
        self.fig.colorbar(p, ax=self.axSpots, cax = self.axCbar, label="V$_{LSR}\,$(km$\,$s$^{-1}$)")
        self.axCbar.set_ylim(spectrum.velocity.min(), spectrum.velocity.max())
        self.__autoscaleSpotPlot(spots)
        return self.spotScatterPlot
    
    def __autoscaleSpotPlot(self, spots):
        diffra = abs(spots.dRA.max() - spots.dRA.min())
        diffdec = abs(spots.dDEC.max() - spots.dDEC.min())

        self.axSpots.set_xlim(spots.dRA.min() - 0.05 * diffra, spots.dRA.max() + 0.05 * diffra)
        self.axSpots.set_ylim(spots.dDEC.min() - 0.05 * diffdec, spots.dDEC.max() + 0.05 * diffdec)
        self.axSpots.invert_xaxis()

    def plotCoudlets(self, data):
        '''
        Plot cloudlets. 
        Designed to be invoked only once, upon startup
        '''
        colors = data.getClJetColors(data.spectrum.velocity.min(), data.spectrum.velocity.max())
        self.cloudletPlots = []
        for index, cloudlet in enumerate(data.cloudlets):
            self.cloudletPlots.append(self.axSpots.scatter(cloudlet.dRA, cloudlet.dDEC, s=np.log(cloudlet.maxFlux*1000.0)**2.0 * 10, edgecolor = 'none', facecolor=colors[index], marker='P', visible=False) )

    
    def __makeFancyTicks(self, ax):
        ax.xaxis.set_tick_params(direction='in', width=1, length = 3, top = True, bottom=True)
        ax.xaxis.set_tick_params(direction='in', width=1, length = 3, which='minor', top = True, bottom=True)
        ax.yaxis.set_tick_params(direction='in', width=1, length = 3, right=True, left=True)
        ax.yaxis.set_tick_params(direction='in', width=1, length = 3, which='minor', right=True, left=True)

    def __declareCustomWidgets(self):
        '''
        This is simple methood that enables beam plot and rectangle selector
        '''
        self.selected_range = [0,0,0,0]
        self.beam_ellipse = Ellipse([0,0], 4.5, 4.5, angle=0.0, fc='none', ec='black', visible=False)
        self.selector = RectangleSelector(self.axSpots, self.line_select_callback, useblit=True, button=[1,3], minspanx=2, minspany=2, spancoords='pixels', interactive=True)
        self.axSpots.add_patch(self.beam_ellipse)

    def line_select_callback(self, eclick, erelease):
        # 'eclick' and 'erelease' are just mouse right button press and release signals
        x1, y1 = eclick.xdata, eclick.ydata
        x2, y2 = erelease.xdata, erelease.ydata

        # setting global table of ranges
        self.selected_range = [x1,x2,y1,y2]
    
    def makeChannelPlot(self, x,y,channels):
        self.spotChannelLabels = []
        for index, channel in enumerate(channels):
            w = self.axSpots.text(x[index], y[index], str(channel), visible=False)
            self.spotChannelLabels.append(w)

    '''
    METHOODS FOR PARTS OF THE PLOT
    '''
    def setBeamProps(self, x_axis, y_axis, posang):
        #self.beam_ellipse.update_from(Ellipse(self.beam_ellipse.get_center(), x_axis, y_axis, angle=posang, ec='black', fc='none', visible=False))
        #self.axSpots.add_patch(self.beam_ellipse)
        self.beam_ellipse.set_width(x_axis)
        self.beam_ellipse.set_height(y_axis)
        self.beam_ellipse.angle = posang
    
    def onClick(self, event):
        x,y = event.xdata, event.ydata
        if x == None or y==None:
            return
        self.beam_ellipse.set_center([x,y])
        self.fig.canvas.draw_idle()
    
    def setChannelLabelsVisible(self, flag):
        [i.set_visible(flag) for i in self.spotChannelLabels]
        self.fig.canvas.draw_idle()
    
    def getMarkedRange(self):
        return self.selector.extents
    
    def plotMarker(self, x, y, flux):
        self.markerPlot.set_offsets([x,y])
        self.markerPlot.set_sizes([np.log(flux*1000.0)**2.0 * 10])
        self.fig.canvas.draw_idle()
    
    def cloudletVisible(self, flag):
        [clPlot.set_visible(flag) for clPlot in self.cloudletPlots]
        self.fig.canvas.draw_idle()
    
    def setNewCloudletPlot(self, data, flag):
        '''
        Destroys previous cloudlet plot and creates new one
        '''
        # - 
        for i in self.cloudletPlots:
            i.remove()
        # -
        colors = data.getClJetColors(data.spectrum.velocity.min(), data.spectrum.velocity.max())
        self.cloudletPlots = []
        for index, cloudlet in enumerate(data.cloudlets):
            self.cloudletPlots.append(self.axSpots.scatter(cloudlet.dRA, cloudlet.dDEC, s=np.log(cloudlet.maxFlux*1000.0)**2.0 * 10, edgecolor = 'none', facecolor=colors[index], marker='P', visible=flag) )
        self.fig.canvas.draw_idle()

if __name__ == '__main__':
    app = cloudletFinder(sys.argv[1])
    app.plotSpots()
    app.exec_()