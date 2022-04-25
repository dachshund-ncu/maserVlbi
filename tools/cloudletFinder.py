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
        super().__init__()
        self.data = maserVlbi(filename)

        self.__declareUIElements()
        self.__placeUIelements()
        self.__connectToSlots()
        self.__setPlotTitle(self.data.project_code)
        self.plot.setBeamProps(self.data.beam_raAxis, self.data.beam_decAxis, self.data.beam_posang)
        self.plot.makeChannelPlot(self.data.spots.dRA, self.data.spots.dDEC, self.data.spots.channels)
        self.cloudletsTab = []
        #self.plot.setBeamProps(25,25,90)
        self.mainWindow.setGeometry(300, 300, 1366, 720)
        self.mainWindow.show()

    
    def __declareUIElements(self):
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
        # -- checkboxes --
        self.showBeam = QtWidgets.QCheckBox(self.window)
        self.showMarkRectangle = QtWidgets.QCheckBox(self.window)
        self.showSpots = QtWidgets.QCheckBox(self.window)
        self.showCloudlets = QtWidgets.QCheckBox(self.window)
        self.showChannels = QtWidgets.QCheckBox(self.window)
        # --
        self.addToCloudlets.setText("Add to cloudlets")
        self.removeFromCloudlets.setText("Remove from cloudlets")
        self.showBeam.setText("Show beam")
        self.showMarkRectangle.setText("Show mark range")
        self.showSpots.setText("Show spots")
        self.showCloudlets.setText("Show cloudlets")
        self.showChannels.setText("Show channels")
        # -
        self.showMarkRectangle.setChecked(True)
        self.showSpots.setChecked(True)
        # --------------
        self.exitButton.setText("Exit")

    def __placeUIelements(self):
        self.frame_spots_vbox.addWidget(self.spotList)
        self.frame_cloudlets_vbox.addWidget(self.cloudletList)
        # --
        self.frame_right_hand_vbox.addWidget(self.addToCloudlets)
        self.frame_right_hand_vbox.addWidget(self.removeFromCloudlets)
        self.frame_right_hand_vbox.addWidget(self.showBeam)
        self.frame_right_hand_vbox.addWidget(self.showMarkRectangle)
        self.frame_right_hand_vbox.addWidget(self.showSpots)
        self.frame_right_hand_vbox.addWidget(self.showCloudlets)
        self.frame_right_hand_vbox.addWidget(self.showCloudlets)
        self.frame_right_hand_vbox.addWidget(self.showChannels)
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
        self.exitButton.clicked.connect(QtWidgets.QApplication.exit)
        self.showBeam.clicked.connect(self.__beam_visible)
        self.kk = self.plot.fig.canvas.mpl_connect('button_press_event', self.plot.onClick)
        self.showMarkRectangle.clicked.connect(self.__rectangle_visible)
        self.showSpots.clicked.connect(self.__spots_visible)
        self.showChannels.clicked.connect(self.__channels_visible)
        self.addToCloudlets.clicked.connect(self.__addToCloudletList)
    def plotSpots(self):
        b = self.plot.addSpotPlot(self.data.spots, self.data.spectrum)
        self.plot.draw()
    

    def __setPlotTitle(self, title):
        self.plot.axSpots.set_title(title)

    def getSpotsFromRange(self, ranges):
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
        self.cloudletList.addItem(f'({round(cloudlet.dRA, 2)},{round(cloudlet.dDEC, 2)}): {round(cloudlet.maxFlux, 2)} Jy/beam')

    '''
    BELOW WE STORE SLOTS
    '''

    def __beam_visible(self):
        if not self.showBeam.isChecked():
            self.plot.beam_ellipse.set_visible(False)
            self.plot.fig.canvas.mpl_disconnect(self.kk)
        else:
            self.plot.beam_ellipse.set_visible(True)
            self.kk = self.plot.fig.canvas.mpl_connect('button_press_event', self.plot.onClick)
            if self.plot.selector.active:
                self.plot.selector.set_active(False)
                self.showMarkRectangle.setChecked(False)
        self.plot.draw()

    def __rectangle_visible(self):
        # setting span selector active / inactive
        if not self.showMarkRectangle.isChecked():
            self.plot.selector.set_active(False)
        else:
            if self.showBeam.isChecked():
                self.showBeam.setChecked(False)
                self.__beam_visible()
            self.plot.selector.set_active(True)
        self.plot.fig.canvas.draw_idle()
    
    def __spots_visible(self):
        self.plot.spotScatterPlot.set_visible(self.showSpots.isChecked())
        self.plot.fig.canvas.draw_idle()
    
    def __channels_visible(self):
        self.plot.setChannelLabelsVisible(self.showChannels.isChecked())

    def __addToCloudletList(self):
        if not self.plot.selector.active:
            print("selector not active!")
            return
        ranges = self.plot.getMarkedRange()
        #dRA, dRA_err, dDEC, dDEC_err, flux, flux_err, channels, velocity = self.getSpotsFromRange(ranges)
        cloudlet = cloudletClass()
        cloudlet.setAttributes(*self.getSpotsFromRange(ranges))
        cloudlet.calcProps()
        self.addToList(cloudlet)
        #cloudlet = self.makeCloudlet(dRA, dRA_err, dDEC, dDEC_err, flux, flux_err, channels, velocity)
        

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
        colors = spots.getJetColors(spectrum.velocity.min(), spectrum.velocity.max())
        self.spotScatterPlot = self.axSpots.scatter(spots.dRA, spots.dDEC, s=np.log(spots.flux*1000.0)**2.0 * 10, edgecolor = colors, facecolor='none')
        p = plt.cm.ScalarMappable(cmap=plt.cm.jet)
        p.set_array(spectrum.velocity)
        self.fig.colorbar(p, ax=self.axSpots, cax = self.axCbar, label="V$_{LSR}\,$(km$\,$s$^{-1}$)")
        self.axCbar.set_ylim(spectrum.velocity.min(), spectrum.velocity.max())
        return self.spotScatterPlot
    
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


if __name__ == '__main__':
    app = cloudletFinder(sys.argv[1])
    app.plotSpots()
    app.exec_()