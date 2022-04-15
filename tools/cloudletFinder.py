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
from matplotlib.figure import Figure
import matplotlib.gridspec as gridspec
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
        # --
        self.addToCloudlets.setText("Add to cloudlets")
        self.removeFromCloudlets.setText("Remove from cloudlets")
        self.showBeam.setText("Show beam")
        self.showMarkRectangle.setText("Show mark range")
        self.showSpots.setText("Show spots")
        self.showCloudlets.setText("Show cloudlets")

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
    

    def plotSpots(self):
        b = self.plot.addSpotPlot(self.data.spots, self.data.spectrum)
        self.plot.draw()


class plotCanvas(FigureCanvas):
    def __init__(self):
        self.fig = plt.figure()
        super(plotCanvas, self).__init__(self.fig)
        self.gs = gridspec.GridSpec(1,2, width_ratios=[30,1], figure=self.fig, wspace=0.0)
        self.__declareNecessaryAxes()
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
        b = self.axSpots.scatter(spots.dRA, spots.dDEC, s=np.log(spots.flux*1000.0)**2.0 * 10, edgecolor = colors, facecolor='none')
        p = plt.cm.ScalarMappable(cmap=plt.cm.jet)
        p.set_array(spectrum.velocity)
        self.fig.colorbar(p, ax=self.axSpots, cax = self.axCbar, label="V$_{LSR}\,$(km$\,$s$^{-1}$)")
        self.axCbar.set_ylim(spectrum.velocity.min(), spectrum.velocity.max())
        return b
    
    def __makeFancyTicks(self, ax):
        ax.xaxis.set_tick_params(direction='in', width=1, length = 3, top = True, bottom=True)
        ax.xaxis.set_tick_params(direction='in', width=1, length = 3, which='minor', top = True, bottom=True)
        ax.yaxis.set_tick_params(direction='in', width=1, length = 3, right=True, left=True)
        ax.yaxis.set_tick_params(direction='in', width=1, length = 3, which='minor', right=True, left=True)


if __name__ == '__main__':
    app = cloudletFinder(sys.argv[1])
    app.plotSpots()
    app.exec_()