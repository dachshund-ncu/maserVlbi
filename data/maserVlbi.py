#! /usr/bin/env python3
# -*- coding: utf-8 -*-
'''
This is main class of the maser data
Mainly it reads the HDF5 files
'''

import h5py
from spotClass import spotsClass
from spectrumClass import spectrumClass
from astropy.time import Time
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.collections import LineCollection
import datetime

class maserVlbi:

    def __init__(self, filename, verbose=False):
        self.verbose = verbose
        self._filename = filename

        if self.verbose:
            print('-------------------------------')
            print(f'---> Loading file \"{self._filename}\"...')

        self._fle = h5py.File(filename, 'r')
        # -- reading spots data --
        self.spots = spotsClass()
        self.spots.readDataFromGroup(self._fle['SPOTS'])
        # -- reading spec data --
        self.spectrum = spectrumClass(self._fle['SPECTRUM'])
        # -- optional: --
        self._getDate(self._fle)
        self._getBeam(self._fle)
        self._getProjectCode(self._fle)
        self._getBandData(self._fle)
        self._getSigmaData(self._fle)
        self._getOrigin(self._fle)

        if self.verbose:
            print('-------------------------------')

    '''
    PRIVATE:
    '''
    def _getDate(self, fle):
        try:
            dset = fle['DATE']
            times = [n.decode("ascii") for n in dset]
            times = Time(times, scale='utc', format='isot')
            self.mjd_begin = times[0].mjd
            self.mjd_end = times[1].mjd
            self.isot_begin = times[0].isot
            self.isot_end = times[1].isot
            if self.verbose:
                print(f'---> begin: {self.isot_begin}, end: {self.isot_end}')
        except:
            if self.verbose:
                print(f"---> no time data found!")
    
    def _getBeam(self, fle):
        try:
            dset = fle['BEAM']
            bmarray = np.array(dset)
            
            self.beam_raAxis = bmarray[0]
            self.beam_decAxis = bmarray[1]
            self.beam_posang = bmarray[2]
            if self.verbose:
                print(f'---> beam: {self.beam_raAxis} x {self.beam_decAxis}, {self.beam_posang}')
        except:
            if self.verbose:
                print(f"---> no beam data found!")
    
    def _getProjectCode(self, fle):
        try:
            dset = fle['TELESCOPE']
            self.array = dset[0].decode("ascii")
            self.project_code = dset[1].decode("ascii")
            self.pi = dset[2].decode("ascii")
            if self.verbose:
                print(f'---> array: {self.array}')
                print(f'---> project coode: {self.project_code}')
                print(f'---> PI: {self.pi}')
        except:
            if self.verbose:
                print(f"---> no telescope data found!")

    def _getBandData(self, fle):
        try:
            dset = fle['BAND']
            self.band_letter = dset[0].decode("ascii")
            self.rest_freq = float(dset[1].decode("ascii"))
            self.molecule = dset[2].decode("ascii")
            if self.verbose:
                print(f'---> band: {self.band_letter}')
                print(f'---> rest frequency: {self.rest_freq} MHz')
                print(f'---> molecule: {self.molecule}')
        except:
            if self.verbose:
                print(f"---> no project code found!")
    
    def _getSigmaData(self, fle):
        try:
            dset = fle['SIGMA']
            self.sigma_level = float(dset[0])
            if self.verbose:
                print(f'---> 1-sigma level in emission-free channel: {self.sigma_level} Jy/beam')
        except:
            if self.verbose:
                print(f"---> no sigma-level info found!")

    def _getOrigin(self, fle):
        dset = fle['ORIGIN']
        arr = np.array(dset)
        self.originRA = arr[0]
        self.originDEC = arr[1]
    
    def __makeFancyTicks(self, ax):
        ax.xaxis.set_tick_params(direction='in', width=1, length = 3, top = True, bottom=True)
        ax.xaxis.set_tick_params(direction='in', width=1, length = 3, which='minor', top = True, bottom=True)
        ax.yaxis.set_tick_params(direction='in', width=1, length = 3, right=True)
        ax.yaxis.set_tick_params(direction='in', width=1, length = 3, which='minor', right=True)

    '''
    PUBLIC:
    '''
    def plot(self):
        # --- figure --- 
        fig = plt.figure(figsize=(5.8,7))
        gs = gridspec.GridSpec(2, 1, height_ratios=[1,3])

        # --- map ---
        vmin = self.spectrum.velocity.min()
        vmax = self.spectrum.velocity.max()
        axMap = fig.add_subplot(gs[1,0])
        
        colors = self.spots.getJetColors(vmin = vmin, vmax = vmax)

        axMap.scatter(self.spots.dRA, self.spots.dDEC, s=np.log(self.spots.flux * 1000)**2.0 * 5, c=colors, edgecolor='black')
        axMap.invert_xaxis()
        axMap.set_xlabel('$\Delta$RA')
        axMap.set_ylabel('$\Delta$DEC')

        # -- spectrum ---
        axSpec = fig.add_subplot(gs[0,0])
        points = np.array([self.spectrum.velocity, self.spectrum.flux]).T.reshape(-1, 1, 2)
        segments = np.concatenate([points[:-1], points[1:]], axis=1)
        norm = plt.Normalize(vmin, vmax)
        lc = LineCollection(segments, cmap='jet', norm=norm)
        lc.set_array(self.spectrum.velocity)
        axSpec.add_collection(lc)
        axSpec.set_xlim(self.spectrum.velocity.min(), self.spectrum.velocity.max())
        dst = self.spectrum.flux.max() - self.spectrum.flux.min()
        axSpec.set_ylim(self.spectrum.flux.min() - 0.05 * dst, self.spectrum.flux.max() + 0.05 * dst)
        axSpec.set_xlabel("V$_{LSR}\,$(km$\,$s$^{-1}$)")
        axSpec.set_ylabel("Flux density (Jy)")

        # -- other --
        axSpec.tick_params(labelbottom=False,labeltop=True)
        axSpec.xaxis.set_label_position('top')

        # -- fancy ticks ---
        self.__makeFancyTicks(axMap)
        self.__makeFancyTicks(axSpec)

        plt.subplots_adjust(top=0.92, bottom=0.085, left=0.095, right=0.950, hspace=0.045)

        tmpPlot, = axMap.plot(np.nan, np.nan, label=self.project_code)
        axMap.legend(handles=[tmpPlot], loc="upper right", handlelength=0, handletextpad=0, framealpha = 0.9)

        plt.show()