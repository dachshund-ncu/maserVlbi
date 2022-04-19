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
        self.spots = spotsClass(self._fle['SPOTS'])
        # -- reading spec data --
        self.spectrum = spectrumClass(self._fle['SPECTRUM'])
        # -- optional: --
        self._getDate(self._fle)
        self._getBeam(self._fle)
        self._getProjectCode(self._fle)
        self._getBandData(self._fle)
        self._getSigmaData(self._fle)

        if self.verbose:
            print('-------------------------------')

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

