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

    def __init__(self, filename):
        self._filename = filename
        self._fle = h5py.File(filename, 'r')
        # -- reading spots data --
        self.spots = spotsClass(self._fle['SPOTS'])
        # -- reading spec data --
        self.spectrum = spectrumClass(self._fle['SPECTRUM'])
        # -- optional: --
        self._getDate(self._fle)
        self._getBeam(self._fle)
        self._getProjectCode(self._fle)

    def _getDate(self, fle):
        try:
            dset = fle['DATE']
            times = [n.decode("ascii") for n in dset]
            times = Time(times, scale='utc', format='isot')
            self.mjd_begin = times[0].mjd
            self.mjd_end = times[1].mjd
            self.isot_begin = times[0].isot
            self.isot_end = times[1].isot
        except:
            print(f"---> \"{self._filename}\": no time data found!")
    
    def _getBeam(self, fle):
        try:
            dset = fle['BEAM']
            bmarray = np.array(dset)
            self.beam_majaxis = bmarray[0]
            self.beam_minaxis = bmarray[1]
            self.beam_posang = bmarray[2]
        except:
            print(f"---> \"{self._filename}\": no beam data found!")
    
    def _getProjectCode(self, fle):
        try:
            dset = fle['PROJECT_CODE']
            self.project_code = dset[0].decode("ascii")
        except:
            print(f"---> \"{self._filename}\": no project code found!")


