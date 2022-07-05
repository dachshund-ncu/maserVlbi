#! /usr/bin/env python3
# -*- coding: utf-8 -*-
'''
This is main class of the maser data
Mainly it reads the HDF5 files
'''

from asyncore import read
import h5py
from spotClass import spotsClass
from spectrumClass import spectrumClass
from cloudletClass import cloudletClass
from astropy.time import Time
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.collections import LineCollection

class maserVlbi:

    def __init__(self, filename, verbose=False):
        self.verbose = verbose
        self._filename = filename

        if self.verbose:
            print('-------------------------------')
            print(f'---> Loading file \"{self._filename}\"...')

        self._fle = h5py.File(filename, 'r+')
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
        self._getCloudlets(self._fle)
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
    
    def _getCloudlets(self, fle):
        self.cloudlets = []
        try:
            dset = fle['CLOUDLETS']
            for clGroup in dset:
                cloud = cloudletClass()
                cloud.loadFromHDF(dset[clGroup]['SPOTS'])
                cloud.calcProps()
                self.cloudlets.append(cloud)
            print(f'---> Succesfully loaded {len(self.cloudlets)} cloudlets!')
        except:
            if (self.verbose):
                print(f"---> no cloudlet information found!")
        
    def __makeFancyTicks(self, ax):
        ax.xaxis.set_tick_params(direction='in', width=1, length = 3, top = True, bottom=True)
        ax.xaxis.set_tick_params(direction='in', width=1, length = 3, which='minor', top = True, bottom=True)
        ax.yaxis.set_tick_params(direction='in', width=1, length = 3, right=True)
        ax.yaxis.set_tick_params(direction='in', width=1, length = 3, which='minor', right=True)
    
    '''
    OTHER:
    '''
    def appendCloudlet(self, cloudletObj):
        self.cloudlets.append(cloudletObj)
    
    def removeCloudlet(self, index):
        self.cloudlets.pop(index)
    
    def printCloudlets(self):
        for cl in self.cloudlets:
            print(cl)

    def saveCloudlets(self):
        '''
        Method to save information about Cloudlets to the opened file
        -> If there is no cloudlets, this method should do sh*t
        -> If there is 'CLOUDLET' key in file, it's contents should be wiped before proceeding
        further
        -> At every pass new 'CLOUDLETS' key should be created
        '''
        if len(self.cloudlets) < 1:
            return # as I said...
        # -
        if 'CLOUDLETS' in self._fle.keys():
            del self._fle['CLOUDLETS']
        cloudletGroup = self._fle.create_group('CLOUDLETS') # <--- ALWAYS CREATE NEW GROUP ON SAVE
        # -
        flag = True
        for index, cloudlet in enumerate(self.cloudlets):
            try:
                oneCl = cloudletGroup.create_group("CLOUDLET#" + str(index))
                self.__saveSpotsToGroups(oneCl, cloudlet)
            except:
                flag = False
                print(f"---> Error writing cloudlet no. {str(index)}")
        if flag:
            print(f'---> Succesfully saved {len(self.cloudlets)} cloudlets!')
        else:
            print(f'---> For some reason, saving cloudlets was not fully succesfull')
    
    def clearCloudletInfo(self):
        if 'CLOUDLETS' in self._fle.keys():
            del self._fle['CLOUDLETS']
            self.cloudlets = []

    def __saveSpotsToGroups(self, group, cloudlet):
        '''
        This method saves information of singular cloudlet
        '''
        spotsGr = group.create_group('SPOTS')
        spotsGr.create_dataset("RA", data = cloudlet.spots.dRA)
        spotsGr.create_dataset("RA_ERR", data = cloudlet.spots.dRA_err)
        spotsGr.create_dataset("DEC", data = cloudlet.spots.dDEC)
        spotsGr.create_dataset("DEC_ERR", data = cloudlet.spots.dDEC_err)
        spotsGr.create_dataset("FLUX", data = cloudlet.spots.flux)
        spotsGr.create_dataset("FLUX_ERR", data = cloudlet.spots.flux_err)
        spotsGr.create_dataset("CHANNELS", data = cloudlet.spots.channels)
        spotsGr.create_dataset("VELOCITY", data = cloudlet.spots.velocity)

    def getClJetColors(self, vMin, vMax):
        '''
        Simply gets cloudlet jet colors
        '''
        if len(self.cloudlets) < 1:
            return False
        veltab = []
        for cloudlet in self.cloudlets:
            veltab.append(cloudlet.velocity)
        veltab = np.asarray(veltab)
        vRange = abs(vMax - vMin)
        scaledVeltab = (veltab - vMin) / vRange
        return plt.cm.jet(scaledVeltab)

    def saveSpotsToAscii(self, filename):
        '''
        This method simply saves spot information into the ASCII file
        '''
        fle = open(filename, 'w+')
        fle.write('# ----------------------------\n')
        fle.write(f'# TIME: {self.isot_begin}\n')
        fle.write(f'# BEAM-SIZE: {self.beam_raAxis} x {self.beam_decAxis} ({self.beam_posang})\n')
        fle.write(f'# ORIGIN: {self.originRA} {self.originDEC}\n')
        fle.write(f'# PROJECT CODE: {self.project_code}\n')
        fle.write(f'# PI: {self.pi}\n')
        fle.write('# ----------------------------\n')
        fle.write('# Channel | Velocity (km/s) | FLux Density (Jy/beam) | Flux error | dRA (mas) | dRA error (mas) | dDEC (mas) | dDEC error (mas)\n')
        for i, vel in enumerate(self.spots.velocity):
            fle.write('%d   %.3f   %.4f   %.4f   %.3f   %.3f   %.3f   %.3f \n' % (self.spots.channels[i], vel, self.spots.flux[i], self.spots.flux_err[i], self.spots.dRA[i], self.spots.dRA_err[i], self.spots.dDEC[i], self.spots.dDEC_err[i]))
        fle.close()

    def plot(self, vmin = None, vmax = None):
        # --- figure --- 
        fig = plt.figure(figsize=(5.8,7))
        gs = gridspec.GridSpec(2, 1, height_ratios=[1,3])

        # --- map ---
        if vmin is None and vmax is None:
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

    def shiftToSpot(self, spotIndex):
        '''
        Shifts spots, cloudlets and changes the origin of the plot:
        '''
        # --- get spot coords (MAS) ---
        shiftRA, shiftDEC = self.__getShiftValue(spotIndex)
        self.spots.shiftTo(shiftRA, shiftDEC)
        [cloudlet.shiftTo(shiftRA, shiftDEC) for cloudlet in self.cloudlets]
        self.__shiftOriginTo(shiftRA, shiftDEC)

    def __getShiftValue(self, spotIndex):
        '''
        Shifts all of the spots to the spot of the given index
        It also changes the Origin of the plot
        '''
        try:
            shiftRA = self.spots.dRA[spotIndex]
            shiftDEC = self.spots.dDEC[spotIndex]
            return shiftRA, shiftDEC
        except:
            raise Exception
        
    def __shiftOriginTo(self, shiftRA, shiftDEC):
        '''
        Shifts the origin by the shiftRA and shiftDEC values
        DEC is easy - we only need to convert shiftDEC from MAS do DEGREES
        RA is harder - we need to compensate for cos(DEC)
        RA is in HR
        '''
        # DEC
        self.originDEC -= (shiftDEC / 3600.0 / 1000.0)
        # RA
        realRaShift = shiftRA / 3600.0 / 1000.0 # mas -> deg
        realRaShift /= 15 # def -> hr
        realRaShift *= 1.0 / np.cos(np.radians(self.originDEC)) # RA shift, compensated for DEC
        self.originRA -= realRaShift
    
    def printOrigin(self):
        '''
        Basically just prints origin in easy to read form
        '''
        print(f"RA: {self.__getRaStr()}, DEC: {self.__getDecStr()}")
    
    def __getRaStr(self):
        '''
        Gets nice-looking string from OriginRa
        '''
        # RA
        RA_hr = int(self.originRA)
        RA_min = (self.originRA % 1) * 60
        RA_sec = (RA_min % 1)*60.0
        RA_min = int(RA_min)
        #RaStr = str(RA_hr).zfill(2) + ":" + str(RA_min).zfill(2) + ":" + str(round(RA_sec,2)).zfill(5)
        return str(RA_hr).zfill(2) + ":" + str(RA_min).zfill(2) + ":" + ('%02.3f' % (round(RA_sec,3))).rjust(6, "0")
    
    def __getDecStr(self):
        '''
        Gets nice-looking string from OriginDEC
        '''
        # DEC
        if self.originDEC < 0:
            decSign = '-'
        else:
            decSign = ''
        tmpdec = abs(self.originDEC)
        DEC_degr = int(tmpdec)
        DEC_min = (tmpdec %1) * 60.0
        DEC_sec = (DEC_min % 1) * 60.0
        DEC_min = int(DEC_min)

        return decSign + str(DEC_degr).zfill(2) + ":" + str(DEC_min).zfill(2) + ":" + ('%02.3f' % (round(DEC_sec,3))).rjust(6, "0")