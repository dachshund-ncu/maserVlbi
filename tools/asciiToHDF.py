#! /usr/bin/env python3
# -*- coding: utf-8 -*-
'''
Simple tool that converts ascii data into fully self-contained HDF files
We need two files:
1. JMFIT file from AIPS
2. ISPEC file from AIPS (with velocities included)
'''

import numpy as np
import sys
import h5py
import configparser

class jmfitFile:
    def __init__(self, filename):
        self.loadJmfitFile(filename)
    
    def loadJmfitFile(self, filename):
        '''
        This is methood for reading .jmfit file, that is output of the aips task jmfit
        We need to focus on few things:
        1 - get most needed info (dRA, dDEC, flux, errors, channels)
        2 - get information about the origin (RA, DEC)
        TODO:
        3 - get other information (position angle, details on the size)
        '''
        # -- reading file --
        fle = open(filename, 'r+')
        a = fle.readlines()
        fle.close()
        # -- getting RA, DEC, FLUX and CHANNEL --
        flux = []
        flux_err = []
        x = []
        x_err = []
        y = []
        y_err = []
        channel = []

        for i in range(len(a)):
            if a[i][0] == '#':
                continue
            
            tmp = a[i].split()

            if tmp[5] == 'm':
                flux.append(float(tmp[3]) / 1000.0)
                flux_err.append(float(tmp[4]) / 1000.0)
                x.append(float(tmp[6]))
                x_err.append(float(tmp[7]))
                y.append(float(tmp[8]))
                y_err.append(float(tmp[9]))
            else:
                flux.append(float(tmp[3]))
                flux_err.append(float(tmp[4]))
                x.append(float(tmp[5]))
                x_err.append(float(tmp[6]))
                y.append(float(tmp[7]))
                y_err.append(float(tmp[8]))
            channel.append(int(tmp[2]))
        
        self.flux = np.asarray(flux)
        self.flux_err = np.asarray(flux_err)
        self.dRA = np.asarray(x)
        self.dRA_err = np.asarray(x_err)
        self.dDEC = np.asarray(y)
        self.dDEC_err = np.asarray(y_err)
        self.chans = np.asarray(channel)

        # -- getting info about DEC and RA --
        if tmp[5] == 'm':
            raDECindex = 10
        else:
            raDECindex = 9
        
        raString = tmp[raDECindex]
        decString = tmp[raDECindex + 1]

        rah = float(raString[0:2])
        ram = float(raString[2:4])
        ras = float(raString[4:])
        self.originRA = rah + ram / 60.0 + ras / 3600.0
        
        if float(decString) >= 0:
            decd = float(decString[0:2])
            decm = float(decString[2:4])
            decs = float(decString[4:])
            self.originDEC = decd + decm / 60.0 + decs / 3600.0
        else:
            decd = float(decString[0:3])
            decm = float(decString[3:5])
            decs = float(decString[5:])
            self.originDEC = -1.0 * (abs(decd) + abs(decm) / 60.0 + abs(decs) / 3600.0)

    def calibrateVelocities(self, ispecFile):
        '''
        We simply take velocities and the corresponding channels from the ISPEC file and
        put them into an array here
        '''
        channelIndexes = [ispecFile.channels.tolist().index(channelNumber) for channelNumber in self.chans]
        self.vel = ispecFile.vlsr[channelIndexes]
        #[print(self.vel[i], self.chans[i]) for i in range(len(self.vel))]

class ispecFile():
    def __init__(self, filename):
        self.loadIspecFile(filename)
    
    def loadIspecFile(self, filename):
        '''
        ISPEC files are simple:
        channel velocity flux_density
        so they can be loaded through numpy.LOADTXT
        '''
        self.channels, self.vlsr, self.fluxDensity = np.loadtxt(filename, unpack=True, usecols=(0,1,2))
        self.vlsr = self.vlsr / 1000.0 # calibrating into km/s

def makeHDF5File(filename, jmfit_file, ispec_file, configFile = None):
    fle = h5py.File(filename, 'w')
    # -- adding SPOTS info --
    # - group
    spotsGroup = fle.create_group("SPOTS")
    # - spot data
    spotsRA = spotsGroup.create_dataset("RA", data = jmfit_file.dRA)
    spotsRA_ERR = spotsGroup.create_dataset("RA_ERR", data = jmfit_file.dRA_err)
    spotsDEC = spotsGroup.create_dataset("DEC", data = jmfit_file.dDEC)
    spotsDEC_ERR = spotsGroup.create_dataset("DEC_ERR", data = jmfit_file.dDEC_err)
    spots_FLUX = spotsGroup.create_dataset("FLUX", data = jmfit_file.flux)
    spots_FLUX_ERR = spotsGroup.create_dataset("FLUX_ERR", data = jmfit_file.flux_err)
    spots_CHANNELS = spotsGroup.create_dataset("CHANNELS", data = jmfit_file.chans)
    spots_VELOCITY = spotsGroup.create_dataset("VELOCITY", data = jmfit_file.vel)
    # - specrum data
    spectrumGroup = fle.create_group("SPECTRUM")
    specChannels = spectrumGroup.create_dataset("CHANNELS", data = ispec_file.channels)
    specVelocity = spectrumGroup.create_dataset("VELOCITY", data = ispec_file.vlsr)
    specFlux = spectrumGroup.create_dataset("FLUX", data = ispec_file.fluxDensity)
    # - other datasets
    origin = fle.create_dataset("ORIGIN", data = np.array([jmfit_file.originRA, jmfit_file.originDEC]))

    if configFile == None:
        if askForPermission("Do you want to include beam info? (y/*)"):
            beam_majaxis, beam_minaxis, pos_ang = getBeamInfo()
            addBeamInfo(fle, [beam_majaxis, beam_minaxis, pos_ang])
        
        if askForPermission("Do you want to include time info? (y/*)"):
            isotBegin, isotEnd = getTimeInfo()
            addDateInfo(fle, [isotBegin, isotEnd])
        
        if askForPermission("Do you want to include project code? (y/*)"):
            telescope = getTelescopeInfo()
            addTelescopeInfo(fle, telescope)
    else:
        beam, time, band, telescope = readConfigFile(configFile)
        addBeamInfo(fle, beam)
        addDateInfo(fle, time)
        addBandInfo(fle, band)
        addTelescopeInfo(fle, telescope)
    

def askForPermission(question):
    print(question)
    fl = input()
    if fl.lower() == 'y':
        return True
    else:
        return False

def getTimeInfo():
    flag = True
    while flag:
        print("Provide isotime for begin of the project: ")
        try:
            time_beg = input('--> ')
            flag = False
        except:
            print("Error! Try again")
    flag = True
    while flag:
        print("Provide isotime for end of the project: ")
        try:
            time_end = input('--> ')
            flag = False
        except:
            print("Error! Try again")
    return time_beg, time_end

def getBeamInfo():
    flag = True
    while flag:
        print("Provide beam major axis: ")
        try:
            beam_majaxis = float(input('--> '))
            flag = False
        except:
            print("Error! Try again")
    flag = True
    while flag:
        print("Provide beam minor axis: ")
        try:
            beam_minaxis = float(input('--> '))
            flag = False
        except:
            print("Error! Try again")
    
    flag = True
    while flag:
        print("Provide position angle: ")
        try:
            posang = float(input('--> '))
            flag = False
        except:
            print("Error! Try again")
    
    return beam_majaxis, beam_minaxis, posang

def getTelescopeInfo():
    flag = True
    while flag:
        print("Write array name: ")
    try:
        array = input('--> ')
        flag = False
    except:
        print("Error! Try again")
    flag = True
    while flag:
        print("Write project code: ")
        try:
            code = input('--> ')
            flag = False
        except:
            print("Error! Try again")
    flag = True
    while flag:
        print("Write PI name: ")
        try:
            pi = input('--> ')
            flag = False
        except:
            print("Error! Try again")
    return [array, code, pi]

def addBandInfo(fle, band):
    strList = [n.encode("ascii", "ignore") for n in band]
    fle.create_dataset("BAND", data=strList)

def addTelescopeInfo(fle, telescope):
    strList = [code.encode("ascii", "ignore") for code in telescope]
    fle.create_dataset("TELESCOPE", data=strList)

def addBeamInfo(fle, beam):
    fle.create_dataset("BEAM", data = np.array(beam) )


def readConfigFile(configFileName):
    confile = configparser.ConfigParser()
    confile.read(configFileName)
    beam = []
    time = []
    band = []
    telescope = []
    if 'BEAM' in confile.sections():
        beam.append(float(confile['BEAM']['beam_majaxis']))
        beam.append(float(confile['BEAM']['beam_minaxis']))
        beam.append(float(confile['BEAM']['beam_posang']))
    if 'TIME' in confile.sections():
        time.append(confile['TIME']['isot_begin'])
        time.append(confile['TIME']['isor_end'])
    if 'BAND' in confile.sections():
        band.append(confile['BAND']['band_letter'])
        band.append(confile['BAND']['rest_freq'])
        band.append(confile['BAND']['molecule'])
    if 'TELESCOPE' in confile.sections():
        telescope.append(confile['TELESCOPE']['array'])
        telescope.append(confile['TELESCOPE']['code'])
        telescope.append(confile['TELESCOPE']['pi'])
    return beam, time, band, telescope



def readConfigJmfitandIspec(configFileName):
    print(configFileName)
    confile = configparser.ConfigParser()
    confile.read(configFileName)
    if 'FILES' in confile.sections():
        jmfitFileName = confile['FILES']['jmfit_file']
        ispecFileName = confile['FILES']['ispec_file']
        jmfit = jmfitFile(jmfitFileName)
        ispec = ispecFile(ispecFileName)
        jmfit.calibrateVelocities(ispec)
        return jmfit, ispec
    else:
        return None



def addDateInfo(fle, date):
    #print(dateIsotBegin)
    #print(dateIsotEnd)
    #wppl = np.array([dateIsotBegin, dateIsotEnd], dtype=str)
    strList = [n.encode("ascii", "ignore") for n in date]
    #wppl = np.void(wppl)
    fle.create_dataset("DATE", data = strList)
    #fle.attrs['DATE_BEGIN'] = dateIsotBegin
    #fle.attrs['DATE_END'] = dateIsotEnd


if __name__ == '__main__':

    if '-o' in sys.argv:
        index = sys.argv.index('-o')
        output_filename = sys.argv[index + 1]
    else:
        output_filename = 'out.hdf5'
    
    if '-conf' in sys.argv:
        index = sys.argv.index('-conf')
        config_file = sys.argv[index + 1]
    else:
        config_file = None
    
    if config_file == None:
        eee = jmfitFile(sys.argv[1])
        eespec = ispecFile(sys.argv[2])
        eee.calibrateVelocities(eespec)
    else:
        eee, eespec = readConfigJmfitandIspec(config_file)

    makeHDF5File(output_filename, eee, eespec, config_file)