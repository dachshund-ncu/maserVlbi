'''
Class that holts spots data inside the maserVlbi object
'''

import numpy as np
import matplotlib.pyplot as plt

class spotsClass:
    def __init__(self):
        pass
    def readDataFromGroup(self, hdf5Group):
        self.dRA = np.array(hdf5Group['RA'])
        self.dRA_err = np.array(hdf5Group['RA_ERR'])
        self.dDEC = np.array(hdf5Group['DEC'])
        self.dDEC_err = np.array(hdf5Group['DEC_ERR'])
        self.flux = np.array(hdf5Group['FLUX'])
        self.flux_err = np.array(hdf5Group['FLUX_ERR'])
        self.channels = np.array(hdf5Group['CHANNELS'])
        self.velocity = np.array(hdf5Group['VELOCITY'])
    
    def __scaledVel(self, vmin, vrange):
        return (self.velocity - vmin) / vrange
    
    def getJetColors(self, vmin = None, vmax = None):
        if vmin == None and vmax == None:
            vmin = self.velocity.min()
            vmax = self.velocity.max()
        vrange = abs(vmax - vmin)
        return plt.cm.jet(self.__scaledVel(vmin, vrange))

    def getBrightestSpot(self):
        spot_index = self.flux.tolist().index(self.flux.max())
        return self.dRA[spot_index], self.dRA_err[spot_index], self.dDEC[spot_index],  self.dDEC_err[spot_index], self.flux[spot_index], self.channels[spot_index], self.velocity[spot_index]
    
    def getSpotsFromRange(self, xrange, yrange):
        return_dRA = []
        return_dDEC = []
        return_flux = []
        return_vel = []
        for ra, dec, flux, vel in zip(self.dRA, self.dDEC, self.flux, self.velocity):
            if ra > xrange[0] and ra < xrange[1] and dec > yrange[0] and dec < yrange[1]:
                return_dRA.append(ra)
                return_dDEC.append(dec)
                return_flux.append(flux)
                return_vel.append(vel)
        return np.asarray(return_dRA), np.asarray(return_dDEC), np.asarray(return_flux), np.asarray(return_vel)
    
    def getPropsFromIndex(self, index):
        return self.dRA[index], self.dDEC[index], self.flux[index], self.velocity[index]
    
    '''
    INITIALIZING THIS CLASS FROM NO HDF5:
    '''
    def set_dra(self, dra, dra_err):
        self.dRA = dra
        self.dRA_err = dra_err
    def set_ddec(self, ddec, ddec_err):
        self.dDEC = ddec
        self.dDEC_err = ddec_err
    def set_flux(self, f, f_er):
        self.flux= f
        self.flux_err = f_er
    def set_channels(self, ch):
        self.channels = ch
    def set_velocity(self, vel):
        self.velocity = vel

    '''
    SPOT SHIFTING
    '''
    def shiftTo(self, shiftRA, shiftDEC):
        '''
        Shifts RA and DEC by a coords, given in args
        '''
        self.dRA -= shiftRA
        self.dDEC -= shiftDEC
    
    def __str__(self):
        '''
        what is printed, when 'print()' is called?
        '''
        printstr = ""
        for index, ra in enumerate(self.dRA):
            printstr += f"{int(self.channels[index])} "
            printstr += f"{str(round(self.velocity[index],2)).zfill(2) } "
            printstr += f"{round(self.flux[index],3)} "
            printstr += f"{round(self.flux_err[index],3)} "
            printstr += f"{round(ra,2)} "
            printstr += f"{round(self.dRA_err[index],2)} "
            printstr += f"{round(self.dDEC[index],2)} "
            printstr += f"{round(self.dDEC_err[index],2)}\n"
        return printstr