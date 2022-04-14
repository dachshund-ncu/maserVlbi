'''
Class that holts spots data inside the maserVlbi object
'''

import numpy as np
import matplotlib.pyplot as plt

class spotsClass:
    def __init__(self, hdf5Group):
        self.__readDataFromGroup(hdf5Group)
    
    def __readDataFromGroup(self, hdf5Group):
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