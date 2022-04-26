'''
Class that holds info about cloudlets
cloudlet has some atrtributes:
-> dRA
-> dDEC
-> dRA_err
-> dDEC_err
-> meanFlux
-> meanFlux_err
-> maxFlux
-> maxFlux_err
And hosts instance of:
-> spots
in "spots" we hold all of the spots that are parts of this cloudlet
'''

from spotClass import spotsClass
import numpy as np
class cloudletClass:

    def __init__(self):
        pass

    def setAttributes(self, dRA, dRA_err, dDEC, dDEC_err, flux, flux_err, channels, velocity):
        '''
        Sets << spots >> attribute of the class. Other properties can be calculated from it later
        '''
        self.spots = spotsClass()
        self.spots.set_dra(dRA, dRA_err)
        self.spots.set_ddec(dDEC, dDEC_err)
        self.spots.set_flux(flux, flux_err)
        self.spots.set_channels(channels)
        self.spots.set_velocity(velocity)
    
    def loadFromHDF(self, group):
        self.spots = spotsClass()
        self.spots.readDataFromGroup(group)

    def calcProps(self):
        '''
        To be used only after 'SPOTS' attribute is created
        '''
        self.dRA = np.average(self.spots.dRA, weights=self.spots.flux)
        self.dDEC = np.average(self.spots.dDEC, weights=self.spots.flux)
        self.dRA_err = np.average(self.spots.dRA_err)
        self.dDEC_err = np.average(self.spots.dDEC_err)
        self.meanFlux = np.average(self.spots.flux)
        self.maxFlux = self.spots.flux.max()
        self.meanFlux_err = np.std(self.spots.flux)
        self.maxFlux_err = self.spots.flux[self.spots.flux.tolist().index(self.maxFlux)]

    def __str__(self):
        return f"x: {round(self.dRA, 2)}, y: {round(self.dDEC, 2)}, flux: {round(self.maxFlux, 2)}"
    
