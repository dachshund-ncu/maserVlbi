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

class cloudletClass:

    def __init__():
        pass

    def setAttributes(self, dRA, dRA_err, dDEC, dDEC_err, flux, flux_err, channels, velocity):
        self.spots = spotsClass()
        self.spots.set_dra(dRA, dRA_err)
        self.spots.set_ddec(dDEC, dDEC_err)
        self.spots.set_flux(flux, flux_err)
        self.spots.set_channels(channels)
        self.spots.set_velocity(velocity)

    
