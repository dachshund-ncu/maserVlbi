'''
Class that holts spots data inside the maserVlbi object
'''

import numpy as np

class spectrumClass:
    def __init__(self, hdf5Group):
        self.__readDataFromGroup(hdf5Group)
    
    def __readDataFromGroup(self, hdf5Group):
        self.channels = np.array(hdf5Group['CHANNELS'])
        self.velocity = np.array(hdf5Group['VELOCITY'])
        self.flux = np.array(hdf5Group['FLUX'])

    def getJetColors():
        pass

    def rotate(self, channels):
        deltaVel = abs(self.velocity[1] - self.velocity[0])
        self.velocity += deltaVel*channels