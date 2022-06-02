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
from scipy.optimize import curve_fit

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
        '''
        Reads << spots >> attribute from HDF group passed in the argument
        '''
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
        self.velocity = np.average(self.spots.velocity, weights=self.spots.flux)

    def __str__(self):
        return f"x: {round(self.dRA, 2)}, y: {round(self.dDEC, 2)}, flux: {round(self.maxFlux, 2)}"
    
    def shiftTo(self, shiftRA, shiftDEC):
        '''
        Shifts RA and DEC by a coords given in args
        '''
        self.dRA -= shiftRA
        self.dDEC -= shiftDEC
        self.spots.shiftTo(shiftRA, shiftDEC)
    
    '''
    =================
    GAUSSIAN ANALYSIS
    =================
    30.05.2022: STARTING DEVELOPING SINGULAR GAUSSIAN FITS
    31.05.2022: STARTING DEVELOPING MULTIPLE GAUSSIAN FITS
    '''

    def fit1DGaussianToData(self, verbose=False):
        '''
        Fits 1D gaussian function to << spots >> data
        '''
        # --- Starting parameters ---
        argtab = self.__determineInputArgs1D()
        bounds = self.__generateBounds1D(argtab)
        # --- returning ---
        return self.__fit1DGauss(argtab, bounds, verbose)
    
    def __fit1DGauss(self, argtab, bounds, verbose=False):
        coeff, varMatrix = curve_fit(self.__gaussian1D, self.spots.velocity, self.spots.flux, p0=argtab, bounds=bounds, method='trf')
        if verbose:
            self.__printFitResults(coeff, varMatrix)
        # --- returning ---
        return coeff, varMatrix

    def __gaussian1D(self, x, amp, vel, fwhm):
        return amp * np.exp(-(x-vel)**2.0 / (2.0 * (fwhm/2.35)**2.0))


    def __determineInputArgs1D(self):
        '''
        Determines imput args based on the stored data
        '''
        vel = np.mean(self.spots.velocity)
        amp = self.spots.flux.max()
        fwhm = 0.3
        return np.asarray([amp, vel, fwhm])
    
    def __generateBounds1D(self, argtab):
        '''
        Simply generates the tuple with bounds, that is easily accepted inside the
        SCIPY.OPTIMIZE.CURVE_FIT method
        Bounds asusmptions:
        - amplitude tolerance is 30%
        - maximum deviation from given in argument velocity is 0.6 km/s
        - maximum feature width is ~ 0.7 km/s
        - minumum is 0.05 km/s
        '''
        # -- tables --
        fluxDown = 0.0
        fluxUp = argtab[0] + 0.3 * argtab[0] # simply add 30% tolerance for amplitude
        # -
        velDown = argtab[1] - 0.3
        velUp = argtab[1] + 0.3 
        # - 
        fwhmDown =  0.05
        fwhmUp = 0.7
        
        return ( [fluxDown, velDown, fwhmDown], [fluxUp, velUp, fwhmUp] )
    
    def __printFitResults(self, coeff, varMatrix):
        '''
        prints fit results
        singular for now
        '''
        print("=============================")
        bl = np.sqrt(np.diag(varMatrix))
        for index, coef in enumerate(coeff):
            if (index) % 3 == 0 and index != 0:
                print("")
            print(f'{coef} +/- {bl[index]} ({bl[index] * 100.0 / coef})')

        print("=============================")
    '''
    ===============
    MULTIPLE GAUSS
    ===============
    '''

    def fitMultipleGauss(self, amp = [], vel = [], fwhm = [], bds = [], verbose=False):
        '''
        This is main public function, that fits gaussian function to the spots
        stored in this cloudlet class
        '''
        if amp == []: # YOU ALWAYS HAVE TO GIVE THE IMPUT PARAMETERS
            return
        
        if len(amp) != len(vel) or len(amp) != len(fwhm) or len(vel) != len(fwhm):
            # always these need to habve same length
            return
        
        if bds == []:
            bds = self.__generateBoundsXD(amp, vel, fwhm) # generates bounds from multiple entries

        # -- actual fitting --
        coeff, varMatrix = self.__fitMultiDimGauss(amp, vel, fwhm, bds)
        # --------------------
        if verbose:
            self.__printFitResults(coeff, varMatrix)
        return coeff, varMatrix
    
    def __fitMultiDimGauss(self, amp, vel, fwhm, bds):
        '''
        fits multi-dimensional gaussian to the data
        '''
        argtab = []
        argtab.extend(amp)
        argtab.extend(vel)
        argtab.extend(fwhm)
        return curve_fit(self.__multipleGaussWrapper, self.spots.velocity, self.spots.flux, p0=argtab, bounds=bds, maxfev=10000, sigma=self.spots.flux_err, method='trf')
    
    def __generateBoundsXD(self, amp, vel, fwhm):
        '''
        generates bounds for multiple gauss fittings
        '''
        try:
            length = len(amp)
        except:
            amp = [amp]
            vel = [vel]
            fwhm = [fwhm]
            length = 1
        down = []
        up = []

        amplitudes_down = []
        amplitudes_up = []

        arg_down = []
        arg_up = []

        fwhm_down = []
        fwhm_up = []

        for i in range(length):
            amplitudes_down.append(0)
            amplitudes_up.append(max(amp) + 0.3*max(amp))

            arg_down.append(vel[i]-0.3)
            arg_up.append(vel[i]+0.3)

            fwhm_down.append(0.05)
            fwhm_up.append(0.7)

        down.extend(amplitudes_down)
        down.extend(arg_down)
        down.extend(fwhm_down)

        up.extend(amplitudes_up)
        up.extend(arg_up)
        up.extend(fwhm_up)

        return (down, up)

    def __multipleGaussWrapper(self, x, *p):
        '''
        This is a wrapper, that splits the p argument table
        into three arrays, that can be easily accepted by self.__generalGauss()
        '''
        amp, arg_max, fwhm = self.__splitCoeffs(p)
        return self.__generalGauss(x,amp,arg_max,fwhm)

    def __splitCoeffs(self, p):
        '''
        This method simply splits coeffs into three arrays
        '''
        amplim = len(p) / 3.0
        arg_makslim = 2.0 * len(p) / 3.0
        amp = []
        arg_max = []
        fwhm = []
        for i in range(len(p)):
            if i < amplim:
                amp.append(p[i])
            elif i >= amplim and i < arg_makslim:
                arg_max.append(p[i])
            elif i >= arg_makslim:
                fwhm.append(p[i])
        return amp, arg_max, fwhm

    def __generalGauss(self, x, amp, arg_max, fwhm):
        '''
        this is a general gauss function
        it returns a sum of the multiple gaussian functions
        based on the length of <amp>
        It is assumed, that amp, arg_max and fwhm are the same lengrh
        '''
        suma = np.zeros(len(x))
        for fluxDensity, maximumVelicity, featureWidth in zip(amp, arg_max, fwhm):
            suma += self.__gaussian1D(x, fluxDensity, maximumVelicity, featureWidth)
        return suma

    def getTableToPlotGauss(self, coeffs, vmin = -1e11, vmax = 1e11):
        '''
        returns a table to easily plot gaussian function
        '''
        if vmin == -1e11:
            vmin = self.spots.velocity.min()
        if vmax == 1e11:
            vmax = self.spots.velocity.max() 
        x = np.linspace(vmin, vmax, 500)
        y = self.__multipleGaussWrapper(x, *coeffs)
        return x,y
    
    def getSingularGausses(self, coeffs, vmin = -1e11, vmax = 1e11):
        '''
        returns a tables of singular gausses
        '''
        if vmin == -1e11:
            vmin = self.spots.velocity.min()
        if vmax == 1e11:
            vmax = self.spots.velocity.max() 
        amp, v_max, fwhm = self.__splitCoeffs(coeffs)

        if len(amp) == 1:
            return [np.nan], [np.nan]

        x = []
        y = []
        for a,b,c in zip(amp, v_max, fwhm):
            print(a,b,c)
            xtmp = np.linspace(vmin, vmax, 500)
            ytmp = self.__gaussian1D(xtmp, a,b,c)
            x.append(xtmp)
            y.append(ytmp)
        return x,y
            
        
