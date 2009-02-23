"""
This is an example script how to use the pybingings for PSInterface

Take a look at the init method defined at the bottom,
you need to set the photonics tables directory, then
run like this:
python -i pybindings_example.py
>>>> psi = init()
>>>> amp, probs = scanprob(psi)

Enjoy! Bernhard <bernhard.voigt@desy.de>
"""

try:
    from icecube import PSInterface
except ImportError, ie:
    print 'You need to build the pybindings of PSInterface to use this script'

import os
import random
import math

class OM(object):

    def __init__(self, x, y, z, orientation=-1):
        '''
        x, y, z in meter
        orientation -1 downwards
        '''
        self.x = x
        self.y = y
        self.z = z
        self.orientation = orientation


class Particle(object):

    def __init__(self, x, y, z, zenith, azimuth, length, energy, type_):
        '''
        x,y,z in meter
        zenith, azimuth in degree
        length in meter
        energy in GeV
        type_ 1=Electron 0=Muon
        '''
        self.x = x
        self.y = y
        self.z = z
        self.zenith = zenith
        self.azimuth = azimuth
        self.length = length
        self.energy = energy
        self.type_ = type_

        

        
class Photonics(object):
    """
    This class provides a conventient interface to the PSI_Photonics class
    of the C++ library PSInterface

    Convinient functions are added to extract
    photonics table data
    """

    def __setPSI__(self):
        self.psi = PSInterface.PSI_Photonics()

    def __init__(self, tablesPath='./', driverFilePath='./',
                 level1File='level1_table.list',
                 level2File='level2_table.list',
                 tableSelection=0,
                 interpol=7, verbosity=0,
                 angleLow=0., angleHigh=180.,
                 zLow=None, zHigh=None):
        

        """
        Returns a newly created Photonics object.

        tablesPath is the top level photonics dir
        driverFilePath is the location of file containing the table list
        level1/2File are the filenames o the files containing the table list
        tableSelection selects level1 (1) or level2 (2) or both (0) table sets
        interpol sets the interpolation mode
        verbosity set the verbosity level of the PSInerface
        angleLow/High reduces the loaded table set to the specfied angular range
        """

        # remember the current working directory
        cwd = os.getcwd()

        self.__setPSI__()
        self.tablePath_ = tablesPath
        self.psi.SetInterpolationMode(interpol)
        if tableSelection == 0 or tableSelection == 1:
            self.psi.SetPhotonicsVerbosity(1, verbosity)
            self.psi.SetAngularSelection(1, angleLow, angleHigh)
            if zLow is not None and zHigh is not None:
                self.psi.SetDepthSelection(1, zLow, zHigh)

            # change into the table directory
            try:
                os.chdir(tablesPath)
            except Exception, e:
                print "Can't change into table dir %s" % tablesPath
                os.chdir(cwd)
                raise e

            if not self.psi.LoadTables(1, driverFilePath + level1File, 0):
                os.chdir(cwd)
                raise Exception("Couldn't load level 1 tables")

        if tableSelection == 0 or tableSelection == 2:
            self.psi.SetPhotonicsVerbosity(2, verbosity)
            self.psi.SetAngularSelection(2, angleLow, angleHigh)
            if zLow is not None and zHigh is not None:
                self.psi.SetDepthSelection(2, zLow, zHigh)


            # change into the table directory
            try:
                os.chdir(tablesPath)
            except Exception, e:
                print "Can't change into table dir %s" % tablesPath
                os.chdir(cwd)
                raise e

            if not self.psi.LoadTables(2, driverFilePath + level2File, 0):
                os.chdir(cwd)
                raise Exception("Couldn't load level 2 tables")
        
        # back to the working dir
        os.chdir(cwd)


    def makeCoordinate(self, om, particle):
        """
        Returns a PSI_Coordinate object for the given OM and Particle objects
        """
        return self.psi.MakeCoordinate(om.x, om.y, om.z, om.orientation,
                                   particle.x, particle.y, particle.z, particle.zenith, particle.azimuth,
                                   particle.length, particle.energy, particle.type_)
    
    def getMeanNPE(self, coordinate):
        """
        Gets the mean number of photo electrons for the given PSI_Coordinate
        """
        return self.psi.GetMeanAmplitude(coordinate)


    def getDelayTimeProbability(self, delay, coordinate):
        """
        Returns the probability of the given delay time for the given PSI_Coordinate

        getMeanNPE must have been called previously for this OM and Paritcle combination!!!
        """
        return self.psi.GetProbability(float(delay), coordinate)


    def getTimeDelay(self, coordinate):
        """
        Returns a random time delay drawn from the arrival time distribution of photons
        emitted by a particle at the om given by the coordinate
        """
        return self.psi.GetTimeDelay(random.random(), coordinate)

# end class

# here's a short example

def init():

    # the photonics folder
    # set it to your needs
    tablesDir = '/afs/ifh.de/group/amanda/icecube/photonicstables/'
    # layerd ice tables
    driverFilePath = '/tables/wl06v210/I3Coord_I3Span/listfiles_wl06v210_z80_a60/'

    # load level 1 table for 0-20 deg zenith
    psi = Photonics(tablesPath = tablesDir,
                    driverFilePath = tablesDir + driverFilePath,
                    level1File = '/level1_shower.list',
                    tableSelection = 1,
                    angleLow=0, angleHigh=20)

    return psi

def scanprob(psi):

    # create a module and particle
    # at 0,0,0, zenith angle=10, azimuth= 120
    # length 0, eneryg 100TeV
    # and an electron
    p = Particle(0, 0, 0, 10, 120, 0, 100e3, 1)
    #om in 100 meter distance below
    om = OM(0, 0, -100)

    #make psi coordinate
    coo = psi.makeCoordinate(om, p)

    # get mean amplitude
    amp = psi.getMeanNPE(coo)

    # get delay time probability for delays from 0 to 2000 ns  in 3 ns steps
    probs = [psi.getDelayTimeProbability(delay, coo) for delay in xrange(0,2000,3)]

    return amp, probs



