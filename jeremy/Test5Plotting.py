###########################################################################
# This file holds the plotting routine for the KP solution from the Schmiesy
# Thesie. There is no timestepping, it plots only the initial condition.
############################################################################

from plotting.plotKPSoln import createPlot
from parseKPData import *
import matplotlib.pyplot as plt
import os


# Choose the local filestructure based on what computer we are using. This
# should point to the main folder for jTEM stuff/jeremy (because of the way the
# following command works it's easier to do this)
#
globalFileStructure = os.path.dirname(os.path.realpath(__file__))
globalFileStructure += '/'
# For example, for the desktop in office, this should give
#    #  globalFileStructure = ('/home/jeremy/Documents/research/'
#    #                    'RiemannSurfaces/jTEM-Jeremy/jeremy/')



## We should only have to change the next two parameters for each test.

# localFileStructure tells us where the data is
localFileStructure = ('plotting/data/Test5/')
# Number of examples
numExamples = 6

def main():

    # Loop over all examples
    for exampleNum in xrange(0,numExamples+1):
        print exampleNum # to keep track of progress
        cFileName, sFileName, gFileName = defFileNames( globalFileStructure, 
                localFileStructure, exampleNum );
        coordData, solnData, groupData = loadData( cFileName, sFileName,
                gFileName )
        plotFilename = ( globalFileStructure + localFileStructure +
                        'ExampleNum' + str(exampleNum) + '.eps')
        createPlot( coordData, solnData, groupData, plotFilename, 
                exampleNum = str(exampleNum) )

    # Show all the plots at the end.
    plt.show()
    

if __name__ == '__main__':
    main()
