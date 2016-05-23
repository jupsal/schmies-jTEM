###########################################################################
# This file holds the plotting routine for the KP solution from the Schmiesy
# Thesie. There is no timestepping, it plots only the initial condition.
############################################################################

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import axes3d
from scipy import interpolate


# Choose the local filestructure based on what computer we are using. 
#
# Desktop in office.
localFileStructure = ('/home/jeremy/Documents/research/'
                        'RiemannSurfaces/jTEM-Jeremy/')

def main():

    # Loop over all examples
    numExamples = 0
    for exampleNum in xrange(0,numExamples+1):
        print exampleNum
        cFileName, sFileName, gFileName = defFileNames( exampleNum );
        coordData, solnData, groupData = loadData(cFileName, sFileName,
                gFileName)
        plotFilename = ( localFileStructure + 'jeremy/plotting/data/'
                    'Test6/ExampleNum' + str(exampleNum) + '.eps')
        createPlot(coordData, solnData, groupData, plotFilename, 
                exampleNum = str(exampleNum))

    # Show all the plots at the end.
    plt.show()
    

def defFileNames( exampleNum ):
    coordfile = ( localFileStructure + 
                'jeremy/plotting/data/Test6/coords' + str(exampleNum)+'.csv' )

    solsfile = ( localFileStructure +   
                'jeremy/plotting/data/Test6/soln' + str(exampleNum)+'.csv' )

    groupfile = ( localFileStructure +
                'jeremy/plotting/data/Test6/group' + str(exampleNum)+'.csv' )

    return coordfile, solsfile, groupfile

def check1D( ):
    # This function checks to see if the data is in fact one dimensional.
    # Needs: x derivative of data, y derivative of data.
