###########################################################################
# This file holds the plotting routine for the KP solution from the Schmiesy
# Thesie. There is no timestepping, it plots only the initial condition.
############################################################################

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import axes3d
from scipy import interpolate
from plotKPSoln import *

# Choose the local filestructure based on what computer we are using. 
#
# Desktop in office.
localFileStructure = ('/home/jeremy/Documents/research/'
                        'RiemannSurfaces/jTEM-Jeremy/')
def main():

    # Loop over all examples
    for exampleNum in xrange(0,1):
        print exampleNum
        cFileName, sFileName, gFileName = defFileNames( exampleNum );
        coordData, solnData, groupData = loadData(cFileName, sFileName,
                gFileName)
        plotFilename = ( localFileStructure + 'jeremy/plotting/data/'
                'SingleTest/ExampleNum' + str(exampleNum) + '.eps')
        createPlot(coordData, solnData, groupData, plotFilename,
                exampleNum=str(exampleNum))

    # Show all the plots at the end.
    plt.show()
    

def defFileNames( exampleNum ):
    #coordfile = "./data/Test4/coords"+str(exampleNum)+".csv"
    coordfile = ( localFileStructure + 
                '/jeremy/plotting/data/SingleTest/coords' + str(exampleNum)+'.csv' )

    solsfile = ( localFileStructure +   
                '/jeremy/plotting/data/SingleTest/soln' + str(exampleNum)+'.csv' )

    groupfile = ( localFileStructure +   
                '/jeremy/plotting/data/SingleTest/group' + str(exampleNum)+'.csv' )

    return coordfile, solsfile, groupfile

if __name__ == '__main__':
    main()
