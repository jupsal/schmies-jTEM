###########################################################################
# This file holds routines for working with the KP Data from Java output. It is
# used for plotting KP solutions as well as other things.
############################################################################

import numpy as np

# Define filenames for data.
def defFileNames( globalFileStructure, localFileStructure, exampleNum ):
    # globalFileStructure defines what we are computer on essentially. It is the
    # home folder for jTEM stuff.
    # localFileStructure defines where the data for this test is.
    # exampleNum tells us which example we are looking at and is added onto
    # localFileStructure
    coordfile = ( globalFileStructure + localFileStructure +
                'coords' + str(exampleNum)+'.csv' )
    solsfile = ( globalFileStructure + localFileStructure +
                'soln' + str(exampleNum)+'.csv' )
    groupfile = ( globalFileStructure + localFileStructure +
                'group' + str(exampleNum)+'.csv' )
            
    return coordfile, solsfile, groupfile

def loadData( cFileName, sFileName, gFileName ):
    # Loads coordinate data, solution data, and group data from filenames
    coordData = np.genfromtxt(cFileName, dtype=float, delimiter=',', 
            names=True); # comes in (t,x,y) tuples
    solnData = np.genfromtxt(sFileName, dtype=float, delimiter=',', 
            names=True); # comes in (Real, Imaginary) pairs.
    groupData = np.genfromtxt(gFileName, dtype=float, delimiter=',', 
            names=True); # comes in (Real, Imaginary) pairs.
    return coordData, solnData, groupData

def parseGroupData( groupData ):
    # groupData is input from loadData. This parses it.
    # Should be called as ReCenters, ImCenters, Radii = parseGroupData( data );
    ReCenters = np.array( [ entry[0] for entry in groupData ] )
    ImCenters = np.array( [ entry[1] for entry in groupData ] )
    Radii = np.array( [ entry[2] for entry in groupData ] )
    return ReCenters, ImCenters, Radii

def parseSolnData( coordData, solnData ):
    # coordData and solnData are to be loaded in from loadData. This parses it.
    # It puts everything on the right grid with the right spacing. We require
    # both solution and coordinate data for this routine.

    # To get data on nxn grid, should be called as 
    #           x,y,z = parseSolnData( coordData, solnData )
    #           where z = function(x,y) which solves KP

    t = np.zeros(coordData.shape) # t vector
    x = np.zeros(coordData.shape) # x vector
    y = np.zeros(coordData.shape) # y vector
    j = 0;
    for txy in coordData:
        # Put the data in the right place.
        t[j] = txy[0]
        x[j] = txy[1]
        y[j] = txy[2]
        j = j+1;

    realSoln = np.array( [ sol[0] for sol in solnData ] )
    imagSoln = np.array( [ sol[1] for sol in solnData ] ) # if we want it.
    # Check that solution is entirely real-valued.
    if checkReSoln( imagSoln ) == False:
        raise ValueError("The solution has nonzero imaginary part.")

    tstep = len(t) # this ONLY works since we are only getting the solution for
                    # one time step

    # Find the xstep by checking for the value at which x changes from 0 to
    # nonzero. Since x changes last
    xstep = np.where(x[:-1] != x[1:])[0][0] + 1
    # Right now this only works for evenly spaced grids. I.e. xstep = ystep.
    ystep = xstep

    z = realSoln[ 0:tstep ]
    # Put it on a grid.
    xx = x.reshape( (xstep, xstep) ) # Put on nxn grid
    yy = y.reshape( (ystep, ystep) ) # Put on nxn grid
    zz = z.reshape( xx.shape )       # Put on same nxn grid

    return xx, yy, zz

def checkReSoln( imagSoln ):
    # Check to see that the solution is real-valued. We won't always want this
    # but sometimes we do.
    # If we want only real-valued solutions then you will want to throw this
    # after the check
    #    raise ValueError("The solution has nonzero imaginary part.")

    realValued = True;
    if sum( imagSoln )>0: 
        # then solution is not real-valued
        realValued = False

    return realValued



