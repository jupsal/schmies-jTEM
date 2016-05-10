###########################################################################
# This file holds the plotting routine for the KP solution from the Schmiesy
# Thesie. There is no timestepping, it plots only the initial condition.
############################################################################

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import axes3d
from scipy import interpolate


def main():

    # Loop over all examples
    for exampleNum in xrange(0,3+1):
        print exampleNum
        cFileName, sFileName, gFileName = defFileNames( exampleNum );
        coordData, solnData, groupData = loadData(cFileName, sFileName,
                gFileName)
        createPlot(coordData, solnData, groupData, str(exampleNum))

    # Show all the plots at the end.
    plt.show()
    

def defFileNames( exampleNum ):
    coordfile = ('/home/jeremy/Documents/research/RiemannSurfaces/jTEM-Jeremy/'
                '/jeremy/plotting/data/Test4/coords' + str(exampleNum)+'.csv' )

    solsfile = ('/home/jeremy/Documents/research/RiemannSurfaces/jTEM-Jeremy/'  
                '/jeremy/plotting/data/Test4/soln' + str(exampleNum)+'.csv' )

    groupfile = ('/home/jeremy/Documents/research/RiemannSurfaces/jTEM-Jeremy/'  
                '/jeremy/plotting/data/Test4/group' + str(exampleNum)+'.csv' )

    return coordfile, solsfile, groupfile


def loadData( cFileName, sFileName, gFileName ):
    coordData = np.genfromtxt(cFileName, dtype=float, delimiter=',', 
            names=True); # comes in (t,x,y) tuples
    solnData = np.genfromtxt(sFileName, dtype=float, delimiter=',', 
            names=True); # comes in (Real, Imaginary) pairs.
    groupData = np.genfromtxt(gFileName, dtype=float, delimiter=',', 
            names=True); # comes in (Real, Imaginary) pairs.
    return coordData, solnData, groupData

def createPlot(coordData, solnData, groupData, exampleNum=' NOT GIVEN'):
    # Create both the KP Soln plot and the Group Data plot side-by-side
    fig = plt.figure()

    # First plot the group data
    ax1 = plt.subplot(2,1,1)
    plotGroupData( groupData, ax1 )

    # Then plot the KP Soln
    ax2 = plt.subplot(2,1,2, projection='3d')
    plotKP( coordData, solnData, ax2 )
    fig.suptitle('Example Number'+str(exampleNum))
    fig.savefig( '/home/jeremy/Documents/research/RiemannSurfaces/jTEM-Jeremy/jeremy/plotting/data/Test4/ExampleNum'
            + str(exampleNum) + '.eps', format = 'eps' )

    
def plotGroupData( groupData, ax ):
    ReCenters = np.array( [ entry[0] for entry in groupData ] )
    ImCenters = np.array( [ entry[1] for entry in groupData ] )
    Radii = np.array( [ entry[2] for entry in groupData ] )
    ReC0 = lambda t,j: ReCenters[j]+Radii[j]*np.cos(t)
    ImC0 = lambda t,j: ImCenters[j]+Radii[j]*np.sin(t)
    teval = np.arange( 0.0, 2*np.pi, 0.1 ) #the range to evaluate over
    ax.plot( ReC0(teval,0), ImC0(teval,0 ) )
    ax.plot( ReC0(teval,1), ImC0(teval,1 ) )
    ax.set_xlim( [-6, 6] )
    ax.set_ylim( [-6, 6] )
    plt.gca().set_aspect('equal')

def plotKP( coordData, solnData, ax ):
    #Define the data
    t = np.zeros(coordData.shape)
    x = np.zeros(coordData.shape)
    y = np.zeros(coordData.shape)
    j = 0;
    for txy in coordData:
        t[j] = txy[0]
        x[j] = txy[1]
        y[j] = txy[2]
        j = j+1;

    # Check that solution is entirely real-valued.
    if sum( [ sol[1] for sol in solnData ] )>0: # then solution is not real-valued
        raise ValueError("The solution has nonzero imaginary part.")

    realSoln = np.array( [ sol[0] for sol in solnData ] )

    tstep = len(t)
    # Find the xstep by checking for the value at which x changes from 0 to
    # nonzero. Since x changes last
    xstep = np.where(x[:-1] != x[1:])[0][0] + 1
    # Right now this only works for evenly spaced grids. I.e. xstep = ystep.
    ystep = xstep

    z = realSoln[ 0:tstep ]
    # Put it on a grid.
    xx = x.reshape( (xstep, xstep) )
    yy = y.reshape( (ystep, ystep) )
    zz = z.reshape( xx.shape )

    surf = ax.plot_surface( xx, yy, zz, rstride=2, cstride=2, cmap=cm.coolwarm,
            linewidth = 0.2)#, antialiased = True )
    ax.set_zlim( np.min(realSoln), np.max(realSoln) )
    ax.set_xlabel('x'); ax.set_ylabel('y');
    #ax.set_title('Example Number'+exampleNum)
    


if __name__ == '__main__':
    main()
