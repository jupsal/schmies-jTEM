###########################################################################
# This file holds the plotting routine for the KP solution from the Schmiesy
# Thesie. There is no timestepping, it plots only the initial condition.
############################################################################

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import axes3d
from scipy import interpolate
from parseKPData import *


def createPlot( x, y, z, groupData, 
        figureLocation, exampleNum=' NOT GIVEN'):

    # Create both the KP Soln plot and the Group Data plot side-by-side
    fig = plt.figure()

    # First plot the group data
    ax1 = plt.subplot(2,2,1)
    plotGroupData( groupData, ax1 )

    ax2 = plt.subplot(2,2,2)
    plotGroupTable( groupData, ax2 )

    # Then plot the KP Soln
    ax3 = plt.subplot(2,2,3, projection='3d')
    ax4 = plt.subplot(2,2,4, projection='3d')
    plotKP( x, y, z, ax3, ax4 )
    fig.suptitle('Example Number'+str(exampleNum))
    fig.savefig( figureLocation, format='eps' )
#    fig.savefig( '/home/jeremy/Documents/research/RiemannSurfaces/jTEM-Jeremy/jeremy/plotting/data/Test4/ExampleNum'
 #           + str(exampleNum) + '.eps', format = 'eps' )
    

def plotGroupData( groupData, ax ):

    ReCenters, ImCenters, Radii = parseGroupData( groupData );

    ReC0 = lambda t,j: ReCenters[j]+Radii[j]*np.cos(t)
    ImC0 = lambda t,j: ImCenters[j]+Radii[j]*np.sin(t)
    teval = np.arange( 0.0, 2*np.pi, 0.1 ) #the range to evaluate over
    ax.plot( ReC0(teval,0), ImC0(teval,0 ) )
    ax.plot( ReC0(teval,1), ImC0(teval,1 ) )
    ax.set_xlim( [-6, 6] )
    ax.set_ylim( [-6, 6] )
    plt.gca().set_aspect('equal')

def plotGroupTable( groupData, ax ):
    colLabels=("Real-Center","Imag-Center","Radius")
    cellText=[]
    for row in groupData:
        cellText.append(row)

    the_table = ax.table( cellText = cellText,
            colLabels = colLabels,
            loc = 'center' )


def plotKP( x, y, z, ax1, ax2 ):
    # Plot with two different angles.
    surf = ax1.plot_surface( x, y, z, rstride=2, cstride=2, cmap=cm.coolwarm,
            linewidth = 0.2)#, antialiased = True )
    #ax1.set_zlim( np.min(realSoln), np.max(realSoln) )
    ax1.set_zlim( np.min(z), np.max(z) )
    ax1.set_xlabel('x'); ax1.set_ylabel('y');
    ax1.view_init(azim = 0, elev = 90)

    surf = ax2.plot_surface( x, y, z, rstride=2, cstride=2, cmap=cm.coolwarm,
            linewidth = 0.2)#, antialiased = True )
    #ax2.set_zlim( np.min(realSoln), np.max(realSoln) )
    ax2.set_zlim( np.min(z), np.max(z) )
    ax2.set_xlabel('x'); ax2.set_ylabel('y');

