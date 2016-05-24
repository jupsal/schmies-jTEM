###########################################################################
# This file just takes a derivative of the data
############################################################################

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d
from matplotlib import cm


def checkStep2( x, y, z, plot=False, plotTitle = 'u_y over u_x' ):
    # Step 2 says that if u_y(x,y,0) = \beta u_x(x,y,0) for some \beta, then the
    # data is one-dimensional. We calculate u_y and u_x then divide pointwise to
    # see if there is a pattern.
    u_x, u_y = np.gradient( z );
    u_y_over_u_x = np.divide( u_y, u_x );

    if plot==True:
        fig = plt.figure()
        ax1 = plt.gca( projection = '3d' )
        surf = ax1.plot_surface(x,y,u_y_over_u_x, cmap = cm.coolwarm);
        ax1.set_xlabel('x'); ax1.set_ylabel('y');
        fig.suptitle(plotTitle)

        #fig = plt.figure()
        #ax2 = plt.gca( projection = '3d' )
        #surf = ax2.plot_surface(x2,y2,u_y_over_u_x2, cmap = cm.coolwarm);
        #ax2.set_xlabel('x'); ax1.set_ylabel('y');
        #fig.suptitle(plotTitle + ' edges removed')

    # Find average, throwing out outliers which are very large numbers
    # First Remove edge cases
    u_y_over_u_x2 = u_y_over_u_x[1:-1, 1:-1]
    u_y_over_u_x2 = u_y_over_u_x2[ np.isfinite(u_y_over_u_x2) ]
    print('The mean of u_y/u_x = ' + str(u_y_over_u_x2.mean()) + 
            ' with standard deviation = ' + str(u_y_over_u_x2.std()) )

    return u_y_over_u_x
