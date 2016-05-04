###########################################################################
# This file holds the plotting routine for the KP solution from the Schmiesy
# Thesie.
############################################################################

import numpy as np
#import matplotlib
#matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib import cm
from mpl_toolkits.mplot3d import axes3d
from scipy import interpolate

# First load in the data
coordData = np.genfromtxt("./coords.csv", dtype=float, delimiter=',', \
        names=True); # comes in (t,x,y) tuples
solnData = np.genfromtxt("./soln.csv", dtype=float, delimiter=',', \
        names=True); # comes in (Real, Imaginary) pairs.

#print(coordData)
#print(solnData)

t = np.zeros(coordData.shape)
x = np.zeros(coordData.shape)
y = np.zeros(coordData.shape)
j = 0;
for txy in coordData:
    t[j] = txy[0]
    x[j] = txy[1]
    y[j] = txy[2]
    j = j+1;

print("Before check")

# Check that solution is entirely real-valued.
if sum( [ sol[1] for sol in solnData ] )>0: # then solution is not real-valued
    raise ValueError("The solution has nonzero imaginary part.")

print("After check")
realSoln = np.array( [ sol[0] for sol in solnData ] )
print("After")

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
zz = z.reshape(xx.shape)

fig = plt.figure()
ax = axes3d.Axes3D(fig)

# Plot for t=0.
surf = ax.plot_surface( xx, yy, zz, rstride=2, cstride=2, cmap=cm.coolwarm )
ax.set_zlim( np.min(realSoln), np.max(realSoln) )

def update(i, ax, fig):
    ax.cla()
    z = realSoln[ i*tstep:(i+1)*tstep ]
    zz = z.reshape( xx.shape )
    surf = ax.plot_surface( xx, yy, zz , rstride=2, cstride=2, cmap=cm.coolwarm )
    ax.set_zlim( np.min(realSoln), np.max(realSoln) )
    ax.set_xlabel('x'); ax.set_ylabel('y');
    return surf,

print len(realSoln)/tstep

ani = animation.FuncAnimation(fig, update, frames=len(realSoln)/tstep,
        fargs = (ax, fig), interval=1, save_count=200 )

# Save as gif. Requires imagemagick
#ani.save('animation2.gif', writer='imagemagick', fps=60)

#mywriter = animation.FFMpegWriter(fps=60)
#ani.save('/data/Documents/research/RiemannSurfaces/jTEM-Jeremy/jeremy/plotting/animation1.mp4',
#        writer=mywriter)

writer=animation.FFMpegFileWriter(fps=5)

ani.save('ANIMATION.mp4', fps=5, writer=writer)

#print("All saved up?")
#plt.show()
