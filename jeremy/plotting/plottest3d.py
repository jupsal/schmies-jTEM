from mpl_toolkits.mplot3d import axes3d
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib import cm

def generate(X, Y, phi):
    R = 1 - np.sqrt(X**2 + Y**2)
    return np.cos(2 * np.pi * X + phi) * R

fig = plt.figure()
ax = axes3d.Axes3D(fig)

xs = np.linspace(-1, 1, 50)
ys = np.linspace(-1, 1, 50)
X, Y = np.meshgrid(xs, ys)
Z = generate(X, Y, 0.0)
wframe = ax.plot_surface(X, Y, Z, rstride=2, cstride=2, cmap=cm.coolwarm )
ax.set_zlim(-1,1)

PHI = [ i * 360 / 2 / np.pi / 100 for i in xrange(10) ]

def update(i, ax, fig):
    ax.cla()
    #phi = i * 360 / 2 / np.pi / 100
    #Z = generate(X, Y, phi)
    Z = generate(X, Y, PHI[i])
    wframe = ax.plot_surface( X, Y, Z, rstride=2,
            cstride=2, cmap=cm.coolwarm )
    ax.set_zlim(-1,1)
    return wframe,

ani = animation.FuncAnimation( fig, update, frames=10, 
        fargs=(ax, fig), interval=100 )
plt.show()
