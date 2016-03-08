import numpy as np
from matplotlib import animation
from fatiando import gridder
from fatiando.seismic import wavefd
from fatiando.vis import mpl
import matplotlib.pyplot as plt
# Set the parameters of the finite difference grid
shape = (601, 601)
area = [0, 60000, 0, 60000]
# Make a density and S wave velocity model
density = 2400 * np.ones(shape)
velocity = 3700
mu = wavefd.lame_mu(velocity, density)

# Make a wave source from a mexican hat wavelet
sources = [wavefd.GaussSource(35000, 35000, area, shape, 100, 1, delay=2)]

# Get the iterator for the simulation
dt = wavefd.maxdt(area, shape, velocity)
duration = 20
maxit = int(duration / dt)
stations = [[30000, 30000], [40000, 40000] ]# x, z coordinate of the seismometer
snapshot=None
simulation = wavefd.elastic_sh(mu, density, area, dt, maxit, sources, stations,
                               snapshot, padding=100, taper=0.01)

for t, u, seismograms in simulation:
    continue

times = np.linspace(0, duration, maxit)
u1=seismograms[0]
u2=seismograms[1]

plt.figure()
plt.plot(times, u1, 'x')
plt.plot(times, u2, 'o')
plt.show()
