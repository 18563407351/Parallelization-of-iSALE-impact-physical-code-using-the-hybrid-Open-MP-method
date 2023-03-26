import pySALEPlot as psp
import matplotlib.pyplot as plt
import numpy as np

# This example script plots a slice of the cell-based pressure
# field and a topography plot of the demo3d example

# Open model file
m = psp.opendatfile('jdata.dat', scale='m')
cg= m.craterGrowth(inc=5)
# Create output directory
outdir = 'CratDim'
psp.mkdir_p(outdir)



# Loop over timesteps
fig, (ax1, ax2) = plt.subplots(nrows=2, figsize=(8, 8), sharex=True)

p1 = ax1.scatter(cg[:, 1], cg[:, 2], color='r')
p2 = ax1.plot(cg[:, 1], cg[:, 2])

q1 = ax2.scatter(cg[:, 1],cg[:, 5], color='r')
q2 = ax2.plot(cg[:, 1],cg[:, 5])

ax1.set_xlabel('time [s]')
ax1.set_ylabel('Crater radius [km]')

ax2.set_xlabel('timestep [s]')
ax2.set_ylabel('Crater volume [m*3]')

fig.savefig('{}/CratDepVol.png'.format(outdir), dpi=300)

# Close figure ready for next timestep
plt.close(fig)
