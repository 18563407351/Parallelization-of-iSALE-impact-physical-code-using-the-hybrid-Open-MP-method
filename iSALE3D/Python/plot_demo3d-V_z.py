import pySALEPlot as psp
import matplotlib.pyplot as plt
import numpy as np

# This example script plots a slice of the cell-based pressure
# field and a topography plot of the demo3d example

# Open model file
m = psp.opendatfile('jdata.dat', scale='km',plottype=['V_z'])

# Create output directory
outdir = 'VsTopo'
psp.mkdir_p(outdir)

# Loop over timesteps
for i in np.arange(m.nsteps):

    # Read pressure for timestep
    s = m.readStep('V_z', i)

    # calculate surface topography
    topo = m.surfaceTopography(step=s, topdown=True)

    # Create matplotlib figure and subplot
    fig, (ax1, ax2) = plt.subplots(nrows=2, figsize=(8, 8), sharex=True)
    ax1.set_aspect('equal')
    ax2.set_aspect('equal')

    # Plot slice along the y=0 boundary
    p = ax1.pcolormesh(m.x[m.ff], m.z[m.ff], s.V_z[m.ff],
                       cmap='viridis', vmin=-1e4, vmax=1e4)

    # Plot topography
    q = ax2.pcolormesh(m.x[m.z0], m.y[m.z0], topo,
                       vmin=-2, vmax=2, cmap='RdYlBu_r')

    # Set axes labels and limits
    ax1.set_xlabel('x [km]')
    ax1.set_ylabel('z [km]')
    ax1.set_xlim(m.xhires)  # The high resolution zone
    ax1.set_ylim(m.zhires)

    ax2.set_xlabel('x [km]')
    ax2.set_ylabel('y [km]')
    ax2.set_xlim(m.xhires)  # The high resolution zone
    ax2.set_ylim(m.yhires)

    # Create a colorbar
    cb1 = fig.colorbar(p, ax=ax1)
    cb1.set_label('Velocity')
    cb2 = fig.colorbar(q, ax=ax2)
    cb2.set_label('Topography [km]')
    #print(s.time)
    ax1.set_title('timestep: {:.5f}'.format(s.time))
    # Save the figure as a png
    fig.savefig('{}/{:03d}.png'.format(outdir, i), dpi=300)

    # Close figure ready for next timestep
    plt.close(fig)

