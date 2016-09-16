import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt
import pyroms

GRID = pyroms.grid.get_ROMS_grid('VIP')

bathy = GRID.vgrid.h
mask_rho = GRID.hgrid.mask_rho

bathy = np.ma.masked_where(mask_rho == 0, bathy)

fig = plt.figure(figsize=(20,10))

plt.pcolor(bathy)

plt.xlim(0,bathy.shape[1])
plt.ylim(0,bathy.shape[0])

plt.colorbar()

plt.tick_params(
    axis='x',          # changes apply to the x-axis
    which='both',      # both major and minor ticks are affected
    bottom='off',      # ticks along the bottom edge are off
    top='off',         # ticks along the top edge are off
    labelbottom='off') # labels along the bottom edge are off

plt.tick_params(
    axis='y',          # changes apply to the y-axis
    which='both',      # both major and minor ticks are affected
    left='off',        # ticks along the left edge are off
    right='off',       # ticks along the right edge are off
    labelleft='off')   # labels along the left edge are off

plt.show()
