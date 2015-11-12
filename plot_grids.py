'''
Read in shelf model and Galveston Bay grid.
Look at spatial resolutions.
Decide how to blend grids for future drifter tracking.
'''

import numpy as np
import tracpy
import tracpy.plotting
import matplotlib as mpl
import matplotlib.pyplot as plt
import netCDF4 as netCDF
from scipy.spatial import ConvexHull
from mpl_toolkits.basemap import Basemap, pyproj
from sunpy import Grid # suntans code
import cmocean


mpl.rcParams.update({'font.size': 14})
mpl.rcParams['font.sans-serif'] = 'Arev Sans, Bitstream Vera Sans, Lucida Grande, Verdana, Geneva, Lucid, Helvetica, Avant Garde, sans-serif'
mpl.rcParams['mathtext.fontset'] = 'custom'
mpl.rcParams['mathtext.cal'] = 'cursive'
mpl.rcParams['mathtext.rm'] = 'sans'
mpl.rcParams['mathtext.tt'] = 'monospace'
mpl.rcParams['mathtext.it'] = 'sans:italic'
mpl.rcParams['mathtext.bf'] = 'sans:bold'
mpl.rcParams['mathtext.sf'] = 'sans'
mpl.rcParams['mathtext.fallback_to_cm'] = 'True'


# Read in shelf model grid
loc = 'http://barataria.tamu.edu:8080/thredds/dodsC/NcML/txla_nesting6.nc'
grid = tracpy.inout.readgrid(loc, usebasemap=True)

# Read in bay model grid
p = pyproj.Proj(proj='utm', zone='15')

# Galveston grid
bayloc = 'CoarseGrid/'
p_utm = np.loadtxt(bayloc + 'points_utm.dat')
lonp, latp = p(p_utm[:,0], p_utm[:,1], inverse=True)
xp, yp = grid['basemap'](lonp, latp)
p_lcc = np.zeros(p_utm.shape)
p_lcc[:,0] = xp
p_lcc[:,1] = yp
np.savetxt(bayloc + 'points.dat', p_lcc)
c_utm = np.loadtxt(bayloc + 'cells_utm.dat')
c_lcc = c_utm.copy()
lonp, latp = p(c_utm[:,0], c_utm[:,1], inverse=True)
xp, yp = grid['basemap'](lonp, latp)
c_lcc[:,0] = xp
c_lcc[:,1] = yp
np.savetxt(bayloc + 'cells.dat', c_lcc)
grd = Grid(bayloc)

# Plot the grids together - but need to zoom to see stuff
fig = plt.figure()
ax = fig.add_subplot(111)
tracpy.plotting.background(grid)
# Shelf grid
dx = 1; dy = 1;
ax.plot(grid['xr'][::dx], grid['yr'][::dy], \
            grid['xr'][::dx].T, grid['yr'][::dy].T, color='darkcyan', alpha=0.4)
# bay grid
grd.plotmesh(edgecolors=('grey',), facecolors=('None',), zorder=9)
Grid.calc_dg(grd)

# # Plot the grid resolution
# fig = plt.figure()
# ax = fig.add_subplot(111)
# tracpy.plotting.background(grid)
# # shelf
# mappable = ax.pcolormesh(grid['xr'], grid['yr'], np.sqrt(grid['pm']**-2 + grid['pn']**-2), cmap=cmocean.cm.gray)#, vmin=50, vmax=2000)
# plt.colorbar(mappable)
# # bay
# # ax.scatter(xp, yp, c=grd.dg, cmap=cmocean.cm.gray, s=22, linewidths=0.)

# NEED TO EXTEND SHELF GRID OUTSIDE OF DOMAIN! JUST GO TO EDGE FOR NOW.
lonB = [-95.29133464, -94.37435276]  # bounding lons for Bay grid, aligned with curvilinear shelf grid
latB = [28.82617513, 29.71269076]  # bounding lats for Bay grid, aligned with curvilinear shelf grid
# Interpolate u grid to be higher resolution in certain ranges
# indices on shelf u grid that bound the Bay grid
xuindsB, yuindsB, _ = tracpy.tools.interpolate2d(lonB, latB, grid, 'd_ll2ij')
# projected locations 
xpuindsB, ypuindsB = grid['basemap'](lonB, latB)
# mean resolution of shelf model in bay region
xres = (xpuindsB[1] - xpuindsB[0])/(xuindsB[1] - xuindsB[0])
yres = (ypuindsB[1] - ypuindsB[0])/(yuindsB[1] - yuindsB[0])
res = 100  # desired resolution in bay interpolation region
# Set up grid space arrays for x and y, starting from shelf model but adding more entries for higher resolution
xuvec = np.concatenate((np.arange(0, xuindsB[0], 1), \
                        np.arange(xuindsB[0], xuindsB[1], xres/res), \
                        np.arange(xuindsB[1], grid['xu'].shape[0], 1)))
yuvec = np.concatenate((np.arange(0, yuindsB[0], 1), \
                        np.arange(yuindsB[0], yuindsB[1], yres/res)))
                        # np.arange(yuindsB[1], grid['yu'].shape[0], 1)))

xugrid = xuvec[:, np.newaxis].repeat(yuvec.size, axis=1)
# xugrid = np.arange(0, grid['xu'].shape[0])[:, np.newaxis].repeat(grid['xu'].shape[1], axis=1)
yugrid = yuvec[np.newaxis, :].repeat(xuvec.size, axis=0)

xpu, ypu, _ = tracpy.tools.interpolate2d(xugrid, yugrid, grid, 'm_ij2xy')


# Plot shelf grid in grid space
fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(grid['xr'], grid['yr'], 'k', grid['xr'].T, grid['yr'].T, 'k')


