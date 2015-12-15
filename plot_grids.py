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
import octant
from matplotlib import delaunay


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

# # Read in bay model grid
# p = pyproj.Proj(proj='utm', zone='15')
# # Galveston grid
# bayloc = 'CoarseGrid/'
# p_utm = np.loadtxt(bayloc + 'points_utm.dat')
# lonp, latp = p(p_utm[:,0], p_utm[:,1], inverse=True)
# xp, yp = grid['basemap'](lonp, latp)
# p_lcc = np.zeros(p_utm.shape)
# p_lcc[:,0] = xp
# p_lcc[:,1] = yp
# np.savetxt(bayloc + 'points.dat', p_lcc)
# c_utm = np.loadtxt(bayloc + 'cells_utm.dat')
# c_lcc = c_utm.copy()
# lonp, latp = p(c_utm[:,0], c_utm[:,1], inverse=True)
# xp, yp = grid['basemap'](lonp, latp)
# c_lcc[:,0] = xp
# c_lcc[:,1] = yp
# np.savetxt(bayloc + 'cells.dat', c_lcc)
# grd = Grid(bayloc)

# # Plot the grids together - but need to zoom to see stuff
fig = plt.figure()
ax = fig.add_subplot(111)
tracpy.plotting.background(grid)
# # Shelf grid
# dx = 1; dy = 1;
# ax.plot(grid['xr'][::dx], grid['yr'][::dy], \
#             grid['xr'][::dx].T, grid['yr'][::dy].T, color='darkcyan', alpha=0.4)
# # bay grid
# grd.plotmesh(edgecolors=('grey',), facecolors=('None',), zorder=9)
# Grid.calc_dg(grd)

lonB = [-95.29133464, -94.37435276]  # bounding lons for Bay grid, aligned with curvilinear shelf grid
latB = [28.82617513, 29.71269076]  # bounding lats for Bay grid, aligned with curvilinear shelf grid

# Calculate grid vertices from shelf rho grid so we have them for later
d = netCDF.Dataset(loc)
angle = d.variables['angle'][:]
# Most things are in tracpy ordering, [x, y] not [y, x]
xverts, yverts = octant.grid.rho_to_vert(grid['xr'], grid['yr'], grid['pm'], grid['pn'], angle.T)
# # cobble together verts mask, this is probably not right but may not matter much. CHECK THIS.
# mask_verts = np.empty((grid['mask'].shape[0]+1, grid['mask'].shape[1]+1))
# mask_verts[:-1, :-1] = grid['mask']
# mask_verts[-1, :-1] = grid['mask'][-1, :]
# mask_verts[:-1, -1] = grid['mask'][:, -1]
# mask_verts[-1, -1] = grid['mask'][-1, -1]

# # apply mask to verts
# xverts = np.ma.MaskedArray(xverts, mask_verts)
# yverts = np.ma.MaskedArray(yverts, mask_verts)

## Interpolate vertices grid to be higher resolution in certain ranges ##
# indices on shelf grid that bound the Bay grid
xindsB, yindsB, _ = tracpy.tools.interpolate2d(lonB, latB, grid, 'd_ll2ij')
yindsB[1] = 190  # setting this manually because of unknown problem
# projected locations 
xpindsB, ypindsB = grid['basemap'](lonB, latB)
# mean resolution of shelf model in bay region
xres = (xverts[xindsB[1], yindsB[0]] - xverts[xindsB[0], yindsB[0]])/(xindsB[1] - xindsB[0])
yres = (yverts[xindsB[1], yindsB[1]] - yverts[xindsB[1], yindsB[0]])/(yindsB[1] - yindsB[0])
# Set up grid space arrays for x and y, starting from shelf model but adding more entries for higher resolution
xvec = np.concatenate((np.arange(0, xindsB[0], 1), \
                        np.arange(xindsB[0], xindsB[1], 1/(xres/res)), \
                        np.arange(xindsB[1], xverts.shape[0], 1)))
yvec = np.concatenate((np.arange(0, yindsB[0], 1), \
                        np.arange(yindsB[0], yverts.shape[1], 1/(yres/res))))
# xgrid, ygrid are in [x, y] order
xgrid = xvec[:, np.newaxis].repeat(yvec.size, axis=1)
ygrid = yvec[np.newaxis, :].repeat(xvec.size, axis=0)
# xverts_blend, yverts_blend are the grid lines for the vertices grid for the blended model output
xverts_blend, yverts_blend, _ = tracpy.tools.interpolate2d(xgrid, ygrid, grid, 'm_ij2xy')

## Calculate lat/lon
lonverts_blend, latverts_blend = grid['basemap'](xverts_blend, yverts_blend, inverse=True)

# Calculate blended grid object
# Changing back to normal ROMS dimensions, [y,x]
grid_blend = octant.grid.CGrid_geo(lonverts_blend.T, latverts_blend.T, grid['basemap'])


## Interpolate bathymetry to new grid
fh = grid['trir'].nn_interpolator(grid['h'].flatten())
h_blend = fh(grid_blend.x_rho, grid_blend.y_rho)


## Save new grid ##
rootgrp = netCDF.Dataset('blended_grid.nc', 'w', format='NETCDF4')
# Set up basic grid information
xl = grid_blend.x_rho.shape[1]
yl = grid_blend.x_rho.shape[0]
N = 3
s_rho = np.linspace(-0.975, -0.025, N)
s_w = np.linspace(-1, 0, N+1)
hc = 0.
Cs_r = np.linspace(-0.975, -0.025, N)
Cs_w = np.linspace(-1, 0, N+1)
theta_b = 0.
hmin = 100.
theta_s = 1e-4
tcline = 0.
# Define dimensions
rootgrp.createDimension('xpsi', xl-1)
rootgrp.createDimension('xr', xl)
rootgrp.createDimension('ypsi', yl-1)
rootgrp.createDimension('yr', yl)
rootgrp.createDimension('zl', N)
rootgrp.createDimension('zlp1', N+1)
# Create variables
xpsis = rootgrp.createVariable('xpsi','f8',('xpsi','ypsi'))
xus = rootgrp.createVariable('xu','f8',('xpsi','yr'))
xvs = rootgrp.createVariable('xv','f8',('xr','ypsi'))
xrs = rootgrp.createVariable('xr','f8',('xr','yr'))
ypsis = rootgrp.createVariable('ypsi','f8',('xpsi','ypsi'))
yus = rootgrp.createVariable('yu','f8',('xpsi','yr'))
yvs = rootgrp.createVariable('yv','f8',('xr','ypsi'))
yrs = rootgrp.createVariable('yr','f8',('xr','yr'))
lonpsis = rootgrp.createVariable('lon_psi','f8',('xpsi','ypsi'))
lonus = rootgrp.createVariable('lon_u','f8',('xpsi','yr'))
lonvs = rootgrp.createVariable('lon_v','f8',('xr','ypsi'))
lonrs = rootgrp.createVariable('lon_rho','f8',('xr','yr'))
latpsis = rootgrp.createVariable('lat_psi','f8',('xpsi','ypsi'))
latus = rootgrp.createVariable('lat_u','f8',('xpsi','yr'))
latvs = rootgrp.createVariable('lat_v','f8',('xr','ypsi'))
latrs = rootgrp.createVariable('lat_rho','f8',('xr','yr'))
pms = rootgrp.createVariable('pm','f8',('xr','yr'))
pns = rootgrp.createVariable('pn','f8',('xr','yr'))
hs = rootgrp.createVariable('h','f8',('xr','yr'))
s_rhos = rootgrp.createVariable('s_rho','f8',('zl'))
s_ws = rootgrp.createVariable('s_w','f8',('zlp1'))
hcs = rootgrp.createVariable('hc','f8')
Cs_rs = rootgrp.createVariable('Cs_r','f8',('zl'))
Cs_ws = rootgrp.createVariable('Cs_w','f8',('zlp1'))
theta_bs = rootgrp.createVariable('theta_b','f8')
theta_ss = rootgrp.createVariable('theta_s','f8')
hmins = rootgrp.createVariable('hmin','f8')
tclines = rootgrp.createVariable('tcline','f8')
Vtransforms = rootgrp.createVariable('Vtransform','f8')
Vstretchings = rootgrp.createVariable('Vstretching','f8')

# Write data to netCDF variables
xpsis[:] = grid_blend.x_psi
xus[:] = grid_blend.x_u
xvs[:] = grid_blend.x_v
xrs[:] = grid_blend.x_rho
ypsis[:] = grid_blend.y_psi
yus[:] = grid_blend.y_u
yvs[:] = grid_blend.y_v
yrs[:] = grid_blend.y_rho
lonpsis[:] = grid_blend.lon_psi
lonus[:] = grid_blend.lon_u
lonvs[:] = grid_blend.lon_v
lonrs[:] = grid_blend.lon_rho
latpsis[:] = grid_blend.lat_psi
latus[:] = grid_blend.lat_u
latvs[:] = grid_blend.lat_v
latrs[:] = grid_blend.lat_rho
pms[:] = grid_blend.pm
pns[:] = grid_blend.pn
hs[:] = h_blend
s_rhos[:] = s_rho
s_ws[:] = s_w
hcs[:] = hc
Cs_rs[:] = Cs_r
Cs_ws[:] = Cs_w
theta_bs[:] = theta_b
theta_ss[:] = theta_s
hmins[:] = hmin
tclines[:] = tcline
Vtransforms[:] = 1
Vstretchings[:] = 1
rootgrp.close()
