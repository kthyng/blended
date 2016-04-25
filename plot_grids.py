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
from matplotlib.mlab import find
import matplotlib.tri as mtri
from pyproj import Proj
from scipy import ndimage
from mesh import edit_mask_mesh


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


basemap = Proj(proj='utm', zone=15)

# Read in shelf model grid
grid_file = '../../grid.nc'
# d = netCDF.Dataset(grid_file)
# using version of tracpy in branch update-grid !!!
grid = tracpy.inout.readgrid(grid_file, basemap)
# Now using normal ROMS convention instead of flipped arrays
# grid = octant.grid.CGrid_geo(lon_vert, lat_vert, basemap)

# loc = 'http://barataria.tamu.edu:8080/thredds/dodsC/NcML/txla_nesting6.nc'
# grid = tracpy.inout.readgrid(loc, usebasemap=True)

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
# fig = plt.figure()
# ax = fig.add_subplot(111)
# tracpy.plotting.background(grid)
# # Shelf grid
# dx = 1; dy = 1;
# ax.plot(grid['xr'][::dx], grid['yr'][::dy], \
#             grid['xr'][::dx].T, grid['yr'][::dy].T, color='darkcyan', alpha=0.4)
# # bay grid
# grd.plotmesh(edgecolors=('grey',), facecolors=('None',), zorder=9)
# Grid.calc_dg(grd)

lonB = [-95.29133464, -94.37435276]  # bounding lons for Bay grid, aligned with curvilinear shelf grid
latB = [28.9, 29.71269076]  # bounding lats for Bay grid, aligned with curvilinear shelf grid
# latB = [28.82617513, 29.71269076]  # bounding lats for Bay grid, aligned with curvilinear shelf grid

# # Calculate grid vertices from shelf rho grid so we have them for later
# d = netCDF.Dataset(loc)

# angle = d.variables['angle'][:]
# # Most things are in tracpy ordering, [x, y] not [y, x]
# xverts, yverts = octant.grid.rho_to_vert(grid['xr'], grid['yr'], grid['pm'], grid['pn'], angle.T)

# Use psi grid points as vertices for this blended product
xverts = grid.x_vert
yverts = grid.y_vert

## Interpolate vertices grid to be higher resolution in certain ranges ##
# indices on shelf grid that bound the Bay grid
xindsB, yindsB, _ = tracpy.tools.interpolate2d(lonB, latB, grid, 'd_ll2ij')
yindsB[0] = 147  # to be an integer
yindsB[1] = 189  # setting this manually because of unknown problem
# xindsB = [246.16400693039532 294.4152779769438] and I want it to be integers instead
xindsB[0] = 248
xindsB[1] = 291
# # projected locations 
# xpindsB, ypindsB = grid['basemap'](lonB, latB)
# mean resolution of shelf model in bay region
xres = (xverts[yindsB[0], xindsB[1]] - xverts[yindsB[0], xindsB[0]])/(xindsB[1] - xindsB[0])
yres = (yverts[yindsB[1], xindsB[1]] - yverts[yindsB[0], xindsB[1]])/(yindsB[1] - yindsB[0])
res = 300  # meters
# Set up grid space arrays for x and y, starting from shelf model but adding more entries for higher resolution
xvec = np.concatenate((np.arange(0, xindsB[0], 1), \
                        np.arange(xindsB[0], xindsB[1], 1./np.round((xres/res))), \
                        np.arange(xindsB[1], xverts.shape[1], 1)))
yvec = np.concatenate((np.arange(0, yindsB[0], 1), \
                        np.arange(yindsB[0], yverts.shape[0]-1, 1/np.round((yres/res)))))

# xgrid, ygrid are in [x, y] order
xgrid = xvec[np.newaxis, :].repeat(yvec.size, axis=0)
ygrid = yvec[:, np.newaxis].repeat(xvec.size, axis=1)
xverts_blend = ndimage.map_coordinates(grid.x_vert.T, np.array([xgrid.flatten(),
                                                      ygrid.flatten()]),
                             order=1, cval=0).reshape(xgrid.shape)
yverts_blend = ndimage.map_coordinates(grid.y_vert.T, np.array([xgrid.flatten(),
                                                      ygrid.flatten()]),
                             order=1, cval=0).reshape(xgrid.shape)

## Calculate lat/lon
lonverts_blend, latverts_blend = grid.proj(xverts_blend, yverts_blend, inverse=True)

# Calculate blended grid object
grid_blend = octant.grid.CGrid_geo(lonverts_blend, latverts_blend, grid.proj)

# # get rho grid in grid space from original vert grid, which we will then use
# # to interpolate the mask for the rho grid
# xgridrho = 0.25*(xgrid[1:, 1:] + xgrid[1:, :-1] +
#                    xgrid[:-1, 1:] + xgrid[:-1, :-1])
# ygridrho = 0.25*(ygrid[1:, 1:] + ygrid[1:, :-1] +
#                    ygrid[:-1, 1:] + ygrid[:-1, :-1])
# maskverts_blend = ndimage.map_coordinates(grid.mask_rho.T, np.array([xgridrho.flatten(),
#                                                       ygridrho.flatten()]),
#                              order=0, cval=0).reshape(xgridrho.shape)

# grid_blend._set_mask_rho(maskverts_blend)  # updates all of the masks

# fix mask to be correct with high res suntans input
# mask_bay_only = np.load('calcs/bayblendmask.npz')['mask']
# edit_mask_mesh(grid_blend.lon_vert, grid_blend.lat_vert, maskverts_blend)
# pcolormesh(grid_blend.lon_vert, grid_blend.lat_vert, mask_bay_only, edgecolors='k', alpha=0.5)
# np.savez('calcs/mask_blend.npz', mask=maskverts_blend)
maskblend = np.load('calcs/mask_blend.npz')['mask']
grid_blend._set_mask_rho(maskblend)

# Calculate blended mask
# set up a new triangulation on the verts instead of rho to be consistent
    # pts = np.column_stack((xr.flatten(), yr.flatten()))
    # tess = Delaunay(pts)
    # tri = mtri.Triangulation(xr.flatten(), yr.flatten(),
    #                          tess.simplices.copy())

# fmask = mtri.LinearTriInterpolator(grid.trir, grid.mask.flatten())
# mask_rho_blend = fmask(grid_blend.x_rho, grid_blend.y_rho)
# grid_blend._set_mask_rho(mask_rho_blend)  # updates all of the masks


### Deal with bathymetry later since it should be a combination of both grids. For now, just send in 1s.
### We might only need it for plotting since we are assuming surface-only drifter modeling.
## Interpolate bathymetry to new grid
# fh = grid['trir'].nn_interpolator(grid['h'].flatten(), mode='nearest')
# h_blend = fh(grid_blend.x_rho, grid_blend.y_rho)
h_blend = np.ones(grid_blend.x_rho.shape)


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
rootgrp.createDimension('xvert', xl+1)
rootgrp.createDimension('ypsi', yl-1)
rootgrp.createDimension('yr', yl)
rootgrp.createDimension('yvert', yl+1)
rootgrp.createDimension('zl', N)
rootgrp.createDimension('zlp1', N+1)
# Create variables
xverts = rootgrp.createVariable('x_vert','f8',('yvert','xvert'))
xpsis = rootgrp.createVariable('x_psi','f8',('ypsi','xpsi'))
xus = rootgrp.createVariable('x_u','f8',('yr','xpsi'))
xvs = rootgrp.createVariable('x_v','f8',('ypsi','xr'))
xrs = rootgrp.createVariable('x_rho','f8',('yr','xr'))
yverts = rootgrp.createVariable('y_vert','f8',('yvert','xvert'))
ypsis = rootgrp.createVariable('y_psi','f8',('ypsi','xpsi'))
yus = rootgrp.createVariable('y_u','f8',('yr','xpsi'))
yvs = rootgrp.createVariable('y_v','f8',('ypsi','xr'))
yrs = rootgrp.createVariable('y_rho','f8',('yr','xr'))
lonverts = rootgrp.createVariable('lon_vert','f8',('yvert','xvert'))
lonpsis = rootgrp.createVariable('lon_psi','f8',('ypsi','xpsi'))
lonus = rootgrp.createVariable('lon_u','f8',('yr','xpsi'))
lonvs = rootgrp.createVariable('lon_v','f8',('ypsi','xr'))
lonrs = rootgrp.createVariable('lon_rho','f8',('yr','xr'))
latverts = rootgrp.createVariable('lat_vert','f8',('yvert','xvert'))
latpsis = rootgrp.createVariable('lat_psi','f8',('ypsi','xpsi'))
latus = rootgrp.createVariable('lat_u','f8',('yr','xpsi'))
latvs = rootgrp.createVariable('lat_v','f8',('ypsi','xr'))
latrs = rootgrp.createVariable('lat_rho','f8',('yr','xr'))
mask_rhos = rootgrp.createVariable('mask_rho','f8',('yr','xr'))
mask_us = rootgrp.createVariable('mask_u','f8',('yr','xpsi'))
mask_vs = rootgrp.createVariable('mask_v','f8',('ypsi','xr'))
mask_psis = rootgrp.createVariable('mask_psi','f8',('ypsi','xpsi'))
pms = rootgrp.createVariable('pm','f8',('yr','xr'))
pns = rootgrp.createVariable('pn','f8',('yr','xr'))
angles = rootgrp.createVariable('angle','f8',('yvert','xvert'))
angle_rhos = rootgrp.createVariable('angle_rho','f8',('yr','xr'))
hs = rootgrp.createVariable('h','f8',('yr','xr'))
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
xverts[:] = grid_blend.x_vert
xpsis[:] = grid_blend.x_psi
xus[:] = grid_blend.x_u
xvs[:] = grid_blend.x_v
xrs[:] = grid_blend.x_rho
yverts[:] = grid_blend.y_vert
ypsis[:] = grid_blend.y_psi
yus[:] = grid_blend.y_u
yvs[:] = grid_blend.y_v
yrs[:] = grid_blend.y_rho
lonverts[:] = grid_blend.lon_vert
lonpsis[:] = grid_blend.lon_psi
lonus[:] = grid_blend.lon_u
lonvs[:] = grid_blend.lon_v
lonrs[:] = grid_blend.lon_rho
latverts[:] = grid_blend.lat_vert
latpsis[:] = grid_blend.lat_psi
latus[:] = grid_blend.lat_u
latvs[:] = grid_blend.lat_v
latrs[:] = grid_blend.lat_rho
mask_rhos[:] = grid_blend.mask_rho
mask_us[:] = grid_blend.mask_u
mask_vs[:] = grid_blend.mask_v
mask_psis[:] = grid_blend.mask_psi
pms[:] = grid_blend.pm
pns[:] = grid_blend.pn
angles[:] = grid_blend.angle
angle_rhos[:] = grid_blend.angle_rho
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
