"""
Make weighting masks for blended model work.
"""

import numpy as np
from matplotlib.path import Path
import netCDF4 as netCDF
from pyproj import Proj
import matplotlib.pyplot as plt
import octant
from mesh import edit_mask_mesh

# basemap
basemap = Proj(proj='utm', zone=15)

# Load blended grid
grid_file = 'blended_grid.nc'
dgrid = netCDF.Dataset(grid_file)
# lon_rho = dgrid.variables['lon_rho'][:]
# lat_rho = dgrid.variables['lat_rho'][:]
# x_rho, y_rho = basemap(lon_rho, lat_rho)

# load grid via octant to get some useful properties
lonvert = dgrid.variables['lon_vert'][:]
latvert = dgrid.variables['lat_vert'][:]
grid = octant.grid.CGrid_geo(lonvert, latvert, basemap)
maskblend = np.load('calcs/mask_blend.npz')['mask']
grid._set_mask_rho(maskblend)

which = 'bay'

if which == 'bay':
    # For SUNTANS

    pass
    ## Bay ##
    # # Save approximate paths to start the process for the 3 masks
    # # pts = plt.ginput(n=0, timeout=0)
    # pts = np.asarray(pts)
    # ptsll = np.empty_like(pts)
    # ptsll[:, 0], ptsll[:, 1] = basemap(pts[:, 0], pts[:, 1], inverse=True)
    # np.savez('calcs/pts-path-bay.npz', pts=pts, ptsll=ptsll)

    # # form path
    # ptsll = np.load('calcs/pts-path-bay.npz')['ptsll']
    # path = Path(ptsll)  # Form path

    # # Use path to approximate weighting mask by finding points from blended grid
    # # that are within the path.
    # ptsgrid = np.vstack((grid.lon_rho.flatten(), grid.lat_rho.flatten()))
    # mask = path.contains_points(ptsgrid.T).reshape(grid.lon_rho.shape)

    # # Use GUI to edit bay mask
    # # Plot the high resolution bay grid info to improve the mask
    # edit_mask_mesh(grid.lon_vert, grid.lat_vert, mask)
    # pcolormesh(grid.lon_vert, grid.lat_vert, grid.mask_rho, edgecolors='k', alpha=0.5)
    # np.savez('calcs/mask_bay.npz', mask=mask)


elif which == 'bayshelf':

    # for SUNTANS near the coastline and increasing ROMS on the shelf

    pass

    ## Bay shelf ##
    # # pts = plt.ginput(n=0, timeout=0)
    # pts = np.asarray(pts)
    # ptsll = np.empty_like(pts)
    # ptsll[:, 0], ptsll[:, 1] = basemap(pts[:, 0], pts[:, 1], inverse=True)
    # np.savez('calcs/pts-path-bayshelf.npz', pts=pts, ptsll=ptsll)

    # # form path
    # ptsll = np.load('calcs/pts-path-bayshelf.npz')['ptsll']
    # path = Path(ptsll)  # Form path

    # # Use path to approximate weighting mask by finding points from blended grid
    # # that are within the path.
    # ptsgrid = np.vstack((grid.lon_rho.flatten(), grid.lat_rho.flatten()))
    # mask = path.contains_points(ptsgrid.T).reshape(grid.lon_rho.shape)

    # # Use GUI to edit bay mask
    # # Plot the high resolution bay grid info to improve the mask
    # edit_mask_mesh(grid.lon_vert, grid.lat_vert, mask)
    # # bay weighting mask
    # mask_bay = np.load('calcs/mask_bay.npz')['mask']
    # pcolormesh(grid.lon_vert, grid.lat_vert, mask_bay, edgecolors='k', alpha=0.5)
    # # suntans output mask, for area on shelf:
    # mask_bay_only = np.load('calcs/bayblendmask.npz')['mask']
    # pcolormesh(grid.lon_vert, grid.lat_vert, mask_bay_only, edgecolors='k', alpha=0.5)
    # # # blended mask:
    # # pcolormesh(grid.lon_vert, grid.lat_vert, grid.mask_rho, edgecolors='k', alpha=0.5)
    # np.savez('calcs/mask_bayshelf.npz', mask=mask)

    # # This leaves a jump in the resultant field. Try radial weighting instead
    # # get weighting for bay shelf region: 1 for suntans at coast and 1 for
    # # roms at edge, and linear in between
    # # first rotate the coordinates
    # theta = -0.6  # radians, found by guess and check and eyeball
    # x_vertrot = np.cos(theta)*grid.x_vert - np.sin(theta)*grid.y_vert
    # y_vertrot = np.sin(theta)*grid.x_vert + np.cos(theta)*grid.y_vert
    # wsun = np.zeros(maskbayshelf.shape)  # weights that go with the mask
    # wroms = np.zeros(maskbayshelf.shape)  # weights that go with the mask
    # for j in range(maskbayshelf.shape[1]):
    #     ind = np.where(maskbayshelf[:,j])[0]  # indices in this column that are masked
    #     # find the distance between each masked point and the first masked
    #     # point in the column
    #     if len(ind)>0:
    #         dist = np.sqrt((x_vertrot[ind,j] - x_vertrot[ind[-1],j])**2 +
    #                            (y_vertrot[ind,j] - y_vertrot[ind[-1],j])**2)
    #         weights = (dist.max() - dist)/dist.max()
    #         wsun[ind, j] = weights  # for suntans
    #         wroms[ind, j] = weights[::-1]  # for roms
    # np.savez('calcs/mask_bayshelf_weights.npz', mask=maskbayshelf, wsun=wsun, wroms=wroms)

    # distance along the coast for the suntans output. Have center of
    # weighting start at the middle of this along the coast.
    # dist = 112825  # meters
    # use grid.x_rho[204,392], grid.y_rho[204,392] at the center point
    xmid = grid.x_rho[204,392]
    ymid = grid.y_rho[204,392]
    r = np.sqrt((grid.x_rho-xmid)**2 + (grid.y_rho-ymid)**2)
    # this makes the zero circle that is within the suntans domain
    wsun = (47000 - r)/47000
    ind = wsun < 0
    wsun[ind] = 0.
    # wsun = (r.max() - r)/r.max()
    wroms = 1-wsun
    np.savez('calcs/mask_bayshelf_weights_radial.npz', mask=maskbayshelf,
             wsun=wsun, wroms=wroms, lon=grid.lon_rho, lat=grid.lat_rho)


elif which == 'restofshelf':

    # for ROMS output

    maskbay = np.load('calcs/mask_bay.npz')['mask']
    maskbayshelf = np.load('calcs/mask_bayshelf.npz')['mask']

    maskshelf = ~(maskbay + maskbayshelf)
    np.savez('calcs/mask_shelf.npz', mask=maskshelf)
