"""
Initialization for running other scripts
"""

import numpy as np
from pyproj import Proj
import netCDF4 as netCDF


# data locations, lon/lat, numbered with increasing distance from channel
ll = np.zeros((3, 2))
ll[0, :] = -(94+44.4/60), 29+20.5/60  # NOAA g06010 at Galveston Entrance channel
ll[1, :] = -(94+47.760/60), 29+8.5516/60  # Steve wind turbine data
ll[2, :] = -(94+53.943/60), 28+58.938/60  # TABS buoy B

def returnbasemap():
    basemap = Proj(proj='utm', zone=15)  # ellps='clrk66',datum='NAD27')
    return basemap


def blended_datalocs():
    """Returns indices of data locations on blended grid."""

    # already ran: but WAS WRONG!
    # loc = 'http://copano.tamu.edu:8080/thredds/dodsC/NcML/txla_hindcast_agg'
    # import tracpy
    # proj = tracpy.tools.make_proj('nwgom-pyproj')
    # grid = tracpy.inout.readgrid(loc, proj)
    # ix, iy, _ = tracpy.tools.interpolate2d(ll[:,0], ll[:,1], grid, 'd_ll2ij')

    # j, i;  along across
    inds = np.zeros((7, 2), dtype=int)
    inds[0, :] = [446, 225]  # NOAA g06010 at Galveston Entrance channel
    inds[1, :] = [405, 161]  # wind
    inds[2, :] = [355, 140]  # TABS buoy B

    # new comparison points within bay
    inds[3, :] = [432, 244]  # (end of bay entrance channel)
    inds[4, :] = [313, 187]  # (small entrance channel)
    inds[5, :] = [435, 290]  # (mid bay)
    inds[6, :] = [470, 330]  # (upper bay)

    return inds


def data_locs():
    """Save data locations and model indices."""

    # ROMS grid only
    locshelfgrid = 'http://barataria.tamu.edu:8080/thredds/dodsC/txla_nesting6_grid/txla_grd_v4_new.nc'
    grid = netCDF.Dataset(locshelfgrid)

    basemap = returnbasemap()

    # blended grid
    # locblendedgrid = 'blended_grid.nc'
    # bgrid = netCDF.Dataset(locblendedgrid)

    # data locations, projected coords
    xy = np.zeros((3, 2))
    for i in range(xy.shape[0]):
        xy[i, :] = basemap(ll[i, 0], ll[i, 1])

    # SUNTANS grid indices
    ibays = np.zeros(3)
    ibays[0] = 9923
    ibays[1] = 14282
    ibays[2] = 11965
    # use the following in theory but too slow to always search
    # for i in range(ibays.size):
    #     ibays[i] = grd.find_nearest(xy[i, :])[1]

    # ROMS shelf grid indices (had to find these by hand)
    ishelfs = np.zeros((3, 2), dtype=int)
    ishelfs[0, :] = 276, 162
    ishelfs[1, :] = 270, 149
    ishelfs[2, :] = 263, 140

    return ll, xy, ibays, ishelfs, basemap, grid
