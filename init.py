"""
Initialization for running other scripts
"""

import numpy as np
from pyproj import Proj


def returnbasemap():
    basemap = Proj(proj='utm', zone=15)  # ellps='clrk66',datum='NAD27')
    return basemap


def data_locs():
    """Save data locations and model indices."""

    basemap = returnbasemap()

    # data locations, lon/lat, numbered with increasing distance from channel
    ll = np.zeros((3, 2))
    ll[0, :] = -(94+44.4/60), 29+20.5/60  # NOAA g06010 at Galveston Entrance channel
    ll[1, :] = -(94+47.760/60), 29+8.5516/60  # Steve wind turbine data
    ll[2, :] = -(94+53.943/60), 28+58.938/60  # TABS buoy B

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

    return ll, xy, ibays, ishelfs
