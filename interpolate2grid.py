'''
Interpolate SUNTANS model output onto blended grid.
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
from sunpy import Spatial


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
shelfgrid = tracpy.inout.readgrid(loc, usebasemap=True)

# Read in blended grid
grid = tracpy.inout.readgrid('blended_grid.nc', usebasemap=True)

# Read in SUNTANS output
sun = Spatial('http://barataria.tamu.edu:8080/thredds/catalog/mrayson_galveston/2007/catalog.html?dataset=mrayson-galveston/2007/GalvCoarse_20071231.nc')
eta = sun.loadData(variable='eta', tstep=0)

# Interpolate model output from SUNTANS grid onto blended grid
eta_interp = sun.interp(eta, my_x, my_y)


