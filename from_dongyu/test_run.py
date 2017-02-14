
import tracpy
import tracpy.calcs
from tracpy.tracpy_class import Tracpy
import os
from datetime import datetime
import numpy as np
import pyproj
import netCDF4 as netCDF
import matplotlib.pyplot as plt
import cartopy
import cartopy.crs as ccrs




nsteps = 25

# Number of steps to divide model output for outputting drifter location
N = 5

# Number of days
ndays = 9

# This is a forward-moving simulation
ff = 1

# Time between outputs
tseas = 3600.  # time between output in seconds
ah = 0.
av = 0.  # m^2/s

# surface drifters
z0 = 's'

# for 3d flag, do3d=0 makes the run 2d and do3d=1 makes the run 3d
do3d = 0
doturb = 0

# for periodic boundary conditions in the x direction
doperiodic = 0

# Flag for streamlines. All the extra steps right after this are for streamlines.
dostream = 0


def init_test():
    time_units = 'seconds since 1970-01-01  00:00:00'
    zpar = 0  # 1 layer

    currents_filename = 'blended_uv_newer.nc'
    grid_filename = '../blended_grid.nc'

    proj = tracpy.tools.make_proj(setup='nwgom-pyproj')

    # Read in grid
    grid = tracpy.inout.readgrid(grid_filename, proj,
                                 vert_filename=currents_filename,
                                 usespherical=True)

    # Initialize Tracpy class
    tp = Tracpy(currents_filename, grid, name='test_blended', tseas=tseas,
                ndays=ndays, nsteps=nsteps, N=N, ff=ff, ah=ah, av=av,
                doturb=doturb, do3d=do3d, z0=z0, zpar=zpar,
                time_units=time_units,
                usespherical=True)
    lontemp = np.linspace(-96, -94, 30)
    lattemp = np.linspace(28, 30, 30)
    lon0, lat0 = np.meshgrid(lontemp, lattemp)

    return tp, lon0, lat0


def init_model():
    time_units = 'hours since 1993-01-01 01:00:00.000 UTC'
    zpar = 29

    currents_filename = 'http://barataria.tamu.edu:8080/thredds/dodsC/fmrc/txla_20yr_obc_his/out/TXLA_20yr_OBC_HIS_best.ncd'
    grid_filename = 'http://barataria.tamu.edu:8080/thredds/dodsC/fmrc/txla_20yr_obc_his/out/TXLA_20yr_OBC_HIS_best.ncd'

    proj = tracpy.tools.make_proj(setup='nwgom-pyproj')

    # Read in grid
    grid = tracpy.inout.readgrid(grid_filename, proj,
                                 vert_filename=currents_filename,
                                 usespherical=True)

    # Initialize Tracpy class
    tp = Tracpy(currents_filename, grid, name='test_model', tseas=tseas,
                ndays=ndays, nsteps=nsteps, N=N, ff=ff, ah=ah, av=av,
                doturb=doturb, do3d=do3d, z0=z0, zpar=zpar,
                time_units=time_units,
                usespherical=True)
    lontemp = np.linspace(-96, -94, 30)
    lattemp = np.linspace(28, 30, 30)
    lon0, lat0 = np.meshgrid(lontemp, lattemp)

    return tp, lon0, lat0


def run():

    # # run for test blended file
    # tp,lon0,lat0 = init_test()
    # date = datetime(2009,4,2, 0,0)
    # lonp, latp, zp, t, T0, U, V = tracpy.run.run(tp, date, lon0, lat0)

    # run for test original model file
    tp,lon0,lat0 = init_model()
    date = datetime(2009,4,2, 0,0)
    lonp, latp, zp, t, T0, U, V = tracpy.run.run(tp, date, lon0, lat0)


def plot():

    d = netCDF.Dataset('tracks/test_model.nc')
    fig = plt.figure()
    ax = plt.axes(projection=ccrs.Mercator(central_longitude=-85.0))
    gl = ax.gridlines(linewidth=0.2, color='gray', alpha=0.5, linestyle='-', draw_labels=True)
    ax.set_extent([-96.5, -92.5, 27, 30], ccrs.PlateCarree())
    ax.coastlines(resolution='10m')  # coastline resolution options are '110m', '50m', '10m'
    ax.plot(d['lonp'][:].T, d['latp'][:].T, 'k', alpha=0.5, transform=ccrs.PlateCarree());
    d = netCDF.Dataset('tracks/test_blended.nc')
    ax.plot(d['lonp'][:].T, d['latp'][:].T, 'r', alpha=0.5, transform=ccrs.PlateCarree());

if __name__ == "__main__":
    run()
