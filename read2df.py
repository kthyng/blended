"""
Given time periods, read in model and data and return as
pandas dataframes.
"""

import numpy as np
from bs4 import BeautifulSoup
import requests
import pandas as pd
import netCDF4 as netCDF
import xarray as xr
import init
from datetime import datetime
from glob import glob


## Parameters ##
per = '30T'  # resampling period

# SUNTANS
# use depth level "3" for now, which is 0.79 meters
# 8 is 2.35672971 meters
# 12 is 4.12763299 meters
iz = 3

# ROMS
il = -1  # which layer to use
# loc = 'http://barataria.tamu.edu:8080/thredds/dodsC/NcML/txla_nesting6.nc'
loc = 'http://barataria.tamu.edu:8080/thredds/dodsC/NcML/txla_hindcast_agg'

# data locations (indices)
ll, xy, ibays, ishelfs, basemap, shelfgrid = init.data_locs()


def rot2d(x, y, ang):
    '''rotate vectors by geometric angle'''
    xr = x*np.cos(ang) - y*np.sin(ang)
    yr = x*np.sin(ang) + y*np.cos(ang)
    return xr, yr


def getbayfiles(years, which):
    """Get list of urls for SUNTANS output available on thredds server.
    code from interpolate2grid.py
    """

    if which == 'coarse':
        name = ''
    elif which == 'fine':
        name = 'FineTri_'
    elif which == 'quad':
        name = 'QuadWide_'

    # years = [2009]  # np.arange(2007, 2012)
    Files = []  # urls for files on thredds server
    for yeart in years:
        url = 'http://barataria.tamu.edu:8080/thredds/catalog/' + name + str(yeart) + '/catalog.html'
        soup = BeautifulSoup(requests.get(url).text)

        # Pull out file name from webpage
        for row in soup.findAll('a'):
            # Only use netcdf files
            if not '.nc' in row.text:
                continue
            Filename = row.text.encode()  # File name from thredds catalog
            if 'Harmonics' in Filename:
                continue  # not regular model output
            Files.append('http://barataria.tamu.edu:8080/thredds/dodsC/' + name + str(yeart) + '/' + Filename)
    return(sorted(Files))


def add2df(df, u, v, ind, i, inds=None):
    """Add u, v to existing df, along with curvilinear version.

    input inds is for rotating the bay model output using the shelf grid
    """

    if inds is not None:  # bay
        # Cartesian
        df['u' + str(i)], df['v' + str(i)] = u[:, ind], v[:, ind]
        # along/across rotation
        df['al' + str(i)], df['ac' + str(i)] = rot2d(df['u' + str(i)], df['v' + str(i)], -shelfgrid['angle'][inds[1], inds[0]])

    elif inds is None:  # shelf
        df['al' + str(i)] = u.isel(xi_u=ind[0], eta_u=ind[1])  # along-coast
        df['ac' + str(i)] = v.isel(xi_v=ind[0], eta_v=ind[1])  # across-coast
        # # Cartesian rotation (commented bc don't have angle in blended grid)
        # df['u' + str(i)], df['v' + str(i)] = rot2d(df['al' + str(i)], df['ac' + str(i)], shelfgrid['angle'][ind[1], ind[0]])  # rotating to be cartesian

    return df


def readblended(dstart, dend, Files):
    """Read in blended model output from loc for data locations.

    Files should be a globbed list of locations."""

    ds = xr.open_dataset(Files)

    inds = init.blended_datalocs()  # indices in blended grid for data locations

    # rename model output
    u = ds['u'].sel(ocean_time=slice(dstart, dend))
    v = ds['v'].sel(ocean_time=slice(dstart, dend))
    t = ds['ocean_time'].sel(ocean_time=slice(dstart, dend))

    u = u.rename({'yr': 'eta_u', 'xpsi': 'xi_u'})
    v = v.rename({'xr': 'xi_v', 'ypsi': 'eta_v'})

    # Initialize DataFrame
    df = pd.DataFrame(index=t)

    # add into dataframe the velocities at each comp location
    for i in range(inds.shape[0]):
        df = add2df(df, u, v, inds[i, :], i)

    # resample all of dataframe
    df = df.resample(per).interpolate()

    # limit date range to day after start until day before end
    # to get only full days
    # dstart = dstart[:8] + str(int(dstart[-2:]) + 1)
    # dend = dend[:8] + str(int(dend[-2:]) - 1)
    # df = df[dstart:dend]

    return df




def readbay(dstart, dend, Files=None):
    """Read in bay model output."""

    if Files is None:
        # Files = getbayfiles()  # this is when they are the thredds server
        # base = '/Volumes/rho.tamu.edu/'  # on tahoma
        base = '/rho/raid/'  # on hafen
        Files = glob(base + 'dongyu/2009/*.nc')
        Files.extend(glob(base + 'dongyu/2010/*.nc'))
        Files.extend(glob(base + 'dongyu/2011/*.nc'))

    dstartsplit = dstart.split('-')
    dstartdt = datetime(int(dstartsplit[0]), int(dstartsplit[1]), int(dstartsplit[2]))
    dendsplit = dend.split('-')
    denddt = datetime(int(dendsplit[0]), int(dendsplit[1]), int(dendsplit[2]))

    # Bay model output
    Filesuse = []
    for i, File in enumerate(Files):
        # print File
        ds = xr.open_dataset(File)
        # see if first time is in desired time range
        if (ds['time'].to_dataframe().index.to_datetime()[0] >= dstartdt
           and ds['time'].to_dataframe().index.to_datetime()[0] < denddt) or \
           (ds['time'].to_dataframe().index.to_datetime()[-1] >= dstartdt
           and ds['time'].to_dataframe().index.to_datetime()[-1] < denddt):
            Filesuse.append(File)
            continue
        ds.close()

    dbay = netCDF.MFDataset(np.asarray(Filesuse))

    # rename model output
    u, v, t = dbay['uc'][:, iz, :], dbay['vc'][:, iz, :], dbay['time']

    # Initialize DataFrame
    dates = netCDF.num2date(t[:], t.units)
    df = pd.DataFrame(index=dates)

    # add into dataframe the velocities at each comp location
    for i, ibay in enumerate(ibays.astype(int)):
        df = add2df(df, u, v, ibay, i, ishelfs[i, :])

    # resample all of dataframe
    df = df.resample(per).interpolate()

    # limit date range to day after start until day before end
    # to get only full days
    dstart = dstart[:8] + str(int(dstart[-2:]) + 1)
    dend = dend[:8] + str(int(dend[-2:]) - 1)
    df = df[dstart:dend]

    return df


def readshelf(dstart, dend):
    """Read in shelf model output."""

    ds = xr.open_dataset(loc)  # need this for all rotations

    # rename model output
    # ts = year + '-' + month  # starting time
    # te = year + '-' + str(int(month)+1) + '-01-00'  # ending time
    u = ds['u'].sel(ocean_time=slice(dstart, dend)).isel(s_rho=il)
    v = ds['v'].sel(ocean_time=slice(dstart, dend)).isel(s_rho=il)
    t = ds['ocean_time'].sel(ocean_time=slice(dstart, dend))

    # Initialize DataFrame
    df = pd.DataFrame(index=t)

    # add into dataframe the velocities at each comp location
    for i in range(ishelfs.shape[0]):
        df = add2df(df, u, v, ishelfs[i, :], i)

    # resample all of dataframe
    df = df.resample(per).interpolate()

    # limit date range to day after start until day before end
    # to get only full days
    dstart = dstart[:8] + str(int(dstart[-2:]) + 1)
    dend = dend[:8] + str(int(dend[-2:]) - 1)
    df = df[dstart:dend]

    return df


def readnoaa(dstart, dend):
    """Read in NOAA g06010 data."""

    # entrance channel
    df = pd.read_table('data/g06010_' + dstart + '_' + dend + '.txt', parse_dates=[[0, 1]], index_col=0, delim_whitespace=True, comment='#',
                       header=None, names=['Date', 'Time', 'Speed', 'Dir'])

    # get as components and convert to m/s
    df['u0'] = (df['Speed']/100.)*np.cos((np.deg2rad(90-df['Dir'])))
    df['v0'] = (df['Speed']/100.)*np.sin((np.deg2rad(90-df['Dir'])))
    df['al0'], df['ac0'] = rot2d(df['u0'], df['v0'], -shelfgrid['angle'][ishelfs[2, 1], ishelfs[2, 0]])  # rotating to be curvilinear

    # resample
    df = df.resample(per)#.interpolate()

    # limit date range to day after start until day before end
    # to get only full days
    dstart = dstart[:8] + str(int(dstart[-2:]) + 1)
    dend = dend[:8] + str(int(dend[-2:]) - 1)
    df = df[dstart:dend]

    # remove known large data gap
    if '2009' in dstart:
        df['2009-08-28 11:50':'2009-09-10 15:25:00'] = np.nan

    return df


def readother(dstart, dend):
    """Read in wind turbine study data."""

    # data from different times are combined in the file. Split up here.
    # date, time, speed, direction (2m depth)
    if '2009' in dstart:
        # 6/26/09-11/13/09
        df = pd.read_table('data/windfarm_rdcp.txt', parse_dates=[[0, 1]],
                           index_col=0, delim_whitespace=True,
                           header=None, names=['Date', 'Time', 'Speed', 'Dir'],
                           usecols=(0, 1, 36, 37), skiprows=0, nrows=2985)
    elif '2010' in dstart:
        # 6/9/10-8/20/10
        df = pd.read_table('data/windfarm_rdcp.txt', parse_dates=[[0, 1]],
                           index_col=0, delim_whitespace=True,
                           header=None, names=['Date', 'Time', 'Speed', 'Dir'],
                           usecols=(0, 1, 36, 37), skiprows=2985, nrows=4147-2986)
    elif '2011' in dstart:
        # 5/27/11-8/2/11
        df = pd.read_table('data/windfarm_rdcp.txt', parse_dates=[[0, 1]],
                           index_col=0, delim_whitespace=True,
                           header=None, names=['Date', 'Time', 'Speed', 'Dir'],
                           usecols=(0, 1, 36, 37), skiprows=4147)

    # get as components and convert to m/s
    df['u1'] = (df['Speed']/100.)*np.cos((np.deg2rad(90-df['Dir'])))
    df['v1'] = (df['Speed']/100.)*np.sin((np.deg2rad(90-df['Dir'])))
    df['al1'], df['ac1'] = rot2d(df['u1'], df['v1'], -shelfgrid['angle'][ishelfs[2, 1], ishelfs[2, 0]])  # rotating to be curvilinear

    # resample
    df = df.resample(per).interpolate()

    # limit date range to day after start until day before end
    # to get only full days
    dstart = dstart[:8] + str(int(dstart[-2:]) + 1)
    dend = dend[:8] + str(int(dend[-2:]) - 1)
    df = df[dstart:dend]

    return df


def readtabs(dstart, dend):
    """Read in TABS B buoy data"""

    # buoy data, buoy B
    df = pd.read_table('data/buoyB_' + dstart + '_' + dend + '.txt', parse_dates=[[0, 1]], index_col=0, delim_whitespace=True, comment='#',
                       header=None, names=['Date', 'Time', 'u2', 'v2', 'Speed', 'Dir', 'WaterT'])
    # convert to m/s
    df['u2'] /= 100.
    df['v2'] /= 100.
    df['al2'], df['ac2'] = rot2d(df['u2'], df['v2'], -shelfgrid['angle'][ishelfs[2, 1], ishelfs[2, 0]])  # rotating to be curvilinear

    # resample
    df = df.resample(per).interpolate()

    # limit date range to day after start until day before end
    # to get only full days
    dstart = dstart[:8] + str(int(dstart[-2:]) + 1)
    dend = dend[:8] + str(int(dend[-2:]) - 1)
    df = df[dstart:dend]

    return df


def readdata(dstart, dend):
    """Read in data into one dataframe."""

    df0 = readnoaa(dstart, dend)
    df1 = readother(dstart, dend)
    df2 = readtabs(dstart, dend)

    df = pd.concat([df0, df1, df2], axis=1)

    return df
