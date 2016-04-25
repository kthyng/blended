import pandas as pd
from pyproj import Proj
import netCDF4 as netCDF
from sunpy import Grid  # suntans code
from BeautifulSoup import BeautifulSoup
import requests
from matplotlib.pyplot import plt
import numpy as np
import xarray as xr
from datetime import datetime

basemap = Proj(proj='utm', zone=15)  # ellps='clrk66',datum='NAD27')
buoylon = -(94+53.943/60)
buoylat = 28+58.938/60
buoyx, buoyy = basemap(buoylon, buoylat)


def rot2d(x, y, ang):
    '''rotate vectors by geometric angle'''
    xr = x*np.cos(ang) - y*np.sin(ang)
    yr = x*np.sin(ang) + y*np.cos(ang)
    return xr, yr

## code from interpolate2grid.py
# Get list of urls for SUNTANS output available on thredds server
years = np.arange(2007, 2012)
Files = []  # urls for files on thredds server
for year in years:
    url = 'http://barataria.tamu.edu:8080/thredds/catalog/mrayson_galveston/' + str(year) + '/catalog.html'
    soup = BeautifulSoup(requests.get(url).text)

    # Pull out file name from webpage
    for row in soup.findAll('a'):
        # Only use netcdf files
        if not '.nc' in row.text:
            continue
        Filename = row.text.encode()  # File name from thredds catalog
        Files.append('http://barataria.tamu.edu:8080/thredds/dodsC/mrayson_galveston/' + str(year) + '/' + Filename)
Files = sorted(Files)  # now in order
##

dbay = netCDF.MFDataset(Files[:30])
# dbay = xr.open_mfdataset(Files[:30])  # January 2007, too slow
uc = dbay['uc']
vc = dbay['vc']
t = dbay['time']
grd = Grid(Files[0])

# Buoy B closest point in SUNTANS:
iBbay = 11964
# Make DataFrame for SUNTANS time series
# # use depth level "2" for now, which is 0.53 meters
# dfbay = pd.DataFrame(data={'uc': uc[:, 2, iBbay], 'vc': vc[:, 2, iBbay]}, index=netCDF.num2date(t[:], t.units))
# use depth level "3" for now, which is 0.79 meters
dfbay = pd.DataFrame(data={'uc': uc[:, 3, iBbay], 'vc': vc[:, 3, iBbay]}, index=netCDF.num2date(t[:], t.units))

# buoy data
dfB = pd.read_table('buoyB_Jan2007.txt', parse_dates=[[0, 1]], index_col=0, delim_whitespace=True, comment='#',
                    header=None, names=['Date', 'Time', 'East', 'North', 'Speed', 'Dir', 'WaterT'])

# shelf output
loc = 'http://barataria.tamu.edu:8080/thredds/dodsC/NcML/txla_nesting6.nc'
d = netCDF.Dataset(loc)
# x_rho, y_rho = proj(d['lon_rho'][:], d['lat_rho'][:])
# BuoyB close point in ROMS:
iBshelfx = 263
iBshelfy = 140
ds = xr.open_dataset(loc)
ucurv = ds['u'].sel(ocean_time='2007-01').isel(s_rho=-1, xi_u=iBshelfx, eta_u=iBshelfy)
vcurv = ds['v'].sel(ocean_time='2007-01').isel(s_rho=-1, xi_v=iBshelfx, eta_v=iBshelfy)
u, v = rot2d(ucurv, vcurv, ds['angle'][iBshelfy, iBshelfx])  # rotating to be cartesian
dfshelf = pd.DataFrame(data={'u': u, 'v': v}, index=ds['ocean_time'].sel(ocean_time='2007-01'))

# plot SUNTANS model
fig = plt.figure()
ax = fig.add_subplot(111)
grd.plotmesh(edgecolors=('grey',), facecolors=('None',), zorder=9)
Grid.calc_dg(grd)

# closest point to buoy B:
# grd.xv[iB], grd.yv[iB] (from grd.find_nearest([buoyx, buoyy]))
ax.plot(buoyx, buoyy, 'go')  # plot buoy B
ax.plot(grd.xv[iB], grd.yv[iB], 'rx')  # plot close SUNTANS point

# plot time series together
# Buoy data
plt.figure(figsize=(14, 6))
ax = dfB['East'].plot(color='k')
# SUNTANS output
ax.plot(dfbay.index, dfbay['uc']*100, 'r')
# Shelf output
ax.plot(dfshelf.index, dfshelf['u']*100, 'b')
plt.legend(['buoy', 'bay', 'shelf'], loc='best')
plt.title('Easting [cm/s]')
plt.tight_layout()
plt.savefig('figures/comparisons/tabs/u.pdf', bbox_inches='tight')

plt.figure(figsize=(14, 6))
ax = dfB['North'].plot(color='k')
# SUNTANS output
ax.plot(dfbay.index, dfbay['vc']*100, 'r')
# Shelf output
ax.plot(dfshelf.index, dfshelf['v']*100, 'b')
plt.legend(['buoy', 'bay', 'shelf'], loc='best')
plt.title('Northing [cm/s]')
plt.tight_layout()
plt.savefig('figures/comparisons/tabs/v.pdf', bbox_inches='tight')
