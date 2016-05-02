import pandas as pd
from pyproj import Proj
import netCDF4 as netCDF
from sunpy import Grid  # suntans code
from BeautifulSoup import BeautifulSoup
import requests
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr


month = '07'  # '01' or '03' or '07'
year = '2009'
w = 0.0  # weight of SUNTANS output

basemap = Proj(proj='utm', zone=15)  # ellps='clrk66',datum='NAD27')
buoylon, buoylat = -(94+53.943/60), 28+58.938/60
Elon, Elat = -(94+44.4/60), 29+20.5/60
Wlon, Wlat = -(94+47.760/60), 29+8.5516/60

def rot2d(x, y, ang):
    '''rotate vectors by geometric angle'''
    xr = x*np.cos(ang) - y*np.sin(ang)
    yr = x*np.sin(ang) + y*np.cos(ang)
    return xr, yr

## code from interpolate2grid.py
# Get list of urls for SUNTANS output available on thredds server
years = np.arange(2007, 2012)
Files = []  # urls for files on thredds server
for yeart in years:
    url = 'http://barataria.tamu.edu:8080/thredds/catalog/mrayson_galveston/' + str(yeart) + '/catalog.html'
    soup = BeautifulSoup(requests.get(url).text)

    # Pull out file name from webpage
    for row in soup.findAll('a'):
        # Only use netcdf files
        if not '.nc' in row.text:
            continue
        Filename = row.text.encode()  # File name from thredds catalog
        Files.append('http://barataria.tamu.edu:8080/thredds/dodsC/mrayson_galveston/' + str(yeart) + '/' + Filename)
Files = sorted(Files)  # now in order
##

# shelf output
loc = 'http://barataria.tamu.edu:8080/thredds/dodsC/NcML/txla_nesting6.nc'
d = netCDF.Dataset(loc)
x_rho, y_rho = basemap(d['lon_rho'][:], d['lat_rho'][:])

# Bay model output
inds = np.where([year + month in File for File in Files])[0]
inds = np.concatenate((inds, [inds[-1]+1]))
dbay = netCDF.MFDataset(np.asarray(Files)[inds])
# dbay = xr.open_mfdataset(Files[:30])  # January 2007, too slow
grd = Grid(Files[0])


# plot SUNTANS model
fig = plt.figure()
ax = fig.add_subplot(111)
grd.plotmesh(edgecolors=('grey',), facecolors=('None',), zorder=9)
Grid.calc_dg(grd)
# plot comparison locations
# closest point to buoy B: (middle)
# grd.xv[iB], grd.yv[iB] (from grd.find_nearest([buoyx, buoyy]))
buoyx, buoyy = basemap(buoylon, buoylat)
ax.plot(buoyx, buoyy, 'go')  # plot buoy B
iBbay = grd.find_nearest([buoyx, buoyy])[1]
ax.plot(grd.xv[iBbay], grd.yv[iBbay], 'gs')  # plot close SUNTANS point
iBshelfy, iBshelfx = 140, 263
ax.plot(x_rho[iBshelfy, iBshelfx], y_rho[iBshelfy, iBshelfx], 'g^')  # plot close shelf point
# medx, medy = 330000, 3226400  # location in projected coords
# imedbay = grd.find_nearest([medx, medy])[1]
# ax.plot(grd.xv[imedbay], grd.yv[imedbay], 'ks')
# imedshelfy, imedshelfx = 149, 272
# ax.plot(x_rho[imedshelfy, imedshelfx], y_rho[imedshelfy, imedshelfx], 'k^')  # plot close shelf point
# nearx, neary = 334500, 3248480
# inearbay = grd.find_nearest([nearx, neary])[1]
# ax.plot(grd.xv[inearbay], grd.yv[inearbay], 'bs')
# inearshelfy, inearshelfx = 162, 277
# ax.plot(x_rho[inearshelfy, inearshelfx], y_rho[inearshelfy, inearshelfx], 'b^')  # plot close shelf point
# coastx, coasty = 327723, 3237728
# icoastbay = grd.find_nearest([coastx, coasty])[1]
# ax.plot(grd.xv[icoastbay], grd.yv[icoastbay], 'rs')
# icoastshelfy, icoastshelfx = 157, 273
# ax.plot(x_rho[icoastshelfy, icoastshelfx], y_rho[icoastshelfy, icoastshelfx], 'r^')  # plot close shelf point
Ex, Ey = basemap(Elon, Elat)
ax.plot(Ex, Ey, 'o', color='darkorange')  # plot entrance channel NOAA buoy g06010
iEbay = grd.find_nearest([Ex, Ey])[1]
ax.plot(grd.xv[iEbay], grd.yv[iEbay], 's', color='darkorange')
iEshelfy, iEshelfx = 162, 276
ax.plot(x_rho[iEshelfy, iEshelfx], y_rho[iEshelfy, iEshelfx], '^', color='darkorange')  # plot close shelf point
Wx, Wy = basemap(Wlon, Wlat)
ax.plot(Wx, Wy, 'o', color='purple')  # plot entrance channel NOAA buoy g06010
iWbay = grd.find_nearest([Wx, Wy])[1]
ax.plot(grd.xv[iWbay], grd.yv[iWbay], 's', color='purple')
iWshelfy, iWshelfx = 162, 276
ax.plot(x_rho[iWshelfy, iWshelfx], y_rho[iWshelfy, iWshelfx], '^', color='purple')  # plot close shelf point
fig.savefig('figures/comparisons/tabs/map_data.png', bbox_inches='tight')



## MAKE DATAFRAMES ##
# shelf
ds = xr.open_dataset(loc)
alongB = ds['u'].sel(ocean_time=slice(year + '-' + month, year + '-' + str(int(month)+1) + '-01-00')).isel(s_rho=-1, xi_u=iBshelfx, eta_u=iBshelfy)
acrossB = ds['v'].sel(ocean_time=slice(year + '-' + month, year + '-' + str(int(month)+1) + '-01-00')).isel(s_rho=-1, xi_v=iBshelfx, eta_v=iBshelfy)
alongmed = ds['u'].sel(ocean_time=slice(year + '-' + month, year + '-' + str(int(month)+1) + '-01-00')).isel(s_rho=-1, xi_u=imedshelfx, eta_u=imedshelfy)
acrossmed = ds['v'].sel(ocean_time=slice(year + '-' + month, year + '-' + str(int(month)+1) + '-01-00')).isel(s_rho=-1, xi_v=imedshelfx, eta_v=imedshelfy)
alongnear = ds['u'].sel(ocean_time=slice(year + '-' + month, year + '-' + str(int(month)+1) + '-01-00')).isel(s_rho=-1, xi_u=inearshelfx, eta_u=inearshelfy)
acrossnear = ds['v'].sel(ocean_time=slice(year + '-' + month, year + '-' + str(int(month)+1) + '-01-00')).isel(s_rho=-1, xi_v=inearshelfx, eta_v=inearshelfy)
uB, vB = rot2d(alongB, acrossB, ds['angle'][iBshelfy, iBshelfx])  # rotating to be cartesian
umed, vmed = rot2d(alongmed, acrossmed, ds['angle'][imedshelfy, imedshelfx])  # rotating to be cartesian
unear, vnear = rot2d(alongnear, acrossnear, ds['angle'][inearshelfy, inearshelfx])  # rotating to be cartesian
dfshelf = pd.DataFrame(data={'uB': uB, 'vB': vB, 'alongB': alongB, 'acrossB': acrossB,
                       'umed': umed, 'vmed': vmed, 'alongmed': alongmed, 'acrossmed': acrossmed,
                       'unear': unear, 'vnear': vnear, 'alongnear': alongnear, 'acrossnear': acrossnear},
                       index=ds['ocean_time'].sel(ocean_time=slice(year + '-' + month, year + '-' + str(int(month)+1) + '-01-00')))
dfshelf['alongcoast'] = ds['u'].sel(ocean_time=slice(year + '-' + month, year + '-' + str(int(month)+1) + '-01-00')).isel(s_rho=-1, xi_u=icoastshelfx, eta_u=icoastshelfy)
dfshelf['acrosscoast'] = ds['v'].sel(ocean_time=slice(year + '-' + month, year + '-' + str(int(month)+1) + '-01-00')).isel(s_rho=-1, xi_v=icoastshelfx, eta_v=icoastshelfy)
dfshelf['ucoast'], dfshelf['vcoast'] = rot2d(dfshelf['alongcoast'], dfshelf['acrosscoast'], ds['angle'][icoastshelfy, icoastshelfx].data)  # rotating to be cartesian
dfshelf['alongE'] = ds['u'].sel(ocean_time=slice(year + '-' + month, year + '-' + str(int(month)+1) + '-01-00')).isel(s_rho=-1, xi_u=iEshelfx, eta_u=iEshelfy)
dfshelf['acrossE'] = ds['v'].sel(ocean_time=slice(year + '-' + month, year + '-' + str(int(month)+1) + '-01-00')).isel(s_rho=-1, xi_v=iEshelfx, eta_v=iEshelfy)
dfshelf['uE'], dfshelf['vE'] = rot2d(dfshelf['alongE'], dfshelf['acrossE'], ds['angle'][iEshelfy, iEshelfx].data)  # rotating to be cartesian
# resample to 30 min
dfshelf = dfshelf.resample('60T').interpolate()

# bay
uB, vB, t = dbay['uc'][:, 3, iBbay], dbay['vc'][:, 3, iBbay], dbay['time']
# Make DataFrame for SUNTANS time series
# # use depth level "2" for now, which is 0.53 meters
# dfbay = pd.DataFrame(data={'u': u[:, 2, iBbay], 'v': v[:, 2, iBbay]}, index=netCDF.num2date(t[:], t.units))
# use depth level "3" for now, which is 0.79 meters
dates = netCDF.num2date(t[:], t.units)
dfbay = pd.DataFrame(data={'uB': uB, 'vB': vB}, index=dates)
dfbay['alongB'], dfbay['acrossB'] = rot2d(dfbay['uB'], dfbay['vB'], -ds['angle'][iBshelfy, iBshelfx].data)  # rotating to be curvilinear
dfbay['umed'], dfbay['vmed'] = dbay['uc'][:, 3, imedbay], dbay['vc'][:, 3, imedbay]
dfbay['unear'], dfbay['vnear'] = dbay['uc'][:, 3, inearbay], dbay['vc'][:, 3, inearbay]
dfbay['ucoast'], dfbay['vcoast'] = dbay['uc'][:, 3, icoastbay], dbay['vc'][:, 3, icoastbay]
dfbay['uE'], dfbay['vE'] = dbay['uc'][:, 3, iEbay], dbay['vc'][:, 3, iEbay]
dfbay['alongmed'], dfbay['acrossmed'] = rot2d(dfbay['umed'], dfbay['vmed'], -ds['angle'][imedshelfy, imedshelfx].data)  # rotating to be curvilinear
dfbay['alongnear'], dfbay['acrossnear'] = rot2d(dfbay['unear'], dfbay['vnear'], -ds['angle'][inearshelfy, inearshelfx].data)  # rotating to be curvilinear
dfbay['alongcoast'], dfbay['acrosscoast'] = rot2d(dfbay['ucoast'], dfbay['vcoast'], -ds['angle'][icoastshelfy, icoastshelfx].data)  # rotating to be curvilinear
dfbay['alongE'], dfbay['acrossE'] = rot2d(dfbay['uE'], dfbay['vE'], -ds['angle'][iEshelfy, iEshelfx].data)  # rotating to be curvilinear
# resample to 30 min
dfbay = dfbay.resample('60T').interpolate()

# buoy data, buoy B
dfB = pd.read_table('buoyB_' + year + month + '.txt', parse_dates=[[0, 1]], index_col=0, delim_whitespace=True, comment='#',
                    header=None, names=['Date', 'Time', 'u', 'v', 'Speed', 'Dir', 'WaterT'])
# convert to m/s
dfB['u'] /= 100.
dfB['v'] /= 100.
dfB['along'], dfB['across'] = rot2d(dfB['u'], dfB['v'], -ds['angle'][iBshelfy, iBshelfx].data)  # rotating to be curvilinear
# resample to 30 min
dfB = dfB.resample('60T').interpolate()
# entrance channel, E
dfE = pd.read_table('g06010_' + year + month + '.txt', parse_dates=[[0, 1]], index_col=0, delim_whitespace=True, comment='#',
                    header=None, names=['Date', 'Time', 'Speed', 'Dir'])
# get as components and convert to m/s
dfE['u'] = (dfE['Speed']/100.)*np.cos((np.deg2rad(90-dfE['Dir'])))
dfE['v'] = (dfE['Speed']/100.)*np.sin((np.deg2rad(90-dfE['Dir'])))
dfE['along'], dfE['across'] = rot2d(dfE['u'], dfE['v'], -ds['angle'][iEshelfy, iEshelfx].data)  # rotating to be curvilinear
# resample to 30 min
dfE = dfE.resample('60T').interpolate()



## Spectrum plots ##

# comparison at buoy B location
y = dfbay['uB'].copy()
n = len(y)  # length of the signal
Y = np.fft.fft(y)/n  # fft computing and normalization
Ybay = Y[range(n/2)]
k = np.arange(n)
Ts = 30*60  # sampling interval (sec)
Fs = 1./Ts
T = n/Fs
frq = k/T  # two sides frequency range
frqbay = frq[range(n/2)]  # one side frequency range

y = dfshelf['uB'].copy()
n = len(y)  # length of the signal
Y = np.fft.fft(y)/n  # fft computing and normalization
Yshelf = Y[range(n/2)]
k = np.arange(n)
Ts = 30*60  # sampling interval (sec)
Fs = 1./Ts
T = n/Fs
frq = k/T  # two sides frequency range
frqshelf = frq[range(n/2)]  # one side frequency range

y = dfB['u'].copy()
n = len(y)  # length of the signal
Y = np.fft.fft(y)/n  # fft computing and normalization
Ybuoy = Y[range(n/2)]
k = np.arange(n)
Ts = 30*60  # sampling interval (sec)
Fs = 1./Ts
T = n/Fs
frq = k/T  # two sides frequency range
frqbuoy = frq[range(n/2)]  # one side frequency range

y = (dfbay['uB']+dfshelf['uB'])[year + '-' + month].copy()
n = len(y)  # length of the signal
Y = np.fft.fft(y)/n  # fft computing and normalization
Ycomb = Y[range(n/2)]
k = np.arange(n)
Ts = 30*60  # sampling interval (sec)
Fs = 1./Ts
T = n/Fs
frq = k/T  # two sides frequency range
frqcomb = frq[range(n/2)]  # one side frequency range

fig, ax = plt.subplots(2, 1, figsize=(12, 8))
ax[0].set_title('Buoy B location')
dfbay['uB'].plot(ax=ax[0], color='m', lw=2, alpha=0.7)
dfshelf['uB'].plot(ax=ax[0], color='b', lw=2, alpha=0.7)
dfB['u'].plot(ax=ax[0], color='g', lw=2, alpha=0.7)
(dfbay['uB']+dfshelf['uB']).plot(ax=ax[0], color='k', lw=2, alpha=0.7)
ax[0].set_xlabel('Time')
ax[0].set_ylabel('Amplitude')
ax[0].set_ylim(-1.0, 1.0)
ax[1].loglog(frqbay, abs(Ybay), 'm', lw=2, alpha=0.7)  # plotting the spectrum
ax[1].loglog(frqshelf, abs(Yshelf), 'b', lw=2, alpha=0.7)  # plotting the spectrum
ax[1].loglog(frqbuoy, abs(Ybuoy), 'g', lw=2, alpha=0.7)  # plotting the spectrum
ax[1].loglog(frqcomb, abs(Ycomb), 'k', lw=2, alpha=0.7)  # plotting the spectrum
ax[1].set_xlabel('Freq (Hz)')
ax[1].set_ylabel('|Y(freq)|')
ax[1].axis('tight')
ax[1].set_ylim(1e-5, 1e0)
fig.tight_layout()
fig.savefig('figures/comparisons/tabs/spectrum_buoyBloc.pdf', bbox_inches='tight')

# comparison at entrance channel buoy location
y = dfbay['uE'].copy()
n = len(y)  # length of the signal
Y = np.fft.fft(y)/n  # fft computing and normalization
Ybay = Y[range(n/2)]
k = np.arange(n)
Ts = 30*60  # sampling interval (sec)
Fs = 1./Ts
T = n/Fs
frq = k/T  # two sides frequency range
frqbay = frq[range(n/2)]  # one side frequency range

y = dfshelf['uE'].copy()
n = len(y)  # length of the signal
Y = np.fft.fft(y)/n  # fft computing and normalization
Yshelf = Y[range(n/2)]
k = np.arange(n)
Ts = 30*60  # sampling interval (sec)
Fs = 1./Ts
T = n/Fs
frq = k/T  # two sides frequency range
frqshelf = frq[range(n/2)]  # one side frequency range

y = dfE['u'].copy()
n = len(y)  # length of the signal
Y = np.fft.fft(y)/n  # fft computing and normalization
Ybuoy = Y[range(n/2)]
k = np.arange(n)
Ts = 30*60  # sampling interval (sec)
Fs = 1./Ts
T = n/Fs
frq = k/T  # two sides frequency range
frqbuoy = frq[range(n/2)]  # one side frequency range

y = (dfbay['uE']+dfshelf['uE'])[year + '-' + month].copy()
n = len(y)  # length of the signal
Y = np.fft.fft(y)/n  # fft computing and normalization
Ycomb = Y[range(n/2)]
k = np.arange(n)
Ts = 30*60  # sampling interval (sec)
Fs = 1./Ts
T = n/Fs
frq = k/T  # two sides frequency range
frqcomb = frq[range(n/2)]  # one side frequency range

fig, ax = plt.subplots(2, 1, figsize=(12, 8))
ax[0].set_title('Entrance channel buoy location')
dfbay['uE'].plot(ax=ax[0], color='m', lw=2, alpha=0.7)
dfshelf['uE'].plot(ax=ax[0], color='b', lw=2, alpha=0.7)
dfE['u'].plot(ax=ax[0], color='g', lw=2, alpha=0.7)
(dfbay['uE']+dfshelf['uE']).plot(ax=ax[0], color='k', lw=2, alpha=0.7)
ax[0].set_xlabel('Time')
ax[0].set_ylabel('Amplitude')
ax[0].set_ylim(-1.0, 1.0)
ax[1].loglog(frqbay, abs(Ybay), 'm', lw=2, alpha=0.7)  # plotting the spectrum
ax[1].loglog(frqshelf, abs(Yshelf), 'b', lw=2, alpha=0.7)  # plotting the spectrum
ax[1].loglog(frqbuoy, abs(Ybuoy), 'g', lw=2, alpha=0.7)  # plotting the spectrum
ax[1].loglog(frqcomb, abs(Ycomb), 'k', lw=2, alpha=0.7)  # plotting the spectrum
ax[1].set_xlabel('Freq (Hz)')
ax[1].set_ylabel('|Y(freq)|')
ax[1].axis('tight')
ax[1].set_ylim(1e-5, 1e0)
fig.tight_layout()
fig.savefig('figures/comparisons/tabs/spectrum_EntranceChannelbuoyloc.pdf', bbox_inches='tight')

# plot frequency of bay model output
y = dfbay['uB'].copy()
n = len(y)  # length of the signal
Y = np.fft.fft(y)/n  # fft computing and normalization
YB = Y[range(n/2)]
y = dfbay['umed'].copy()
n = len(y)  # length of the signal
Y = np.fft.fft(y)/n  # fft computing and normalization
Ymed = Y[range(n/2)]
y = dfbay['unear'].copy()
n = len(y)  # length of the signal
Y = np.fft.fft(y)/n  # fft computing and normalization
Ynear = Y[range(n/2)]
y = dfbay['ucoast'].copy()
n = len(y)  # length of the signal
Y = np.fft.fft(y)/n  # fft computing and normalization
Ycoast = Y[range(n/2)]
k = np.arange(n)
Ts = 30*60  # sampling interval (sec)
Fs = 1./Ts
T = n/Fs
frq = k/T  # two sides frequency range
frq = frq[range(n/2)]  # one side frequency range
periods = (1./frq)/(3600.*24)  # in days

fig, ax = plt.subplots(2, 1, figsize=(12, 8))
ax[0].set_title('Bay model u velocity')
dfbay['uB'].plot(ax=ax[0], color='g', lw=2, alpha=0.7)
dfbay['umed'].plot(ax=ax[0], color='k', lw=2, alpha=0.7)
# dfbay['unear'].plot(ax=ax[0], color='r')
dfbay['ucoast'].plot(ax=ax[0], color='r', lw=2, alpha=0.7)
ax[0].set_xlabel('Time')
ax[0].set_ylabel('Amplitude')
ax[0].set_ylim(-0.7, 0.7)
ax[1].loglog(frq, abs(YB), 'g', lw=2, alpha=0.7)  # plotting the spectrum
ax[1].loglog(frq, abs(Ymed), 'k', lw=2, alpha=0.7)  # plotting the spectrum
# ax[1].loglog(frq, abs(Ynear), 'r')  # plotting the spectrum
ax[1].loglog(frq, abs(Ycoast), 'r', lw=2, alpha=0.7)  # plotting the spectrum
ax[1].set_xlabel('Freq (Hz)')
ax[1].set_ylabel('|Y(freq)|')
ax[1].axis('tight')
ax[1].set_ylim(1e-5, 1e-1)
fig.tight_layout()
fig.savefig('figures/comparisons/tabs/spectrum_bay.pdf', bbox_inches='tight')

# plot frequency of shelf model output
y = dfshelf['uB'].copy()
n = len(y)  # length of the signal
Y = np.fft.fft(y)/n  # fft computing and normalization
YB = Y[range(n/2)]
y = dfshelf['umed'].copy()
n = len(y)  # length of the signal
Y = np.fft.fft(y)/n  # fft computing and normalization
Ymed = Y[range(n/2)]
y = dfshelf['unear'].copy()
n = len(y)  # length of the signal
Y = np.fft.fft(y)/n  # fft computing and normalization
Ynear = Y[range(n/2)]
y = dfshelf['ucoast'].copy()
n = len(y)  # length of the signal
Y = np.fft.fft(y)/n  # fft computing and normalization
Ycoast = Y[range(n/2)]
k = np.arange(n)
Ts = 30*60  # sampling interval (sec)
Fs = 1./Ts
T = n/Fs
frq = k/T  # two sides frequency range
frq = frq[range(n/2)]  # one side frequency range
periods = (1./frq)/(3600.*24)  # in days

fig, ax = plt.subplots(2, 1, figsize=(12, 8))
ax[0].set_title('Shelf model u velocity')
dfshelf['uB'].plot(ax=ax[0], color='g', lw=2, alpha=0.7)
dfshelf['umed'].plot(ax=ax[0], color='k', lw=2, alpha=0.7)
# dfshelf['unear'].plot(ax=ax[0], color='g', lw=2, alpha=0.7)
dfshelf['ucoast'].plot(ax=ax[0], color='r', lw=2, alpha=0.7)
ax[0].set_xlabel('Time')
ax[0].set_ylabel('Amplitude')
ax[0].set_ylim(-0.7, 0.7)
ax[1].loglog(frq, abs(YB), 'g', lw=2, alpha=0.7)  # plotting the spectrum
ax[1].loglog(frq, abs(Ymed), 'k', lw=2, alpha=0.7)  # plotting the spectrum
# ax[1].loglog(frq, abs(Ynear), 'g')  # plotting the spectrum
ax[1].loglog(frq, abs(Ycoast), 'r', lw=2, alpha=0.7)  # plotting the spectrum
ax[1].set_xlabel('Freq (Hz)')
ax[1].set_ylabel('|Y(freq)|')
ax[1].axis('tight')
ax[1].set_ylim(1e-5, 1e-1)
fig.tight_layout()
fig.savefig('figures/comparisons/tabs/spectrum_shelf.pdf', bbox_inches='tight')

# plot frequency of TABS data
y = dfB['u'].copy()
n = len(y)  # length of the signal
Y = np.fft.fft(y)/n  # fft computing and normalization
Y = Y[range(n/2)]
k = np.arange(n)
Ts = 30*60  # sampling interval (sec)
Fs = 1./Ts
T = n/Fs
frq = k/T  # two sides frequency range
frq = frq[range(n/2)]  # one side frequency range
periods = (1./frq)/(3600.*24)  # in days

fig, ax = plt.subplots(2, 1, figsize=(12, 8))
ax[0].set_title('Buoy u velocity')
dfB['u'].plot(ax=ax[0], color='g', lw=2, alpha=0.7)
ax[0].set_xlabel('Time')
ax[0].set_ylabel('Amplitude')
ax[0].set_ylim(-0.7, 0.7)
ax[1].loglog(frq, abs(Y), color='g', lw=2, alpha=0.7)  # plotting the spectrum
ax[1].set_xlabel('Freq (Hz)')
ax[1].set_ylabel('|Y(freq)|')
ax[1].axis('tight')
ax[1].set_ylim(1e-5, 1e-1)
fig.tight_layout()
fig.savefig('figures/comparisons/tabs/spectrum_buoyB.pdf', bbox_inches='tight')

# plot frequency of channel data
y = dfE['u'].copy()
n = len(y)  # length of the signal
Y = np.fft.fft(y)/n  # fft computing and normalization
Y = Y[range(n/2)]
k = np.arange(n)
Ts = 30*60  # sampling interval (sec)
Fs = 1./Ts
T = n/Fs
frq = k/T  # two sides frequency range
frq = frq[range(n/2)]  # one side frequency range
periods = (1./frq)/(3600.*24)  # in days

fig, ax = plt.subplots(2, 1, figsize=(12, 8))
ax[0].set_title('Buoy u velocity')
dfE['u'].plot(ax=ax[0], color='darkorange', lw=2, alpha=0.7)
ax[0].set_xlabel('Time')
ax[0].set_ylabel('Amplitude')
ax[0].set_ylim(-0.7, 0.7)
ax[1].loglog(frq, abs(Y), color='darkorange', lw=2, alpha=0.7)  # plotting the spectrum
ax[1].set_xlabel('Freq (Hz)')
ax[1].set_ylabel('|Y(freq)|')
ax[1].axis('tight')
ax[1].set_ylim(1e-5, 1e-1)
fig.tight_layout()
fig.savefig('figures/comparisons/tabs/spectrum_buoyE.pdf', bbox_inches='tight')


####

## CARTESIAN ##
# Sum of SUNTANS and shelf
# Combined output
# dfshelfbayu = 0.0*dfbay['u'].resample('60T').interpolate() + dfshelf['u'].resample('60T').interpolate()
dfshelfbayu = w*dfbay['u'] + dfshelf['u']
# skill for u, bay
subay = 1-np.sum((dfB['u']-dfbay['u'])**2)/np.sum(dfB['u']**2)
# skill for u, shelf
sushelf = 1-np.sum((dfB['u']-dfshelf['u'])**2)/np.sum(dfB['u']**2)
# skill for u, bay+shelf
sushelfbay = 1-np.sum((dfB['u']-dfshelfbayu)**2)/np.sum(dfB['u']**2)
# plot time series together
# Buoy data
plt.figure(figsize=(14, 6))
# ax = pd.rolling_mean(dfB['u'], window=3).plot(color='k', lw=4)
ax = dfB['u'].plot(color='k', lw=4)
# SUNTANS output
ax.plot(dfbay.index, dfbay['u'], 'r', lw=2, alpha=0.7)
# Shelf output
ax.plot(dfshelf.index, dfshelf['u'], 'b', lw=2, alpha=0.7)
ax.plot(dfshelfbayu.index, dfshelfbayu, 'g.', lw=4)
plt.legend(['buoy', 'bay', 'shelf', 'bay+shelf'], loc='best', frameon=False)
plt.title('u [m/s]')
plt.tight_layout()
ax.text(0.16, 0.18, 'buoy vs. bay: %1.2f' % subay, transform=ax.transAxes)
ax.text(0.16, 0.12, 'buoy vs. shelf: %1.2f' % sushelf, transform=ax.transAxes)
ax.text(0.16, 0.06, 'buoy vs. bay+shelf: %1.2f' % sushelfbay, transform=ax.transAxes)
# ax.text(0.16, 0.18, 'buoy vs. bay: %1.2f' % dfB['u'].corr(dfbay['u']), transform=ax.transAxes)
# ax.text(0.16, 0.12, 'buoy vs. shelf: %1.2f' % dfB['u'].corr(dfshelf['u']), transform=ax.transAxes)
# ax.text(0.16, 0.06, 'buoy vs. bay+shelf: %1.2f' % dfB['u'].corr(dfshelfbayu), transform=ax.transAxes)
plt.savefig('figures/comparisons/tabs/u_w' + str(w) + 'month' + month + '.pdf', bbox_inches='tight')

# Sum of SUNTANS and shelf
# Combined output
dfshelfbayv = w*dfbay['v'] + dfshelf['v']
# dfshelfbayv = 0.0*dfbay['v'].resample('60T').interpolate() + dfshelf['v'].resample('60T').interpolate()
# skill for u, bay
svbay = 1-np.sum((dfB['v']-dfbay['v'])**2)/np.sum(dfB['v']**2)
# skill for u, shelf
svshelf = 1-np.sum((dfB['v']-dfshelf['v'])**2)/np.sum(dfB['v']**2)
# skill for u, bay+shelf
svshelfbay = 1-np.sum((dfB['v']-dfshelfbayv)**2)/np.sum(dfB['v']**2)
# # skill for u, bay
# svbay = 1-np.sum((dfB['v']-dfbay['v'])**2)/np.sum((dfB['v'] - dfbay['v']*0)**2)
# # skill for u, shelf
# svshelf = 1-np.sum((dfB['v']-dfshelf['v'])**2)/np.sum((dfB['v'] - dfshelf['v']*0)**2)
# # skill for u, bay+shelf
# svshelfbay = 1-np.sum((dfB['v']-dfshelfbayv)**2)/np.sum((dfB['v'] - dfshelfbayv*0)**2)
plt.figure(figsize=(14, 6))
ax = dfB['v'].plot(color='k', lw=4)
# SUNTANS output
ax.plot(dfbay.index, dfbay['v'], 'r', lw=2, alpha=0.7)
# Shelf output
ax.plot(dfshelf.index, dfshelf['v'], 'b', lw=2, alpha=0.7)
ax.plot(dfshelfbayv.index, dfshelfbayv, 'g.', lw=4)
plt.legend(['buoy', 'bay', 'shelf', 'bay+shelf'], loc='best', frameon=False)
plt.title('v [m/s]')
plt.tight_layout()
ax.text(0.16, 0.18, 'buoy vs. bay: %1.2f' % svbay, transform=ax.transAxes)
ax.text(0.16, 0.12, 'buoy vs. shelf: %1.2f' % svshelf, transform=ax.transAxes)
ax.text(0.16, 0.06, 'buoy vs. bay+shelf: %1.2f' % svshelfbay, transform=ax.transAxes)
# ax.text(0.16, 0.18, 'buoy vs. bay: %1.2f' % dfB['v'].corr(dfbay['v']), transform=ax.transAxes)
# ax.text(0.16, 0.12, 'buoy vs. shelf: %1.2f' % dfB['v'].corr(dfshelf['v']), transform=ax.transAxes)
# ax.text(0.16, 0.06, 'buoy vs. bay+shelf: %1.2f' % dfB['v'].corr(dfshelfbayv), transform=ax.transAxes)
plt.savefig('figures/comparisons/tabs/v_w' + str(w) + 'month' + month + '.pdf', bbox_inches='tight')

# testing
# # svshelf
# shelf1 = (dfB['v']-dfshelf['v'])
# shelf2 = (dfB['v']-dfshelfbayv)
1-np.sum((dfB['v']-dfB['v'])**2)/np.sum(dfB['v']**2)

# ## CURVILINEAR ##
# # plot time series together
# # Buoy data
# plt.figure(figsize=(14, 6))
# ax = dfB['along'].plot(color='k', lw=4)
# # SUNTANS output
# ax.plot(dfbay.index, dfbay['along'], 'r', lw=2, alpha=0.7)
# # Shelf output
# ax.plot(dfshelf.index, dfshelf['along'], 'b', lw=2, alpha=0.7)
# # Sum of SUNTANS and shelf
# # Combined output
# dfshelfbayalong = dfbay['along'].resample('60T').interpolate() + dfshelf['along'].resample('60T').interpolate()
# ax.plot(dfshelfbayalong.index, dfshelfbayalong, 'g', lw=4)
# plt.legend(['buoy', 'bay', 'shelf', 'bay+shelf'], loc='best', frameon=False)
# plt.title('Along [m/s]')
# plt.tight_layout()
# ax.text(0.16, 0.18, 'buoy vs. bay: %1.2f' % dfB['along'].corr(dfbay['along']), transform=ax.transAxes)
# ax.text(0.16, 0.12, 'buoy vs. shelf: %1.2f' % dfB['along'].corr(dfshelf['along']), transform=ax.transAxes)
# ax.text(0.16, 0.06, 'buoy vs. bay+shelf: %1.2f' % dfB['along'].corr(dfshelfbayalong), transform=ax.transAxes)
# plt.savefig('figures/comparisons/tabs/along' + month + '.pdf', bbox_inches='tight')

# plt.figure(figsize=(14, 6))
# ax = dfB['across'].plot(color='k', lw=4)
# # SUNTANS output
# ax.plot(dfbay.index, dfbay['across'], 'r', lw=2, alpha=0.7)
# # Shelf output
# ax.plot(dfshelf.index, dfshelf['across'], 'b', lw=2, alpha=0.7)
# # Sum of SUNTANS and shelf
# # Combined output
# dfshelfbayacross = dfbay['across'].resample('60T').interpolate() + dfshelf['across'].resample('60T').interpolate()
# ax.plot(dfshelfbayacross.index, dfshelfbayacross, 'g', lw=4)
# plt.legend(['buoy', 'bay', 'shelf', 'bay+shelf'], loc='best', frameon=False)
# plt.title('Across [m/s]')
# plt.tight_layout()
# ax.text(0.16, 0.18, 'buoy vs. bay: %1.2f' % dfB['across'].corr(dfbay['across']), transform=ax.transAxes)
# ax.text(0.16, 0.12, 'buoy vs. shelf: %1.2f' % dfB['across'].corr(dfshelf['across']), transform=ax.transAxes)
# ax.text(0.16, 0.06, 'buoy vs. bay+shelf: %1.2f' % dfB['across'].corr(dfshelfbayacross), transform=ax.transAxes)
# plt.savefig('figures/comparisons/tabs/across' + month + '.pdf', bbox_inches='tight')

# # Plot variance
# uvar = np.var(u[:, 3, :], axis=0)
# vvar = np.var(v[:, 3, :], axis=0)
# # salt = dbay['salt']
# # saltvar = np.var(salt[:, 3, :], axis=0)

# fig = plt.figure(figsize=(8, 8))
# ax = fig.add_subplot(111)
# mappable = ax.scatter(grd.xv, grd.yv, c=vvar, s=100, cmap='viridis', vmax=0.1)
# ax.set_xlim(310000, 355000)
# ax.set_ylim(3230000, 3265000)
# fig.colorbar(mappable)
# # from make_masks.py
# cs = ax.contour(grid.x_rho, grid.y_rho, wsun, [0.01, 0.25, 0.5, 0.75, 0.99], colors='g')
# plt.clabel(cs)
# fig.savefig('figures/weights_vvar.png', bbox_inches='tight')
