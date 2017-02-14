import cartopy
import matplotlib.pyplot as plt
import netCDF4 as netCDF
import cmocean.cm as cmo
import xarray as xr


# blended model output
d = netCDF.Dataset('blended_uv_newer.nc')
lon = d['lon'][:]
lat = d['lat'][:]
u = d['u'][:]
v = d['v'][:]
t = d['ocean_time']
dates = netCDF.num2date(t[:], t.units)

l10 = cartopy.feature.NaturalEarthFeature('physical', 'land', '10m')

# zoomed out
fig = plt.figure(figsize=(10,8))
ax = fig.add_subplot(111, projection=cartopy.crs.LambertConformal())
mappable = ax.pcolormesh(lon, lat, u[0,:,:], transform=cartopy.crs.PlateCarree(), cmap=cmo.delta)
ax.add_feature(l10, color='0.8')
fig.colorbar(mappable)
fig.savefig('figures/blended.png', bbox_inches='tight')

# zoomed in
fig = plt.figure(figsize=(10,8))
ax = fig.add_subplot(111, projection=cartopy.crs.LambertConformal())
mappable = ax.pcolormesh(lon, lat, u[0,:,:], transform=cartopy.crs.PlateCarree(), cmap=cmo.delta)
ax.add_feature(l10, color='0.8')
ax.set_extent([-95.4, -94.4, 28.9, 29.8], cartopy.crs.PlateCarree())
fig.colorbar(mappable)
fig.savefig('figures/blended-zoomed.png', bbox_inches='tight')
