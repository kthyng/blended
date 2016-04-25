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
from BeautifulSoup import BeautifulSoup
import requests
from matplotlib.pylab import find
import matplotlib.tri as mtri
from scipy.spatial import Delaunay
import pdb
from pyproj import Proj
import octant


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

def rot2d(x, y, ang):
    '''rotate vectors by geometric angle'''
    xr = x*np.cos(ang) - y*np.sin(ang)
    yr = x*np.sin(ang) + y*np.cos(ang)
    return xr, yr


basemap = Proj(proj='utm', zone=15)#, ellps='clrk66',datum='NAD27')

# Read in shelf model grid
# USE NEW MODEL OUTPUT?
loc = 'http://barataria.tamu.edu:8080/thredds/dodsC/NcML/txla_nesting6.nc'
# shelfgrid = tracpy.inout.readgrid(loc, usebasemap=False, proj='utm', zone='15')
dshelf = netCDF.Dataset(loc)
tshelf = dshelf.variables['ocean_time']
datesshelf = netCDF.num2date(tshelf[:], tshelf.units)
# make u and v grid triangulations for interpolation
xrshelf, yrshelf = basemap(dshelf.variables['lon_rho'][:], dshelf.variables['lat_rho'][:])
xushelf, yushelf = basemap(dshelf.variables['lon_u'][:], dshelf.variables['lat_u'][:])
xvshelf, yvshelf = basemap(dshelf.variables['lon_v'][:], dshelf.variables['lat_v'][:])

mask_rhoshelf = dshelf.variables['mask_rho'][:]
iactrshelf = mask_rhoshelf.astype(bool)  # indices for active cells

# find u and v masks, code from octant.grid
mask_ushelf = mask_rhoshelf[:, 1:]*mask_rhoshelf[:, :-1]
mask_vshelf = mask_rhoshelf[1:, :]*mask_rhoshelf[:-1, :]
iactushelf = mask_ushelf.astype(bool)  # indices for active cells
iactvshelf = mask_vshelf.astype(bool)  # indices for active cells

# only send unmasked values to avoid have masked values interpolated into domain
pts = np.column_stack((xrshelf[iactrshelf].flatten(), yrshelf[iactrshelf].flatten()))
tess = Delaunay(pts)
trir = mtri.Triangulation(xrshelf[iactrshelf].flatten(), yrshelf[iactrshelf].flatten(), tess.simplices.copy())

pts = np.column_stack((xushelf[iactushelf].flatten(), yushelf[iactushelf].flatten()))
tess = Delaunay(pts)
triu = mtri.Triangulation(xushelf[iactushelf].flatten(), yushelf[iactushelf].flatten(), tess.simplices.copy())

pts = np.column_stack((xvshelf[iactvshelf].flatten(), yvshelf[iactvshelf].flatten()))
tess = Delaunay(pts)
triv = mtri.Triangulation(xvshelf[iactvshelf].flatten(), yvshelf[iactvshelf].flatten(), tess.simplices.copy())


# Read in blended grid
# grid_file = 'blended_grid.nc'
# # gridblend = tracpy.inout.readgrid(grid_file, usebasemap=False, proj='utm', zone='15')
# # this gives somewhat different differences between projected coordinates as compared with previous basemap
# # definition for the default values.
# gridblend = netCDF.Dataset(grid_file)
# xrblend, yrblend = basemap(gridblend.variables['lon_rho'][:], gridblend.variables['lat_rho'][:])
# mask_rhoblend = gridblend.variables['mask_rho'][:]
# # iactiver = gridblend.variables['mask_rho'][:].astype(bool)  # indices for active cells
# imaskrblend = ~mask_rhoblend.astype(bool)  # indices for masked cells
# xublend, yublend = basemap(gridblend.variables['lon_u'][:], gridblend.variables['lat_u'][:])
# xvblend, yvblend = basemap(gridblend.variables['lon_v'][:], gridblend.variables['lat_v'][:])
# xpsiblend, ypsiblend = basemap(gridblend.variables['lon_psi'][:], gridblend.variables['lat_psi'][:])

# # find u and v masks, code from octant.grid
# mask_ublend = mask_rhoblend[:, 1:]*mask_rhoblend[:, :-1]
# mask_vblend = mask_rhoblend[1:, :]*mask_rhoblend[:-1, :]
# imaskublend = ~mask_ublend.astype(bool)  # indices for masked cells
# imaskvblend = ~mask_vblend.astype(bool)  # indices for masked cells

grid_file = 'blended_grid.nc'
dgrid = netCDF.Dataset(grid_file)
lonvert = dgrid.variables['lon_vert'][:]
latvert = dgrid.variables['lat_vert'][:]
gridblend = octant.grid.CGrid_geo(lonvert, latvert, basemap)
maskblend = np.load('calcs/mask_blend.npz')['mask']
gridblend._set_mask_rho(maskblend)

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

# # Set up spatial weighting masks for the SUNTANS blended product
# # and the ROMS blended product.
# wrbayblend = np.zeros(mask_rhoblend.shape)
# wrshelfblend = np.ones(mask_rhoblend.shape)
# # only bay information should be used in the bay
# wrbayblend[162:, 246:551] = 1

# Read in the masks created in make_masks.py
maskbay = np.load('calcs/mask_bay.npz')['mask']
maskbayu = maskbay[:, 1:]*maskbay[:, :-1]  # code from octant _get_mask_u
maskbayv = maskbay[1:, :]*maskbay[:-1, :]

maskbayshelf = np.load('calcs/mask_bayshelf.npz')['mask']
maskbayshelfu = maskbayshelf[:, 1:]*maskbayshelf[:, :-1]  # code from octant _get_mask_u
maskbayshelfv = maskbayshelf[1:, :]*maskbayshelf[:-1, :]

maskshelf = np.load('calcs/mask_shelf.npz')['mask']
# other masks are exclusive at mask boundaries but need one to be
# inclusive so there are no missing sections of points coverage
maskshelfu = maskshelf[:, 1:]+maskshelf[:, :-1]  # code from octant _get_mask_u
maskshelfv = maskshelf[1:, :]+maskshelf[:-1, :]

# weights for bayshelf region
wsun = np.load('calcs/mask_bayshelf_weights_radial.npz')['wsun']
wsunu = wsun[:, 1:]  # don't want to interpolate this
wsunv = wsun[1:, :]

wroms = np.load('calcs/mask_bayshelf_weights_radial.npz')['wroms']
wromsu = wroms[:, 1:]  # don't want to interpolate this
wromsv = wroms[1:, :]

# # line for weighting across the bay shelf: y = m*x+b
# m = -1.229711996684854
# b = 3621774.9274811796

etablend = np.empty_like(gridblend.x_rho)

# Read in SUNTANS output
for File in Files:

    # klayer = 0 is surface, tstep=-99 reads in all time steps
    # Do u, v, eta
    sun = Spatial(File, klayer=2, tstep=-99, variable='eta')#, clim=clim)
    # sun = Spatial(File, klayer=0, tstep=-99, variable='eta')#, clim=clim)
    datesbay = sun.time  # get dates available in file
    eta = sun.loadData()

    # Loop through times in bay file
    for i, date in enumerate(datesbay):

        # Find the same time from the shelf model output
        # Assumes that there is a time in the shelf model that perfectly aligns with the
        # times in the bay model at some point
        tindshelf = find(datesshelf == date)
        print date
        # pdb.set_trace()
        if not find(datesshelf == date):  # if a shelf time doesn't align, don't use this bay model output
            continue

        # have to do u and v by time step for some reason
        u = Spatial(File, klayer=['surface'], tstep=i, variable='uc').loadData()
        v = Spatial(File, klayer=['surface'], tstep=i, variable='vc').loadData()

        # Interpolate model output from SUNTANS grid onto blended grid
        etabayblend = sun.interpolate(eta[i, :], gridblend.x_rho, gridblend.y_rho)
        ubayblend = sun.interpolate(u[:], gridblend.x_u, gridblend.y_u)
        vbayblend = sun.interpolate(v[:], gridblend.x_v, gridblend.y_v)

        # Rotate SUNTANS u,v velocities to be along/across
        # need to rotate negative to go onto curvilinear grid
        # make shifted array of vs in order to be able to rotate u
        vfake = np.empty_like(ubayblend)  # need v's to match u array size for rotating
        vfake[:vbayblend.shape[0], :] = tracpy.op.resize(vbayblend, 1)
        vfake[-1, :] = tracpy.op.resize(vbayblend[-1, :], 0)
        ubayblend, _ = rot2d(ubayblend, vfake, -tracpy.op.resize(gridblend.angle[:, 1:-1], 0))
        # make shifted array of us in order to be able to rotate v
        ufake = np.empty_like(vbayblend)  # need v's to match u array size for rotating
        ufake[:, :ubayblend.shape[1]] = tracpy.op.resize(ubayblend, 0)
        ufake[:, -1] = tracpy.op.resize(ubayblend[:, -1], 0)
        _, vbayblend = rot2d(ufake, vbayblend, -tracpy.op.resize(gridblend.angle[1:-1, :], 1))

        # Interpolate model output from TXLA onto blended grid
        etashelf = np.squeeze(dshelf.variables['zeta'][tindshelf, :, :])
        # function for interpolating from data on the shelf grids
        # only send unmasked values to avoid have masked values interpolated into domain
        fetashelf = mtri.LinearTriInterpolator(trir, etashelf[iactrshelf].flatten())
        # eta on the blended grid
        etashelfblend = fetashelf(gridblend.x_rho, gridblend.y_rho)
        # now mask out masked area of the blended grid
        etashelfblend = np.ma.masked_where(~gridblend.mask_rho.astype(bool), etashelfblend)

        ushelf = np.squeeze(dshelf.variables['u'][tindshelf, -1, :, :])  # at surface
        fushelf = mtri.LinearTriInterpolator(triu, ushelf[iactushelf].flatten())
        ushelfblend = fushelf(gridblend.x_u, gridblend.y_u)
        ushelfblend = np.ma.masked_where(~gridblend.mask_u.astype(bool), ushelfblend)

        vshelf = np.squeeze(dshelf.variables['v'][tindshelf, -1, :, :])  # at surface
        fvshelf = mtri.LinearTriInterpolator(triv, vshelf[iactvshelf].flatten())
        vshelfblend = fvshelf(gridblend.x_v, gridblend.y_v)
        vshelfblend = np.ma.masked_where(~gridblend.mask_v.astype(bool), vshelfblend)

        # Blended two model outputs together spatially
        import pdb; pdb.set_trace()
        etabayblend.fill_value = 0.0
        etashelfblend.fill_value = 0.0
        etablend = maskbay.astype(int)*etabayblend.filled() + \
            maskbayshelf.astype(int)*wsun*etabayblend.filled() + \
            maskbayshelf.astype(int)*wroms*etashelfblend.filled() + \
            maskshelf.astype(int)*etashelfblend.filled()
        etablend = np.ma.masked_where(~gridblend.mask_rho.astype(bool), etablend)

        # Examine all pieces going into eta blend calculation
        fig, axes = plt.subplots(2, 4, figsize=(18, 14), sharex=True, sharey=True)

        # column 0
        mappable = axes[0, 0].pcolormesh(gridblend.x_rho, gridblend.y_rho, maskbay.astype(int)*etabayblend.filled(), cmap=cmocean.cm.eta, vmin=-0.1, vmax=0.1)
        axes[0, 0].axis([270000, 370000, 3150000, 3320000])
        axes[0, 0].set_title('maskbay*etabayblend')
        fig.colorbar(mappable, ax=axes[0, 0])

        mappable = axes[1, 0].pcolormesh(gridblend.x_rho, gridblend.y_rho, etabayblend.filled(), cmap=cmocean.cm.eta, vmin=-0.1, vmax=0.1)
        axes[1, 0].contour(gridblend.x_rho, gridblend.y_rho, maskbay.astype(int), [0], colors='b')
        axes[1, 0].axis([270000, 370000, 3150000, 3320000])
        axes[1, 0].set_title('etabayblend')
        fig.colorbar(mappable, ax=axes[1, 0])

        # column 1
        mappable = axes[0, 1].pcolormesh(gridblend.x_rho, gridblend.y_rho, maskbayshelf.astype(int)*wsun*etabayblend.filled(), cmap=cmocean.cm.eta, vmin=-0.1, vmax=0.1)
        axes[0, 1].axis([270000, 370000, 3150000, 3320000])
        axes[0, 1].set_title('maskbayshelf*wsun*etabayblend')
        fig.colorbar(mappable, ax=axes[0, 1])

        mappable = axes[1, 1].pcolormesh(gridblend.x_rho, gridblend.y_rho, etabayblend.filled(), cmap=cmocean.cm.eta, vmin=-0.1, vmax=0.1)
        cs = axes[1, 1].contour(gridblend.x_rho, gridblend.y_rho, wsun, [0.0, 0.25, 0.5, 0.75, 1.0], colors='r')
        plt.clabel(cs, fmt='%1.2f')
        axes[1, 1].contour(gridblend.x_rho, gridblend.y_rho, maskbayshelf.astype(int), [0], colors='b')
        axes[1, 1].axis([270000, 370000, 3150000, 3320000])
        axes[1, 1].set_title('etabayblend')
        fig.colorbar(mappable, ax=axes[1, 1])

        # column 2
        mappable = axes[0, 2].pcolormesh(gridblend.x_rho, gridblend.y_rho, maskbayshelf.astype(int)*wroms*etashelfblend.filled(), cmap=cmocean.cm.eta, vmin=-0.1, vmax=0.1)
        axes[0, 2].axis([270000, 370000, 3150000, 3320000])
        axes[0, 2].set_title('maskbayshelf*wroms*etashelfblend')
        fig.colorbar(mappable, ax=axes[0, 2])

        mappable = axes[1, 2].pcolormesh(gridblend.x_rho, gridblend.y_rho, etashelfblend.filled(), cmap=cmocean.cm.eta, vmin=-0.1, vmax=0.1)
        cs = axes[1, 2].contour(gridblend.x_rho, gridblend.y_rho, wroms, [0, 0.25, 0.5, 0.75, 0.99], colors='r')
        plt.clabel(cs, fmt='%1.2f')
        axes[1, 2].contour(gridblend.x_rho, gridblend.y_rho, maskbayshelf.astype(int), [0], colors='b')
        axes[1, 2].axis([270000, 370000, 3150000, 3320000])
        axes[1, 2].set_title('etashelfblend')
        fig.colorbar(mappable, ax=axes[1, 2])

        # col 3
        mappable = axes[0, 3].pcolormesh(gridblend.x_rho, gridblend.y_rho, etashelfblend.filled(), cmap=cmocean.cm.eta, vmin=-0.1, vmax=0.1)
        axes[0, 3].contour(gridblend.x_rho, gridblend.y_rho, maskshelf.astype(int), [0], colors='b')
        axes[0, 3].axis([270000, 370000, 3150000, 3320000])
        axes[0, 3].set_title('etashelfblend')
        fig.colorbar(mappable, ax=axes[0, 3])
        # TOTAL
        mappable = axes[1, 3].pcolormesh(gridblend.x_rho, gridblend.y_rho, etablend.filled(), cmap=cmocean.cm.eta, vmin=-0.1, vmax=0.1)
        axes[1, 3].axis([270000, 370000, 3150000, 3320000])
        axes[1, 3].set_title('Total: etablend')
        fig.colorbar(mappable, ax=axes[1, 3])

        plt.tight_layout()
        fig.savefig('figures/comparisons/eta_pieces.png', bbox_inches='tight')

        etamax = max(etashelfblend.max(), etablend.max())
        etamin = min(etashelfblend.min(), etablend.min())
        fig = plt.figure(figsize=(16,8.375))
        ax1 = fig.add_subplot(1, 2, 1)
        ax1.pcolormesh(gridblend.x_rho, gridblend.y_rho, etashelfblend, vmin=etamin, vmax=etamax)
        ax1.axis([250000,400000,3150000,3350000])
        ax1.set_title('etashelfblend only')
        ax2 = fig.add_subplot(1, 2, 2)
        ax2.pcolormesh(gridblend.x_rho, gridblend.y_rho, etashelfblend, vmin=etamin, vmax=etamax)
        mappable = ax2.pcolormesh(gridblend.x_rho, gridblend.y_rho, etablend, vmin=etamin, vmax=etamax)
        cs = ax2.contour(gridblend.x_rho, gridblend.y_rho, wsun, [0, 0.25, 0.5, 0.75, 1], colors='k')
        plt.clabel(cs)
        cs = ax2.contour(gridblend.x_rho, gridblend.y_rho, wroms, [0, 0.25, 0.5, 0.75, 1], colors='b')
        plt.clabel(cs)
        ax2.axis([250000,400000,3150000,3350000])
        ax2.set_title('etashelfblend and etablend')
        fig.colorbar(mappable)
        fig.savefig('figures/comparisons/etacomp.png', bbox_inches='tight')



        ubayblend.fill_value = 0.0
        ushelfblend.fill_value = 0.0
        ublend = maskbayu.astype(int)*ubayblend.filled() + \
            maskbayshelfu.astype(int)*wsunu*ubayblend.filled() + \
            maskbayshelfu.astype(int)*wromsu*ushelfblend.filled() + \
            maskshelfu.astype(int)*ushelfblend.filled()
        ublend = np.ma.masked_where(~gridblend.mask_u.astype(bool), ublend)

        fig, axes = plt.subplots(2, 4, figsize=(18, 14), sharex=True, sharey=True)

        # column 0
        mappable = axes[0, 0].pcolormesh(gridblend.x_u, gridblend.y_u, maskbayu.astype(int)*ubayblend.filled(), cmap=cmocean.cm.u, vmin=-0.1, vmax=0.1)
        axes[0, 0].axis([270000, 370000, 3150000, 3320000])
        axes[0, 0].set_title('maskbayu*ubayblend')
        fig.colorbar(mappable, ax=axes[0, 0])

        mappable = axes[1, 0].pcolormesh(gridblend.x_u, gridblend.y_u, ubayblend.filled(), cmap=cmocean.cm.u, vmin=-0.1, vmax=0.1)
        axes[1, 0].contour(gridblend.x_u, gridblend.y_u, maskbayu.astype(int), [0], colors='b')
        axes[1, 0].axis([270000, 370000, 3150000, 3320000])
        axes[1, 0].set_title('ubayblend')
        fig.colorbar(mappable, ax=axes[1, 0])

        # column 1
        mappable = axes[0, 1].pcolormesh(gridblend.x_u, gridblend.y_u, maskbayshelfu.astype(int)*wsunu*ubayblend.filled(), cmap=cmocean.cm.u, vmin=-0.1, vmax=0.1)
        axes[0, 1].axis([270000, 370000, 3150000, 3320000])
        axes[0, 1].set_title('maskbayshelfu*wsunu*ubayblend')
        fig.colorbar(mappable, ax=axes[0, 1])

        mappable = axes[1, 1].pcolormesh(gridblend.x_u, gridblend.y_u, ubayblend.filled(), cmap=cmocean.cm.u, vmin=-0.1, vmax=0.1)
        cs = axes[1, 1].contour(gridblend.x_u, gridblend.y_u, wsunu, [0.0, 0.25, 0.5, 0.75, 1.0], colors='r')
        plt.clabel(cs, fmt='%1.2f')
        axes[1, 1].contour(gridblend.x_u, gridblend.y_u, maskbayshelfu.astype(int), [0], colors='b')
        axes[1, 1].axis([270000, 370000, 3150000, 3320000])
        axes[1, 1].set_title('ubayblend')
        fig.colorbar(mappable, ax=axes[1, 1])

        # column 2
        mappable = axes[0, 2].pcolormesh(gridblend.x_u, gridblend.y_u, maskbayshelfu.astype(int)*wromsu*ushelfblend.filled(), cmap=cmocean.cm.u, vmin=-0.1, vmax=0.1)
        axes[0, 2].axis([270000, 370000, 3150000, 3320000])
        axes[0, 2].set_title('maskbayshelfu*wromsu*ushelfblend')
        fig.colorbar(mappable, ax=axes[0, 2])

        mappable = axes[1, 2].pcolormesh(gridblend.x_u, gridblend.y_u, ushelfblend.filled(), cmap=cmocean.cm.u, vmin=-0.1, vmax=0.1)
        cs = axes[1, 2].contour(gridblend.x_u, gridblend.y_u, wromsu, [0, 0.25, 0.5, 0.75, 0.99], colors='r')
        plt.clabel(cs, fmt='%1.2f')
        axes[1, 2].contour(gridblend.x_u, gridblend.y_u, maskbayshelfu.astype(int), [0], colors='b')
        axes[1, 2].axis([270000, 370000, 3150000, 3320000])
        axes[1, 2].set_title('ushelfblend')
        fig.colorbar(mappable, ax=axes[1, 2])

        # col 3
        mappable = axes[0, 3].pcolormesh(gridblend.x_u, gridblend.y_u, ushelfblend.filled(), cmap=cmocean.cm.u, vmin=-0.1, vmax=0.1)
        axes[0, 3].contour(gridblend.x_u, gridblend.y_u, maskshelfu.astype(int), [0], colors='b')
        axes[0, 3].axis([270000, 370000, 3150000, 3320000])
        axes[0, 3].set_title('ushelfblend')
        fig.colorbar(mappable, ax=axes[0, 3])
        # TOTAL
        mappable = axes[1, 3].pcolormesh(gridblend.x_u, gridblend.y_u, ublend.filled(), cmap=cmocean.cm.u, vmin=-0.1, vmax=0.1)
        axes[1, 3].axis([270000, 370000, 3150000, 3320000])
        axes[1, 3].set_title('Total: ublend')
        fig.colorbar(mappable, ax=axes[1, 3])

        plt.tight_layout()
        fig.savefig('figures/comparisons/u_pieces.png', bbox_inches='tight')


        umax = max(ushelfblend.max(), ublend.max())*0.5
        umin = min(ushelfblend.min(), ublend.min())*0.5
        fig = plt.figure(figsize=(16,8.375))
        ax1 = fig.add_subplot(1, 2, 1)
        ax1.pcolormesh(gridblend.x_u, gridblend.y_u, ushelfblend, vmin=umin, vmax=umax)
        ax1.axis([250000,400000,3150000,3350000])
        ax1.set_title('ushelfblend only')
        ax2 = fig.add_subplot(1, 2, 2)
        ax2.pcolormesh(gridblend.x_u, gridblend.y_u, ushelfblend, vmin=umin, vmax=umax)
        mappable = ax2.pcolormesh(gridblend.x_u, gridblend.y_u, ublend, vmin=umin, vmax=umax)
        ax2.contour(gridblend.x_u, gridblend.y_u, wsunu, colors='k')
        ax2.axis([250000,400000,3150000,3350000])
        ax2.set_title('ushelfblend and ublend')
        fig.colorbar(mappable)
        fig.savefig('figures/comparisons/ucomp.png', bbox_inches='tight')


        vbayblend.fill_value = 0.0
        vshelfblend.fill_value = 0.0
        vblend = maskbayv.astype(int)*vbayblend.filled() + \
            maskbayshelfv.astype(int)*wsunv*vbayblend.filled() + \
            maskbayshelfv.astype(int)*wromsv*vshelfblend.filled() + \
            maskshelfv.astype(int)*vshelfblend.filled()
        vblend = np.ma.masked_where(~gridblend.mask_v.astype(bool), vblend)

        fig, axes = plt.subplots(2, 4, figsize=(18, 14), sharex=True, sharey=True)

        # column 0
        mappable = axes[0, 0].pcolormesh(gridblend.x_v, gridblend.y_v, maskbayv.astype(int)*vbayblend.filled(), cmap=cmocean.cm.v, vmin=-0.1, vmax=0.1)
        axes[0, 0].axis([270000, 370000, 3150000, 3320000])
        axes[0, 0].set_title('maskbayv*vbayblend')
        fig.colorbar(mappable, ax=axes[0, 0])

        mappable = axes[1, 0].pcolormesh(gridblend.x_v, gridblend.y_v, vbayblend.filled(), cmap=cmocean.cm.v, vmin=-0.1, vmax=0.1)
        axes[1, 0].contour(gridblend.x_v, gridblend.y_v, maskbayv.astype(int), [0], colors='b')
        axes[1, 0].axis([270000, 370000, 3150000, 3320000])
        axes[1, 0].set_title('vbayblend')
        fig.colorbar(mappable, ax=axes[1, 0])

        # column 1
        mappable = axes[0, 1].pcolormesh(gridblend.x_v, gridblend.y_v, maskbayshelfv.astype(int)*wsunv*vbayblend.filled(), cmap=cmocean.cm.v, vmin=-0.1, vmax=0.1)
        axes[0, 1].axis([270000, 370000, 3150000, 3320000])
        axes[0, 1].set_title('maskbayshelfv*wsunv*vbayblend')
        fig.colorbar(mappable, ax=axes[0, 1])

        mappable = axes[1, 1].pcolormesh(gridblend.x_v, gridblend.y_v, vbayblend.filled(), cmap=cmocean.cm.v, vmin=-0.1, vmax=0.1)
        cs = axes[1, 1].contour(gridblend.x_v, gridblend.y_v, wsunv, [0.0, 0.25, 0.5, 0.75, 1.0], colors='r')
        plt.clabel(cs, fmt='%1.2f')
        axes[1, 1].contour(gridblend.x_v, gridblend.y_v, maskbayshelfv.astype(int), [0], colors='b')
        axes[1, 1].axis([270000, 370000, 3150000, 3320000])
        axes[1, 1].set_title('vbayblend')
        fig.colorbar(mappable, ax=axes[1, 1])

        # column 2
        mappable = axes[0, 2].pcolormesh(gridblend.x_v, gridblend.y_v, maskbayshelfv.astype(int)*wromsv*vshelfblend.filled(), cmap=cmocean.cm.v, vmin=-0.1, vmax=0.1)
        axes[0, 2].axis([270000, 370000, 3150000, 3320000])
        axes[0, 2].set_title('maskbayshelfv*wromsv*vshelfblend')
        fig.colorbar(mappable, ax=axes[0, 2])

        mappable = axes[1, 2].pcolormesh(gridblend.x_v, gridblend.y_v, vshelfblend.filled(), cmap=cmocean.cm.v, vmin=-0.1, vmax=0.1)
        cs = axes[1, 2].contour(gridblend.x_v, gridblend.y_v, wromsv, [0, 0.25, 0.5, 0.75, 0.99], colors='r')
        plt.clabel(cs, fmt='%1.2f')
        axes[1, 2].contour(gridblend.x_v, gridblend.y_v, maskbayshelfv.astype(int), [0], colors='b')
        axes[1, 2].axis([270000, 370000, 3150000, 3320000])
        axes[1, 2].set_title('vshelfblend')
        fig.colorbar(mappable, ax=axes[1, 2])

        # col 3
        mappable = axes[0, 3].pcolormesh(gridblend.x_v, gridblend.y_v, vshelfblend.filled(), cmap=cmocean.cm.v, vmin=-0.1, vmax=0.1)
        axes[0, 3].contour(gridblend.x_v, gridblend.y_v, maskshelfv.astype(int), [0], colors='b')
        axes[0, 3].axis([270000, 370000, 3150000, 3320000])
        axes[0, 3].set_title('vshelfblend')
        fig.colorbar(mappable, ax=axes[0, 3])
        # TOTAL
        mappable = axes[1, 3].pcolormesh(gridblend.x_v, gridblend.y_v, vblend.filled(), cmap=cmocean.cm.v, vmin=-0.1, vmax=0.1)
        axes[1, 3].axis([270000, 370000, 3150000, 3320000])
        axes[1, 3].set_title('Total: vblend')
        fig.colorbar(mappable, ax=axes[1, 3])

        plt.tight_layout()
        fig.savefig('figures/comparisons/v_pieces.png', bbox_inches='tight')


        vmax = max(vshelfblend.max(), vblend.max())*0.5
        vmin = min(vshelfblend.min(), vblend.min())*0.5
        fig = plt.figure(figsize=(16,8.375))
        ax1 = fig.add_subplot(1, 2, 1)
        ax1.pcolormesh(gridblend.x_v, gridblend.y_v, vshelfblend, vmin=vmin, vmax=vmax)
        ax1.axis([250000,400000,3150000,3350000])
        ax1.set_title('vshelfblend only')
        ax2 = fig.add_subplot(1, 2, 2)
        ax2.pcolormesh(gridblend.x_v, gridblend.y_v, vshelfblend, vmin=vmin, vmax=vmax)
        mappable = ax2.pcolormesh(gridblend.x_v, gridblend.y_v, vblend, vmin=vmin, vmax=vmax)
        ax2.axis([250000,400000,3150000,3350000])
        ax2.set_title('vshelfblend and vblend')
        fig.colorbar(mappable)
        fig.savefig('figures/comparisons/vcomp.png', bbox_inches='tight')

# NOTE THAT THERE ARE SOME HOLES IN THE BAY BLENDED MODEL OUTPUT BECAUSE OF
# THE UNSTRUCTURED GRID OCCASIONALLY MISALIGNING WITH THE BLENDED GRID.
# THESE WILL PROBABLY NEED TO BE DEALT WITH.



# # Saved masking from SUNTANS model output, shows high res coast info
# np.savez('calcs/bayblendmask.npz', mask=~etabayblend.mask, x=xrblend, y=yrblend, lon=gridblend.variables['lon_rho'][:], lat=gridblend.variables['lat_rho'][:])

# ## Plot SUNTANS example: eta
# plt.figure()
# # Original SUNTANS output
# plt.scatter(sun.xv, sun.yv, c=eta[i,:], edgecolors='none')
# # Gridded, interpolated
# plt.scatter(xrblend.flatten(), yrblend.flatten(), c=etabayblend.flatten(), edgecolors='none', marker='x')

# ## Plot SUNTANS example: u
# plt.figure()
# plt.scatter(sun.xv, sun.yv, c=u, edgecolors='none')
# plt.scatter(xublend.flatten(), yublend.flatten(), c=ubayblend.flatten(), edgecolors='none', marker='x')

# ## Plot SUNTANS example: v
# plt.figure()
# plt.scatter(sun.xv, sun.yv, c=v, edgecolors='none')
# plt.scatter(xvblend.flatten(), yvblend.flatten(), c=vbayblend.flatten(), edgecolors='none', marker='x')

# ## Plot shelf example: eta
# plt.figure()
# plt.scatter(xrshelf.flatten(), yrshelf.flatten(), c=etashelf.flatten(), edgecolors='none', marker='o')
# plt.scatter(xrblend.flatten(), yrblend.flatten(), c=etashelfblend.flatten(), edgecolors='none', marker='x')

# ## Plot shelf example: u
# plt.figure()
# plt.scatter(xushelf.flatten(), yushelf.flatten(), c=ushelf.flatten(), edgecolors='none', marker='o')
# plt.scatter(xublend.flatten(), yublend.flatten(), c=ushelfblend.flatten(), edgecolors='none', marker='x')

# ## Plot shelf example: v
# plt.figure()
# plt.scatter(xvshelf.flatten(), yvshelf.flatten(), c=vshelf.flatten(), edgecolors='none', marker='o')
# plt.scatter(xvblend.flatten(), yvblend.flatten(), c=vshelfblend.flatten(), edgecolors='none', marker='x')

## Plot grids
# plt.figure()
# plt.plot(xpsiblend, ypsiblend, 'k', xpsiblend.T, ypsiblend.T, 'k')  # blended grid boxes (psi grid)
# plt.scatter(xrblend.flatten(), yrblend.flatten(), c= (~np.isnan(etabayblend)).astype(int).flatten(), edgecolors='none', marker='.', zorder=5)
# plt.plot(xrshelf, yrshelf, 'bo', ms=10)  # TXLA grid
# plt.plot(xrblend, yrblend, 'rx', ms=8)  # blended grid
# plt.plot(sun.xv, sun.yv, 'g.', ms=8)  # SUNTANS grid
