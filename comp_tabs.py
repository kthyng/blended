import matplotlib.pyplot as plt
import numpy as np
from sunpy import Grid  # suntans code
import netCDF4 as netCDF
import read2df
import init
import os
import pandas as pd

# w = 0.0  # weight of SUNTANS output

locb = 'http://barataria.tamu.edu:8080/thredds/dodsC/2009/GalvCoarse_2009_AVG_0179.nc'

#
# # dates limitations here come from wind data, but also now have them from SUNTANS
# dstart = '2011-05-27'  # '2009-06-26', '2010-06-11', '2011-05-27'
# dend = '2011-08-02'  # '2009-11-13', '2010-08-20', '2011-08-02'
#
#
# plotting colors for the different data locations
# colors = ['k', 'purple', 'g']
colors = ['#76CA19', '#1976CA', '#CA1976']  # green, blue, pink; for data
# gshades = ['#76CA19', '#4A7E10', '#1D3306']  # green shades: data, bay, shelf
# bshades = ['#1976CA', '#104A7E', '#061D33']  # blue shades
# pshades = ['#CA1976', '#7E104A', '#33061D']  # pink shades
locnames = ['channel', 'med', 'far']
lws = {'data': 3, 'shelf': 2, 'bay': 1}
#
# # variable to plot
# var = 'al'


def plot_map(colors=colors):
    """Plot SUNTANS domain and data comp locations."""

    # SUNTANS bay
    grd = Grid(locb)

    # data locations (indices)
    ll, xy, ibays, ishelfs, basemap, shelfgrid = init.data_locs()

    x_rho, y_rho = basemap(shelfgrid['lon_rho'][:], shelfgrid['lat_rho'][:])

    fig = plt.figure(figsize=(6, 8))
    ax = fig.add_subplot(111)

    # plot SUNTANS model
    grd.plotmesh(edgecolors=('grey',), facecolors=('None',), zorder=9)
    # Grid.calc_dg(grd)

    # data locs
    for i in xrange(xy.shape[0]):

        ax.plot(xy[i, 0], xy[i, 1], 'o', color=colors[i], zorder=5)  # data
        ax.plot(grd.xv[ibays[i]], grd.yv[ibays[i]], 's', color=colors[i], zorder=5)  # plot close SUNTANS point
        ax.plot(x_rho[ishelfs[i, 1], ishelfs[i, 0]], y_rho[ishelfs[i, 1], ishelfs[i, 0]], '^', color=colors[i], zorder=5)  # plot close shelf point

    ax.set_title('Data locations', fontsize=14)

    fig.savefig('figures/comparisons/tabs/map.png', bbox_inches='tight')


def plot_spectra(dfs, var, colors=colors, lw=2, ls='-', legend=None, doss=False, dor2=False):
    """Input list of dataframes to plot.

    var is the keyword name of the column variable to plot.
    """

    lwsave = lw

    fig, ax = plt.subplots(2, 1, figsize=(8, 8))

    # loop over the types of data (bay, shelf, data)
    j = 0  # counter
    for i, df in enumerate(dfs):

        if isinstance(lwsave, list):
            lw = lwsave[i]
        else:
            lw = lwsave

        key = df.keys()[[var in key for key in df.keys()]][0]
        # i = np.where([str(num) in key for num in np.array([0,1,2])])[0][0]  # 0, 1, or 2

        # plot time series
        df[key].plot(ax=ax[0], color=colors[i], lw=lw, alpha=0.7, ls=ls)

        # calculate spectrum
        y = df[key]
        n = len(y)  # length of the signal
        if np.isnan(y).sum() > 0:
            igaps = np.where(np.isnan(y))[0][0]
            igape = np.where(np.isnan(y))[0][-1]
            y1 = y[:igaps]
            n1 = len(y1)
            y2 = y[igape+1:]
            n2 = len(y2)
            Y = [(np.fft.fft(y1)/n1)[range(n1/2)]]
            Y.append((np.fft.fft(y2)/n2)[range(n2/2)])
            k = [np.arange(n1)]
            k.append(np.arange(n2))
            Ts = 30*60  # sampling interval (sec)
            Fs = 1./Ts
            T = [n1/Fs]
            T.append(n2/Fs)
            frq = [(k[0]/T[0])[range(n1/2)]]  # two sides frequency range
            frq.append((k[0]/T[0])[range(n2/2)])
        else:
            Y = np.fft.fft(y)/n  # fft computing and normalization
            Y = Y[range(n/2)]
            k = np.arange(n)
            Ts = 30*60  # sampling interval (sec)
            Fs = 1./Ts
            T = n/Fs
            frq = k/T  # two sides frequency range
            frq = frq[range(n/2)]  # one side frequency range

        # plotting the spectrum
        if isinstance(Y, list):
            for frqs, Ys in zip(frq, Y):
                ax[1].loglog(frqs, abs(Ys), color=colors[i], lw=lw, alpha=0.7, ls=ls)
        else:
            ax[1].loglog(frq, abs(Y), color=colors[i], lw=lw, alpha=0.7, ls=ls)

        # legend
        if legend is not None:
            plt.legend(legend, loc='best', frameon=False)

        # skill score
        if doss:
            if len(dfs) == 1:  # need more than one entry to compare
                continue
            if i == 0:  # assume first entry is the data to be compared with
                continue
            data = dfs[0][key]  # just renaming
            model = df[key]  # renaming
            # labels and skill scores
            if np.isnan(data).sum() > 0:
                igaps = np.where(np.isnan(data))[0][0]
                igape = np.where(np.isnan(data))[0][-1]
                ax[0].text(0.05 + j*0.2, 0.03, legend[i] + ', ss = %1.2f, %1.2f' % (ss(data.iloc[:igaps], model.iloc[:igaps]), ss(data.iloc[igape+1:], model[igape+1:])), transform=ax[0].transAxes, color=colors[i])
            else:
                ax[0].text(0.05 + j*0.2, 0.03, legend[i] + ', ss = %1.2f' % ss(data, model), transform=ax[0].transAxes, color=colors[i])

        # r^2
        if dor2:
            # if len(dfs) == 1:  # need more than one entry to compare
            #     continue
            # if i == 0:  # assume first entry is the data to be compared with
            #     continue
            data = dfs[0][key]  # just renaming
            model = df[key]  # renaming
            # # labels and skill scores
            # if np.isnan(data).sum() > 0:
            #     igaps = np.where(np.isnan(data))[0][0]
            #     igape = np.where(np.isnan(data))[0][-1]
            #     ax[0].text(0.05 + i*0.2, 0.03, legend[i] + ', ss = %1.2f, %1.2f' % (ss(data.iloc[:igaps], model.iloc[:igaps]), ss(data.iloc[igape+1:], model[igape+1:])), transform=ax[0].transAxes, color=colors[i])
            # else:
            ax[0].text(0.5 + j*0.2, 0.03, legend[i] + ', r$^2$ = %1.2f' % r2(data, model), transform=ax[0].transAxes, color=colors[i])

        if doss or dor2:
            j += 1

    # make plot nice
    ax[0].set_title(key)
    ax[0].set_xlabel('Time')
    ax[0].set_ylabel('Amplitude')
    ax[0].set_ylim(-1.3, 1.3)
    ax[1].set_xlabel('Freq (Hz)')
    ax[1].set_ylabel('|Y(freq)|')
    ax[1].axis('tight')
    ax[1].set_ylim(1e-5, 0.5)
    ax[1].set_xlim(1e-7, 0.5e-4)
    fig.tight_layout()
    fig.savefig('figures/comparisons/tabs/spectrum.pdf', bbox_inches='tight')
    # fig.savefig('figures/comparisons/tabs/spectrum_' + var + key + dstart + '_' + dend + '.pdf', bbox_inches='tight')
    # plt.close(fig)


def ss(data, model):
    """Calculate skill score."""

    return 1 - np.sum((data - model)**2)/np.sum(data**2)


def r2(data, model):
    """Calculate correlation coefficient."""

    return (np.corrcoef(data, model)**2)[0,1]


    # ax.text(0.16, 0.18, 'buoy vs. bay: %1.2f' % subay, transform=ax.transAxes)
    # ax.text(0.16, 0.12, 'buoy vs. shelf: %1.2f' % sushelf, transform=ax.transAxes)
    # ax.text(0.16, 0.06, 'buoy vs. bay+shelf: %1.2f' % sushelfbay, transform=ax.transAxes)
    # # ax.text(0.16, 0.18, 'buoy vs. bay: %1.2f' % dfB['u'].corr(dfbay['u']), transform=ax.transAxes)
    # # ax.text(0.16, 0.12, 'buoy vs. shelf: %1.2f' % dfB['u'].corr(dfshelf['u']), transform=ax.transAxes)
    # # ax.text(0.16, 0.06, 'buoy vs. bay+shelf: %1.2f' % dfB['u'].corr(dfshelfbayu), transform=ax.transAxes)
    # plt.savefig('figures/comparisons/tabs/u_w' + str(w) + 'month' + month + '.pdf', bbox_inches='tight')
    #
    # # Sum of SUNTANS and shelf
    # # Combined output
    # dfshelfbayv = w*dfbay['v'] + dfshelf['v']
    # # dfshelfbayv = 0.0*dfbay['v'].resample('60T').interpolate() + dfshelf['v'].resample('60T').interpolate()
    # # skill for u, bay
    # svbay = 1-np.sum((dfB['v']-dfbay['v'])**2)/np.sum(dfB['v']**2)
    # # skill for u, shelf
    # svshelf = 1-np.sum((dfB['v']-dfshelf['v'])**2)/np.sum(dfB['v']**2)
    # # skill for u, bay+shelf
    # svshelfbay = 1-np.sum((dfB['v']-dfshelfbayv)**2)/np.sum(dfB['v']**2)
    # # # skill for u, bay
    # # svbay = 1-np.sum((dfB['v']-dfbay['v'])**2)/np.sum((dfB['v'] - dfbay['v']*0)**2)
    # # # skill for u, shelf
    # # svshelf = 1-np.sum((dfB['v']-dfshelf['v'])**2)/np.sum((dfB['v'] - dfshelf['v']*0)**2)
    # # # skill for u, bay+shelf
    # # svshelfbay = 1-np.sum((dfB['v']-dfshelfbayv)**2)/np.sum((dfB['v'] - dfshelfbayv*0)**2)
    # plt.figure(figsize=(14, 6))
    # ax = dfB['v'].plot(color='k', lw=4)
    # # SUNTANS output
    # ax.plot(dfbay.index, dfbay['v'], 'r', lw=2, alpha=0.7)
    # # Shelf output
    # ax.plot(dfshelf.index, dfshelf['v'], 'b', lw=2, alpha=0.7)
    # ax.plot(dfshelfbayv.index, dfshelfbayv, 'g.', lw=4)
    # plt.legend(['buoy', 'bay', 'shelf', 'bay+shelf'], loc='best', frameon=False)
    # plt.title('v [m/s]')
    # plt.tight_layout()
    # ax.text(0.16, 0.18, 'buoy vs. bay: %1.2f' % svbay, transform=ax.transAxes)
    # ax.text(0.16, 0.12, 'buoy vs. shelf: %1.2f' % svshelf, transform=ax.transAxes)
    # ax.text(0.16, 0.06, 'buoy vs. bay+shelf: %1.2f' % svshelfbay, transform=ax.transAxes)
    # # ax.text(0.16, 0.18, 'buoy vs. bay: %1.2f' % dfB['v'].corr(dfbay['v']), transform=ax.transAxes)
    # # ax.text(0.16, 0.12, 'buoy vs. shelf: %1.2f' % dfB['v'].corr(dfshelf['v']), transform=ax.transAxes)
    # # ax.text(0.16, 0.06, 'buoy vs. bay+shelf: %1.2f' % dfB['v'].corr(dfshelfbayv), transform=ax.transAxes)
    # plt.savefig('figures/comparisons/tabs/v_w' + str(w) + 'month' + month + '.pdf', bbox_inches='tight')

    # # testing
    # # # svshelf
    # # shelf1 = (dfB['v']-dfshelf['v'])
    # # shelf2 = (dfB['v']-dfshelfbayv)
    # 1-np.sum((dfB['v']-dfB['v'])**2)/np.sum(dfB['v']**2)

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

# Run code

# fname = 'calcs/df/bay_' + dstart + '_' + dend + '.pkl'
# if os.path.exists(fname):
#     dfbay = pd.read_pickle(fname)
# else:
#     dfbay = read2df.readbay(dstart, dend)
#     dfbay.to_pickle(fname)
#
# fname = 'calcs/df/shelf_' + dstart + '_' + dend + '.pkl'
# if os.path.exists(fname):
#     dfshelf = pd.read_pickle(fname)
# else:
#     dfshelf = read2df.readshelf(dstart, dend)
#     dfshelf.to_pickle(fname)
#
# fname = 'calcs/df/data_' + dstart + '_' + dend + '.pkl'
# if os.path.exists(fname):
#     dfdata = pd.read_pickle(fname)
# else:
#     dfdata = read2df.readdata(dstart, dend)
#     dfdata.to_pickle(fname)
#
# # Make a dictionary of dataframes
# df = {'bay': dfbay, 'shelf': dfshelf, 'data': dfdata}
#
# # make plots
# # plot_map()
#
# # plot data loc 0
# plot_spectra(df)
