import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
# from netCDF4 import Dataset
import os
import numpy as np
import numpy.ma as ma
import math
import pandas as pd
import time
from matplotlib import animation
import sys

sys.path.insert(0, '../CE')
from ModifyGeneratorCapacityWithWaterTData import getAllGridCellLatLongsInSpatFile, createAverageTFilename


def discrete_cmap(N, base_cmap=None):
    """Create an N-bin discrete colormap from the specified input map"""

    # Note that if base_cmap is a string or None, you can simply do
    #    return plt.cm.get_cmap(base_cmap, N)
    # The following works for string, None, or a colormap instance:

    base = plt.cm.get_cmap(base_cmap)
    color_list = base(np.linspace(0, 1, N))
    cmap_name = base.name + str(N)
    return base.from_list(cmap_name, color_list, N)


def draw_map_background(m, ax):
    ax.set_facecolor('#729FCF')
    m.fillcontinents(color='#FAFAFA', ax=ax, zorder=0, lake_color='#729FCF')
    m.drawcounties(ax=ax, color='#c5c5c5')
    m.drawrivers(ax=ax, color='#729FCF')
    m.drawstates(ax=ax)
    m.drawcountries(ax=ax)
    m.drawcoastlines(ax=ax)
    # draw parallels and meridians.
    # label parallels on right and top
    # meridians on bottom and left
    parallels = np.arange(m.boundarylats[0], m.boundarylats[1], 10.)
    # labels = [left,right,top,bottom]
    m.drawparallels(parallels, labels=[False, True, False, False], size='xx-small', color='gray')
    meridians = np.arange(m.boundarylonmin, m.boundarylonmax, 10.)
    m.drawmeridians(meridians, labels=[False, False, False, True], size='xx-small', color='gray')


def population_map():
    dataset = Dataset(os.path.expanduser('~/Downloads/gpw-v4-population-count-rev10_totpop_2pt5_min_nc/gpw_v4_e_atotpopbt_cntm_2pt5_min.nc'))

    # Extract data from NetCDF file
    lats = dataset.variables['latitude'][:]
    lons = dataset.variables['longitude'][:]
    time = dataset.variables['raster'][:]
    values = dataset.variables['Population Count, v4.10 (2000, 2005, 2010, 2015, 2020): 2.5 arc-minutes'][:]

    years = [2000, 2005, 2010, 2015, 2020]

    ll_lat = 24
    ll_lon = -125
    ur_lat = 50
    ur_lon = -60

    #ll_lat = 29.65
    #ll_lon = -94.143
    #ur_lat = 38.77
    #ur_lon = -76.56

    for i, y in enumerate(years):
        print('Creating plot for year {0:4d}'.format(y))
        pop = ma.array(values[i, :, :])

        # get mask
        mm = pop.mask

        masklon = (lons > ll_lon) & (lons < ur_lon)
        masklat = (lats > ll_lat) & (lats < ur_lat)

        m = Basemap(llcrnrlon=ll_lon, llcrnrlat=ll_lat, urcrnrlon=ur_lon, urcrnrlat=ur_lat,
                    resolution='i', area_thresh=2500.)
        fig = plt.figure()
        ax = fig.add_subplot(111)
        draw_map_background(m, ax)

        xx, yy = np.meshgrid(lons, lats)

        bins = np.array([0.0, 1.0, 5, 25, 250, 1000, math.inf])

        pop = np.digitize(pop, bins)
        pop = np.ma.array(pop, mask=mm)

        # mask values outside the plotting box
        x = np.outer(masklat, masklon)
        pop[~x] = ma.masked

        my_cmap = discrete_cmap(6, base_cmap=plt.get_cmap('YlGnBu'))

        my_cmap.set_bad(color='white', alpha=0)

        ax.text(ll_lon, ll_lat, y, fontsize=10, fontweight='bold', ha='left', va='bottom', color='k')

        im1 = m.pcolormesh(xx, yy, pop, vmin=0, vmax=np.max(pop), latlon=True, cmap=my_cmap)
        cbar = plt.colorbar(im1, orientation='horizontal')

        cbar.ax.get_xaxis().set_ticks([])
        for j, lab in enumerate(['< 1', ' < 5', '< 25', '< 250', '< 1000', '> 1000']):
            cbar.ax.text((2 * j + 1) / 12, -0.1, lab, ha='center', va='top')

        plt.savefig('./maps/example{0:4d}.png'.format(y), bbox_inches='tight')
        plt.close(fig)
        print('Done!')

    gif_name = './maps/outputName'
    os.system('convert -delay 45 -loop 0 ./maps/*.png {}.gif'.format(gif_name))


def powerplants_map():

    ll_lat = 29.65
    ll_lon = -94.143
    ur_lat = 38.77
    ur_lon = -76.56

    fname='~/Documents/CE/AreaTVACellsallCurtailEnvRegsCnoneSnor/genFleetInitial.csv'
    genfleet = pd.read_csv(filepath_or_buffer=os.path.expanduser(fname))

    df = genfleet.groupby(by='ORIS Plant Code').agg(({'Plant Name': 'first', 'PlantType': 'first',
                                                      'Capacity (MW)': 'sum',
                                                      'Latitude': 'mean', 'Longitude': 'mean'}))

    ptypes = np.unique(genfleet['PlantType'])
    shapestype = ["o", "v", "^", "<", ">", "s", "p", "P", "*", "+", "x", "D", "8"]
    df_shapes = pd.DataFrame({'ptypes': ptypes, 'shapes': shapestype})

    df = pd.merge(df, df_shapes, how='inner', left_on='PlantType', right_on='ptypes')

    m = Basemap(llcrnrlon=ll_lon, llcrnrlat=ll_lat, urcrnrlon=ur_lon, urcrnrlat=ur_lat, resolution='i', area_thresh=2500.)
    fig = plt.figure()
    ax = fig.add_subplot(111)
    draw_map_background(m, ax)

    list_series = []
    for i, pt in enumerate(ptypes):
        df_plot = df[df.PlantType == pt].reset_index(drop=True)
        sp = m.scatter(x=df_plot['Longitude'], y=df_plot['Latitude'], marker=df_plot.shapes[0], s=15,
                       label=pt, c='white', edgecolors='black')
        list_series.append(sp)

    plt.legend(handles=list_series, loc='upper center',
               bbox_to_anchor=(0.5, -0.05), ncol=4, fontsize='x-small')

    plt.savefig('./example.png')
    plt.close(fig)


def powerplants_map_2():
    rbmDataDir = '/Volumes/KIKO64/Data/DatafromUW/RBMRawWaterTData10Aug2016'
    rbmOutputDir = '/Volumes/KIKO64/Data/DatafromUW/RBMProcessedWaterT25Aug2016/bcc-csm1-1-m_rcp45_r1i1p1'
    spat = 'bcc-csm1-1-m_rcp45_r1i1p1'
    listcells = list(getAllGridCellLatLongsInSpatFile(rbmDataDir, spat))

    lat = np.array([t[0] for t in listcells])
    lon = np.array([t[1] for t in listcells])

    lat = np.sort(np.unique(lat))
    lon = np.sort(np.unique(lon))

    range_lat = (np.min(lat), np.max(lat))
    range_lon = (np.min(lon), np.max(lon))

    a = np.sort(np.unique(lon))
    delta = a[1] - a[0]

    lats_array = np.arange(start=range_lat[0], stop=range_lat[1], step=delta)
    lons_array = np.arange(start=range_lon[0], stop=range_lon[1], step=delta)

    xx, yy = np.meshgrid(lons_array, lats_array)

    values = ma.array(-1*np.ones(shape=xx.shape))

    locPrecision = 4
    for i, la in enumerate(lats_array):
        for j, lo in enumerate(lons_array):
            outputDir = os.path.join(rbmOutputDir, '{lat:.{p}f}_{long:.{p}f}'.format(p=locPrecision, lat=la, long=lo))
            fname = os.path.join(outputDir, createAverageTFilename(locPrecision, la, lo))
            if os.path.exists(outputDir):
                a = pd.read_csv(filepath_or_buffer=fname)
                values[i, j] = a['AverageWaterT(degC)'][0]

    empty_cells = np.array((values < 0))
    values[empty_cells] = ma.masked

    ll_lat = 29.65
    ll_lon = -94.143
    ur_lat = 38.77
    ur_lon = -76.56

    my_cmap = plt.get_cmap('YlGnBu')

    my_cmap.set_bad(color='white', alpha=0)

    m = Basemap(llcrnrlon=ll_lon, llcrnrlat=ll_lat, urcrnrlon=ur_lon, urcrnrlat=ur_lat,
                resolution='i', area_thresh=2500.)
    fig = plt.figure()
    ax = fig.add_subplot(111)
    draw_map_background(m, ax)

    im1 = m.pcolormesh(xx, yy, values, vmin=0, vmax=np.max(values), latlon=True, cmap=my_cmap)

    plt.savefig('./map.png')
    plt.close(fig)

