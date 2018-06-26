import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import netCDF4 as nc
import os
import numpy as np
import numpy.ma as ma
import math
import pandas as pd
import time
import datetime as dt
from matplotlib import animation
import sys
import time

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

    fname='/Users/kiko/Documents/CE/AreaTVACellsallCurtailEnvRegsCnoneSnor/genFleetInitial.csv'
    genfleet = pd.read_csv(filepath_or_buffer=os.path.expanduser(fname))

    df = genfleet.groupby(by='ORIS Plant Code').agg(({'Plant Name': 'first', 'PlantType': 'first',
                                                      'Capacity (MW)': 'sum',
                                                      'Latitude': 'mean', 'Longitude': 'mean'}))
    df = df.replace(to_replace={'PlantType': {'Biomass|Landfill Gas|Non-Fossil Waste|O/G Steam': 'Other'}}, regex=True)

    df = df.replace(to_replace={'PlantType': {'IGCC|Combined Cycle': 'Combined Cycle'}}, regex=True)

    ptypes = np.unique(df['PlantType'])
    shapestype = ["o", "v", "^", "<", ">", "s", "D", "H", "P"]

    df_shapes = pd.DataFrame({'ptypes': ptypes, 'shapes': shapestype[:len(ptypes)]})

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

    plt.savefig('./example.png', bbox_inches='tight')
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


def read_netcdf_data():

    path_data = '/Users/kiko/Downloads/bcc-csm1-1-m'

    # Extract data from NetCDF file
    #lats = dataset.variables['latitude'][:]
    #lons = dataset.variables['longitude'][:]
    #time = dataset.variables['raster'][:]
    #values = dataset.variables['Population Count, v4.10 (2000, 2005, 2010, 2015, 2020): 2.5 arc-minutes'][:]

    years = np.arange(2015, 2051, step=5)

    #ll_lat = 24
    #ll_lon = -125
    #ur_lat = 50
    #ur_lon = -60

    ll_lat = 29.65
    ll_lon = -94.143
    ur_lat = 38.77
    ur_lon = -76.56

    for i, y in enumerate(years):
        print('Reading data for year {0:4d}'.format(y))
        dataset = nc.Dataset(os.path.join(path_data, 'forcing_maca_bcc-csm1-1-m_{0:4d}.nc'.format(y)))

        # Extract data from NetCDF file
        lats = dataset.variables['lat'][:]
        lons = dataset.variables['lon'][:]

        temp = dataset.variables['temp'][:]
        air_pressure = dataset.variables['air_pressure'][:]
        rel_humid = dataset.variables['rel_humid'][:]

        # create array with hourly time stamps for year
        start = dt.datetime(y, 1, 1)
        end = dt.datetime(y, 12, 31, 23, 00, 00)
        date_array = pd.date_range(start=start, end=end, freq='H')

        for m in np.arange(1, 13):
            print('Creating plot for month {0:2d}...'.format(m), end='', flush=True)
            # compute monthly temp
            mean_temp = np.mean(temp[date_array.month == m, :, :], axis=0)

            bm = Basemap(llcrnrlon=ll_lon, llcrnrlat=ll_lat, urcrnrlon=ur_lon, urcrnrlat=ur_lat, resolution='i', area_thresh=2500.)
            fig = plt.figure()
            ax = fig.add_subplot(111)
            draw_map_background(bm, ax)

            xx, yy = np.meshgrid(lons, lats)

            my_cmap = plt.get_cmap('RdYlBu_r')

            my_cmap.set_bad(color='white', alpha=0)

            ax.text(ll_lon, ll_lat, '{0:02d}-{1:04d}'.format(m, y), fontsize=10, fontweight='bold', ha='left',
                    va='bottom', color='k')
            im1 = bm.pcolormesh(xx, yy, mean_temp, vmin=-20, vmax=50, latlon=True, cmap=my_cmap)
            cbar = plt.colorbar(im1, orientation='horizontal')

            plt.savefig(os.path.expanduser('~/Documents/maps/example{1:4d}_{0:02d}.png'.format(m, y)),
                        bbox_inches='tight')
            plt.close(fig)
            print('Done!')

    gif_name = os.path.expanduser('~/Documents/maps/outputName')
    os.system('convert -delay 45 -loop 0 ~/Documents/maps/*.png {}.gif'.format(gif_name))

#    time = dataset.variables['time'][:]
#
#    if len(time) != len(date_array):
#        print('Netcdf files have data for {0:4d} hours, but year {1:4d} has {2:4d} hours. '
#              'Check data sources'.format(len(time), year, len(date_array)))


def process_netcdf(cellLat, cellLon):

    t1 = time.time()

    path_data = '/Users/kiko/Downloads/bcc-csm1-1-m'
    years = np.arange(2015, 2051, step=5)

    for i, y in enumerate(years):
        print('Reading data for year {0:4d}'.format(y))
        dataset = nc.Dataset(os.path.join(path_data, 'forcing_maca_bcc-csm1-1-m_{0:4d}.nc'.format(y)))

        # Extract data from NetCDF file
        lats = dataset.variables['lat'][:]
        lons = dataset.variables['lon'][:]
        temp = dataset.variables['temp'][:]
        air_pressure = dataset.variables['air_pressure'][:]
        rel_humid = dataset.variables['rel_humid'][:]

        ix = np.argwhere(lats == cellLat).flatten()[0]
        iy = np.argwhere(lons == cellLon).flatten()[0]

        # create array with hourly time stamps for year
        start = dt.datetime(y, 1, 1)
        end = dt.datetime(y, 12, 31, 23, 00, 00)
        date_array = pd.date_range(start=start, end=end, freq='H')

        a = pd.DataFrame({'date': date_array, 'temp': temp[:, ix, iy], 'air_pressure': air_pressure[:, ix, iy],
                          'rel_humid': rel_humid[:, ix, iy]}, columns=['date', 'temp', 'rel_humid', 'air_pressure'])

        if i == 0:
            df_out = a
        else:
            df_out = pd.concat([df_out, a])

    print('Elapsed time: {0:10d} s'.format(int(time.time() - t1)))

#    time = dataset.variables['time'][:]
#
#    if len(time) != len(date_array):
#        print('Netcdf files have data for {0:4d} hours, but year {1:4d} has {2:4d} hours. '
#              'Check data sources'.format(len(time), year, len(date_array)))


def curtailment_map(pathin, pathout):
    #pathin='/Users/kiko/Documents/CE/test.nc'
    #pathout='/Users/kiko/Documents/maps/'

    dataset = nc.Dataset(os.path.expanduser(pathin))

    # Extract data from NetCDF file
    lats = dataset.variables['latitude'][:]
    lons = dataset.variables['longitude'][:]
    time = dataset.variables['time'][:]
    values = dataset.variables['Coal Steam+OT'][:]

    # ll_lat = 24
    # ll_lon = -125
    # ur_lat = 50
    # ur_lon = -60

    ll_lat = 29.65
    ll_lon = -94.143
    ur_lat = 38.77
    ur_lon = -76.56

    year = 2020

    start = dt.datetime(year, 1, 1)
    end = dt.datetime(year, 12, 31, 23, 00, 00)
    date_array = pd.date_range(start=start, end=end, freq='H')

    start = dt.datetime(year, 1, 1)
    end = dt.datetime(year, 12, 31)
    array_days = pd.date_range(start=start, end=end, freq='D')

    for i, d in enumerate(array_days):
        print('Creating plot for {0:04d}/{1:02d}/{2:02d}'.format(d.year, d.month, d.day))

        # get all hours for day d and compute daily mean for each cell
        irows = np.where(date_array.strftime('%Y-%m-%d') == d.strftime('%Y-%m-%d'))[0]
        cap = values.data[irows, :, :]
        cap = np.mean(cap, axis=0)
        cap = ma.array(cap, mask=values.mask[0, :, :])

        # get mask
        mm = cap.mask

        masklon = (lons > ll_lon) & (lons < ur_lon)
        masklat = (lats > ll_lat) & (lats < ur_lat)

        m = Basemap(llcrnrlon=ll_lon, llcrnrlat=ll_lat, urcrnrlon=ur_lon, urcrnrlat=ur_lat,
                    resolution='i', area_thresh=2500.)

        fig = plt.figure()
        ax = fig.add_subplot(111)
        draw_map_background(m, ax)

        xx, yy = np.meshgrid(lons, lats)

        bins = np.array([0, 65, 130, 195, 260, 325, 390, 455, 520, 585, 650])

        cap = np.digitize(cap, bins)
        cap = np.ma.array(cap, mask=mm)

        # mask values outside the plotting box
        x = np.outer(masklat, masklon)
        cap[~x] = ma.masked

        my_cmap = discrete_cmap(10, base_cmap=plt.get_cmap('RdBu'))

        my_cmap.set_bad(color='white', alpha=0)

        label_date = '{0:04d}/{1:02d}/{2:02d}'.format(d.year, d.month, d.day)

        ax.text(ll_lon, ll_lat, label_date, fontsize=10, fontweight='bold', ha='left', va='bottom', color='k')

        im1 = m.pcolormesh(xx, yy, cap, vmin=1, vmax=11, latlon=True, cmap=my_cmap)

        cbar = plt.colorbar(im1, orientation='horizontal')

        plt.savefig(os.path.join(os.path.expanduser(pathout), 'example{0:4d}.png'.format(i)), bbox_inches='tight')
        plt.close(fig)
        print('Done!')

    gif_name = os.path.join(os.path.expanduser(pathout), 'outputName')
    os.system('convert -delay 45 -loop 0 {0}/*.png {1}.gif'.format(pathout, gif_name))

