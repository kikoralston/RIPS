import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from netCDF4 import Dataset
import os
import numpy as np
import numpy.ma as ma
import math
import time
from matplotlib import animation


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