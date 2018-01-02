import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap


def draw_map_background(m, ax):
    ax.set_facecolor('#729FCF')
    m.fillcontinents(color='#FAFAFA', ax=ax, zorder=0, lake_color='#729FCF')
    m.drawcounties(ax=ax, color='#c5c5c5')
    m.drawrivers(ax=ax, color='#729FCF')
    m.drawstates(ax=ax)
    m.drawcountries(ax=ax)
    m.drawcoastlines(ax=ax)


ll_lat = 29.65
ll_lon = -94.143
ur_lat = 38.77
ur_lon = -76.56
m = Basemap(llcrnrlon=ll_lon, llcrnrlat=ll_lat, urcrnrlon=ur_lon, urcrnrlat=ur_lat,
            resolution='i', area_thresh=2500.)
fig = plt.figure()
ax = fig.add_subplot(111)
draw_map_background(m, ax)

plt.savefig('./maps/example.png', bbox_inches='tight')
plt.close(fig)