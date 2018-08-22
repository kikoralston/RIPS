#Michael Craig, July 13 2017
#Imports

from shapely.geometry import MultiPoint, Point, Polygon
import shapefile
import pickle as pk

from AuxFuncs import *

#Load IPM zone shape file and return a dict mapping FIP # to IPM zone and
#a dict of FIP # to polygon. Note that FIP # is the only unique identifier in the
#IPM shape file.
def getIPMPolys():
    sf = shapefile.Reader("C:\\Users\\mtcraig\\Desktop\\EPP Research\\Databases\\IPMRegionsShapeFile\\IPM_Regions_20121023_US")
    shapes = sf.shapes()
    fields = sf.fields
    records = sf.records()
    #Get FIPS and zone idx in records
    for idx in range(1,len(fields)):
        if 'FIPS' in fields[idx]: fipsIdx = idx-1 #need to subtract 1 b/c first entry in fields is metaentry
        elif 'IPM_Region' in fields[idx]: zoneIdx = idx-1
    fipsToPolys,fipsToZones = dict(),dict()
    for i, record in enumerate(records):
        fips,zone = record[fipsIdx],record[zoneIdx]
        if fips in fipsToZones: print('FIPS ' + fips + ' already in fipsToZones!')
        fipsToZones[fips] = zone
        points = shapes[i].points
        poly = Polygon(points)
        fipsToPolys[fips] = poly
    return fipsToZones,fipsToPolys

#Returns FIPS that a given lat/lon is in
def getFIPSOfPt(polygons, lat, lon):
    p = Point(lon, lat)
    for key in polygons:
        if polygons[key].contains(p):
            return key
    return None

#Inputs: set of cells formatted as lat_lon
#Outputs: dict mapping cells to IPM zones
def assignCellsToZones(cells):
    fipsToZones,fipsToPolys = getIPMPolys()
    cellsToZones = dict()
    for cell in cells:
        (lat,lon) = cell.split('_')
        fips = getFIPSOfPt(fipsToPolys,lat,lon)
        cellsToZones[cell] = fipsToZones[fips]
    return cellsToZones


def create_dict_cells_ipm(file_name):
    """
    Creates a dictionary mapping lat long grod cell to a IPM zone

    This function reads a netcdf file with the full set of grid cells being used in the SERC study and assigns
    an IPM zone to each of them. It saves the results to a dictionary and saves this dictionary to a pickle file so
    it can be read and used in the future. The objective is to save time when running the full SERC simulation

    :param file_name: string with path to pickle file to save dictionary
    :return: nothing
    """

    genparam = Generalparameters.Generalparameters()
    genparam.load(fname='./generalparameters.txt')

    curtailparam = Curtailmentparameters.Curtailmentparameters()
    curtailparam.load(fname='./curtailmentparameters.txt')

    dataset = nc.Dataset('/Users/kiko/Downloads/serc.NorESM1-M.RCP85.stream_T.nc')

    # Extract data from NetCDF file
    lats = dataset.variables['lat'][:]
    lons = dataset.variables['lon'][:]

    cells2zones = dict()

    t_year = time.time()
    for i, la in enumerate(lats):
        for j, lo in enumerate(lons):
            cell = '{}_{}'.format(la, lo)
            print('{}_{}'.format(la, lo))
            if getFIPSOfPt(genparam.fipsToPolys, la, lo) is not None:
                z = genparam.fipsToZones[getFIPSOfPt(genparam.fipsToPolys, la, lo)]
                cells2zones[cell] = z
            else:
                # some cells in the original grid may be in the sea, so they would have no zone assigned to them
                cells2zones[cell] = None
    print('Elapsed Time: ' + str_elapsedtime(t_year))

    with open(file_name, 'wb') as f:
        pk.dump(obj=cells2zones, file=f)


def test_read_dict(file_name):

    with open(file_name, 'rb') as f:
        cells2zones = pk.load(f)

    t_start = time.time()
    for i, c in enumerate(cells2zones):
        print('{}. ({}) : {}'.format(i+1, c, cells2zones[c]))
    print('Elapsed Time: ' + str_elapsedtime(t_start))


#TEST CODE
# fipsToZones,fipsToPolys = getIPMPolys()
# b = getFIPSOfPt(fipsToPolys,33.821,-86.976)
# c = getFIPSOfPt(fipsToPolys,36.423,-85.423)
# d = getFIPSOfPt(fipsToPolys,35.022,-77.54)
# print(fipsToZones[b],fipsToZones[c],fipsToZones[d])

