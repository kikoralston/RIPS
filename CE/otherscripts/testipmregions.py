#Michael Craig, July 13 2017
#Imports 

from shapely.geometry import MultiPoint, Point, Polygon
import shapefile

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


#TEST CODE
# fipsToZones,fipsToPolys = getIPMPolys()
# b = getFIPSOfPt(fipsToPolys,33.821,-86.976)
# c = getFIPSOfPt(fipsToPolys,36.423,-85.423)
# d = getFIPSOfPt(fipsToPolys,35.022,-77.54)
# print(fipsToZones[b],fipsToZones[c],fipsToZones[d])

