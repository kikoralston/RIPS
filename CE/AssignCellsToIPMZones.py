# Michael Craig, July 13 2017
# Imports IPM shapefile and maps lat/lon values to IPM zones

from shapely.geometry import MultiPoint, Point, Polygon
import shapefile
import os
import numpy as np


# Load IPM zone shape file and return a dict mapping FIP # to IPM zone and
# a dict of FIP # to polygon. Note that FIP # is the only unique identifier in the
# IPM shape file.
def getIPMPolys(dataRoot, ipmZones):

    dataDir = os.path.join(dataRoot, "IPMRegionsShapeFile")

    sf = shapefile.Reader(os.path.join(dataDir, 'IPM_Regions_20121023_US'))
    shapes = sf.shapes()
    fields = sf.fields
    records = sf.records()
    # Get FIPS and zone idx in records
    for idx in range(1, len(fields)):
        if 'FIPS' in fields[idx]:
            fipsIdx = idx - 1  # need to subtract 1 b/c first entry in fields is metaentry
        elif 'IPM_Region' in fields[idx]:
            zoneIdx = idx - 1
    fipsToPolys, fipsToZones = dict(), dict()
    for i, record in enumerate(records):
        fips, zone = record[fipsIdx], record[zoneIdx]
        # if zone in ['S_SOU','S_VACA','S_C_KY','S_C_TVA']: #limit to SERC zones for computational speed
        if fips in fipsToZones: print('FIPS ' + fips + ' already in fipsToZones!')
        fipsToZones[fips] = zone
        points = shapes[i].points
        poly = Polygon(points)
        fipsToPolys[fips] = poly
    return fipsToZones, fipsToPolys


# Returns FIPS that a given lat/lon is in
def getFIPSOfPt(polygons, lat, lon):
    p = Point(lon, lat)
    for key in polygons:
        if polygons[key].contains(p):
            return key
    return None


# Inputs: set of cells formatted as lat_lon
# Outputs: dict mapping cells to IPM zones
def assignCellsToIPMZones(cells, fipsToZones, fipsToPolys):
    cellsToZones = dict()
    for cell in cells: cellsToZones[cell] = mapCellToIPMZone(cell, fipsToZones, fipsToPolys)
    return cellsToZones


# Input: 1 cell. Output: zone that cell belongs in.
def mapCellToIPMZone(cell, fipsToZones, fipsToPolys):
    (lat, lon) = cell.split('_')
    fips = getFIPSOfPt(fipsToPolys, float(lat), float(lon))

    if fips is None:
        out = 'NA'
    else:
        out = fipsToZones[fips]

    return out


def locInZone(lat, lon, tgtZone, fipsToZones, fipsToPolys):
    """Returns True if given lat/lon is in tgtZone (using IPM zones)

    :param lat:
    :param lon:
    :param tgtZone:
    :param fipsToZones:
    :param fipsToPolys:
    :return:
    """

    fips = getFIPSOfPt(fipsToPolys, lat, lon)

    if fips in fipsToZones:
        zone = fipsToZones[fips]
    else:
        zone = None  # several lat/lons on coast are not included in IPM zones

    return zone == tgtZone


def get_centroid_zone(fipsToPolys, fipsToZones, zone):

    fips_polys_zone = [fipsToPolys[p] for p in fipsToPolys if fipsToZones[p] == zone]

    lat_centroid = 0
    lon_centroid = 0

    n = len(fips_polys_zone)

    for p in fips_polys_zone:
        lon_centroid = lon_centroid + p.centroid.coords[0][0]/n
        lat_centroid = lat_centroid + p.centroid.coords[0][1]/n

    return np.round(lat_centroid, 4), np.round(lon_centroid, 4)


# # TEST CODE
# fipsToZones,fipsToPolys = getIPMPolys('pc','S_C_TVA')
# b = getFIPSOfPt(fipsToPolys,35.4375,-89.4375)
# c = getFIPSOfPt(fipsToPolys,35.6875,-89.4375)
# d = getFIPSOfPt(fipsToPolys,36.3125,-86.4375)
# print(fipsToZones[b],fipsToZones[c],fipsToZones[d])
# e = getFIPSOfPt(fipsToPolys,34.1875,-86.0625)
# print(e)
# print(fipsToZones[e])

# print(locInZone(35.05,-77.32,'S_VACA',fipsToZones,fipsToPolys))
# print(locInZone(35.05,-77.32,'S_SOU',fipsToZones,fipsToPolys))
# print(locInZone(35.608,-87.87,'S_C_TVA',fipsToZones,fipsToPolys))
