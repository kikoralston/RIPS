# Michael Craig, July 13 2017
# Imports IPM shapefile and maps lat/lon values to IPM zones

from shapely.geometry import MultiPoint, Point, Polygon
import shapefile
import os
import numpy as np


def getIPMPolys(dataRoot, ipmZones):
    """Get IPM zones polygons

    Load IPM zone shape file and return a dict mapping FIPS number to IPM zone name and a dict of FIPS number to
    polygon of counties. Note that FIP # is the only unique identifier in the IPM shape file.

    :param dataRoot: (string) full path to data folder
    :param ipmZones: (list) list with names of IPM zones included in study
    :return: tuple with two dictionaries
             * dict {fips: ipm zone name}
             * dict {fips: polygon}
    """
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


def getFIPSOfPt(polygons, lat, lon):
    """Get FIPS code of lat/lon point

    Returns the FIPS code that a given lat/lon is in

    :param polygons: (dict) dictionary {fips: polygon of counties}
    :param lat: (numeric) latitude
    :param lon: (numeric) longitude
    :return: (string) FIPS code of county
    """
    p = Point(lon, lat)
    for key in polygons:
        if polygons[key].contains(p):
            return key
    return None


# Inputs: set of cells formatted as lat_lon
# Outputs: dict mapping cells to IPM zones
def assignCellsToIPMZones(cells, fipsToZones, fipsToPolys):
    """Map all grid cells to IPM zones

    :param cells: (list) list with grid cells. Each element of the list is a string in the format '{lat}_{long}'
    :param fipsToZones: dict {fips: ipm zone name}
    :param fipsToPolys: dict {fips: polygon}
    :return: (dict) {grid cell: IPM zone}
    """
    cellsToZones = dict()
    for cell in cells: cellsToZones[cell] = mapCellToIPMZone(cell, fipsToZones, fipsToPolys)
    return cellsToZones


# Input: 1 cell. Output: zone that cell belongs in.
def mapCellToIPMZone(cell, fipsToZones, fipsToPolys):
    """Map one single grid cell to IPM zones

    :param cell: (string) string with grid cell in the format '{lat}_{long}'
    :param fipsToZones: dict {fips: ipm zone name}
    :param fipsToPolys: dict {fips: polygon}
    :return: (string) ipm zone name that cell belongs in.
    """
    (lat, lon) = cell.split('_')
    fips = getFIPSOfPt(fipsToPolys, float(lat), float(lon))

    if fips is None:
        out = 'NA'
    else:
        out = fipsToZones[fips]

    return out


def locInZone(lat, lon, tgtZone, fipsToZones, fipsToPolys):
    """Check if given cell belongs to IPM zone

    :param lat: (numeric) latitude
    :param lon: (numeric) longitude
    :param tgtZone: (string) name of ipm zone
    :param fipsToZones: dict {fips: ipm zone name}
    :param fipsToPolys: dict {fips: polygon}
    :return: (boolean) Returns True if given lat/lon is in tgtZone (using IPM zones)
    """

    fips = getFIPSOfPt(fipsToPolys, lat, lon)

    if fips in fipsToZones:
        zone = fipsToZones[fips]
    else:
        zone = None  # several lat/lons on coast are not included in IPM zones

    return zone == tgtZone


def get_centroid_zone(fipsToPolys, fipsToZones, zone):
    """Get centroid of ipm zone

    :param fipsToZones: dict {fips: ipm zone name}
    :param fipsToPolys: dict {fips: polygon}
    :param zone: (string) name of ipm zone
    :return: (tuple) latitude and longitude of centroid (numeric values)
    """
    fips_polys_zone = [fipsToPolys[p] for p in fipsToPolys if fipsToZones[p] == zone]

    lat_centroid = 0
    lon_centroid = 0

    n = len(fips_polys_zone)

    for p in fips_polys_zone:
        lon_centroid = lon_centroid + p.centroid.coords[0][0]/n
        lat_centroid = lat_centroid + p.centroid.coords[0][1]/n

    return np.round(lat_centroid, 4), np.round(lon_centroid, 4)


def getStatePolys(dataRoot, states):
    """

    Load IPM zone shape file and return a dict mapping FIPS number to states and a dict of FIPS number to states polygon.
    Note that FIPS number is the only unique identifier in the IPM shape file.

    :param dataRoot: (string) full path to data folder
    :param states: (list) list of states
    :return: (dict) dictionary mapping state name to polygon
    """
    dataDir = os.path.join(dataRoot, "StateShapeFiles")

    sf = shapefile.Reader(os.path.join(dataDir, 'cb_2016_us_state_20m.shp'))
    shapes = sf.shapes()
    fields = sf.fields
    records = sf.records()
    # Get state idx in records
    for idx in range(1, len(fields)):
        if 'NAME' in fields[idx]: nameIdx = idx - 1  # need to subtract 1 b/c first entry in fields is metaentry
    # Map states to polys
    stateToPoly = dict()
    for i, record in enumerate(records):
        state = record[nameIdx]
        if state in states:
            points = shapes[i].points
            poly = Polygon(points)
            stateToPoly[state] = poly
    return stateToPoly


def getStateOfPt(stateToPoly, lat, lon):
    """Get name of state that lat/long point is at

    :param stateToPoly: (dict) dictionary mapping state name to polygon
    :param lat: (numeric) latitude
    :param lon: (numeric) longitude
    :return: (string) name of state that lat/long point is at
    """
    p = Point(lon, lat)
    for key in stateToPoly:
        if stateToPoly[key].contains(p):
            return key
    return None

# test = getStatePolys('pc',['Alabama','Georgia','South Carolina','North Carolina'])
# # print(getStateOfPt(test,37.523,-78.867))
# print(getStateOfPt(test,34.302,-81.922))
# print(getStateOfPt(test,35.481,-83.5055))


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
