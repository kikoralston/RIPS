# Michael Craig, 9 August 2017
# Determine what state a given lat/long is in

import os
from shapely.geometry import MultiPoint, Point, Polygon
import shapefile


# Load IPM zone shape file and return a dict mapping FIP # to IPM zone and
# a dict of FIP # to polygon. Note that FIP # is the only unique identifier in the
# IPM shape file.
def getStatePolys(dataRoot, states):

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


# Returns state that a given lat/lon is in
def getStateOfPt(stateToPoly, lat, lon):
    p = Point(lon, lat)
    for key in stateToPoly:
        if stateToPoly[key].contains(p):
            return key
    return None

# test = getStatePolys('pc',['Alabama','Georgia','South Carolina','North Carolina'])
# # print(getStateOfPt(test,37.523,-78.867))
# print(getStateOfPt(test,34.302,-81.922))
# print(getStateOfPt(test,35.481,-83.5055))
