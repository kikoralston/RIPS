#Michael Craig, 18 May 2017
#Code from: http://stackoverflow.com/questions/42501349/given-latitude-longitude-say-if-the-coordinate-is-within-continental-us-or-not

from shapely.geometry import MultiPoint, Point, Polygon
import shapefile

#return a polygon for each state in a dictionary
def getStateBordersPolygon():
    sf = shapefile.Reader("C:\\Users\\mtcraig\\Desktop\\EPP Research\\Databases\\Shape files\\cb_2016_us_state_20m\\cb_2016_us_state_20m")
    shapes = sf.shapes()
    #shapes[i].points
    fields = sf.fields
    records = sf.records()
    state_polygons = {}
    for i, record in enumerate(records):
        state = record[5]
        points = shapes[i].points
        poly = Polygon(points)
        state_polygons[state] = poly
    return state_polygons

#Return state that a lat,lon is in
def inState(statePolygons,lat, lon):
    p = Point(lon, lat)
    for state in statePolygons:
        if statePolygons[state].contains(p):
            return state
    return None

#State polygons
statePolygons = getStateBordersPolygon()   
print(inState(statePolygons,30.7732,-84.9806))
print(inState(statePolygons,35.7058,-76.0342))
print(inState(statePolygons,34.7461,-89.884))
print(inState(statePolygons,30,-80))