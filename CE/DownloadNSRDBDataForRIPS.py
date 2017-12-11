#Michael Craig
#November 17, 2017
#Download NSRDB data through their API. Need to have NSRDB key in same directory
#as this code in a text file.

import os, time
import pandas as pd
from AuxFuncs import *

def master():
    #Set years to download (integer or tmy) (tmy = typical meteorology year)
    years = ['tmy']
    #Set directory to download to
    outputDir = 'C:\\Users\\mtcraig\\Desktop\\EPP Research\\Databases\\NSRDBRIPS'
    if not os.path.exists(outputDir): os.makedirs(outputDir)    
    allFiles = [f for f in os.listdir(outputDir) if os.path.isfile(os.path.join(outputDir, f))]
    #Whose API key to use
    whoseApiKey = 'michael' #michael or paulina
    #Create boxes via lat/lons that you will download from, and the lat/lon step
    #(note that NSRDB has much more data than what I download here!)
    latLonInt = 0.25 #degrees (.5 first, then .25)
    southBox = [(34.98,-90.01),(30.91,-81.73)] #upper left, lower right
    northBox = [(38.22,-89.16),(35.25,-77.16)]
    scBox = [(35.07,-81.44),(33.10,-79.51)]
    #Set up 2d list to hold metadata
    if 'nsrdbMetadata.csv' in [f for f in os.listdir(outputDir) if os.path.isfile(os.path.join(outputDir, f))]:
        elevations = readCSVto2dList(os.path.join(outputDir,'nsrdbMetadata.csv'))
    else: elevations = [['queriedLat','queriedLon','nsrdbLat','nsrdbLon','elevation(m)']]
    #Iterate through each box, save NSRDB data, & put metadata into new file
    for box in [southBox,northBox,scBox]:
        currLat,currLon = box[0][0],box[0][1]
        endLat,endLon = box[1][0],box[1][1]
        while currLat > endLat:
            while currLon < endLon:
                if createNSRDBFilename(currLat,currLon,years) not in allFiles:
                    elev,siteLat,siteLon = downloadLatLonNSRDBAPI(currLat,currLon,years,outputDir,whoseApiKey)
                    elevations.append([currLat,currLon,siteLat,siteLon,elev])
                    write2dListToCSV(elevations,os.path.join(outputDir,'nsrdbMetadata.csv'))
                currLon += latLonInt
            currLon = box[0][1] #reset longitude
            currLat -= latLonInt
        print('Finished box!')

#Download NSRDB data using API for given lat/lon for each year passed in 
#(can be integer or 'tmy' in list), and save file to outputDir
def downloadLatLonNSRDBAPI(lat,lon,years,outputDir,whoseApiKey):
    #Set attributes, whether want leap year data, mins of interval data (TMY is only
    #hourly), and whether convert to UTC (false means local tz)
    attributes = 'ghi,dhi,dni,wind_speed_10m_nwp,surface_air_temperature_nwp'
    leap_year = 'true' #CSI Interval data includes leap years
    interval = '60'
    utc = 'false'
    #Set other API inputs (use + not spaces)
    if whoseApiKey == 'michael':
        with open(os.path.join('nsrdbApiKey.txt')) as f: api_key = f.read()
        your_name = 'Michael+Craig'
        your_email = 'mtcraig@andrew.cmu.edu'
    elif whoseApiKey == 'paulina':
        with open(os.path.join('nsrdbApiKeyPauli.txt')) as f: api_key = f.read()
        your_name = 'Paulina+Jaramillo'
        your_email = 'paulina@cmu.edu'
    reason_for_use = 'research'
    your_affiliation = 'carnegie+mellon'
    mailing_list = 'true'
    #For each year, read in data, then concatenate dataframes
    print('Downloading point: ',(lat,lon))
    for idx in range(len(years)):
        year = years[idx]
        #Declare url string
        url = 'http://developer.nrel.gov/api/solar/nsrdb_0512_download.csv?wkt=POINT({lon}%20{lat})&names={year}&leap_day={leap}&interval={interval}&utc={utc}&full_name={name}&email={email}&affiliation={affiliation}&mailing_list={mailing_list}&reason={reason}&api_key={api}&attributes={attr}'.format(year=year, lat=lat, lon=lon, leap=leap_year, interval=interval, utc=utc, name=your_name, email=your_email, mailing_list=mailing_list, affiliation=your_affiliation, reason=reason_for_use, api=api_key, attr=attributes)
        #Get metadata (first 2 lines of file) and save lat & lon
        if idx == 0:
            info = pd.read_csv(url, nrows=1)
            siteLat,siteLon = round(info['Latitude'].values[0],2),round(info['Longitude'].values[0],2)
            siteElev = info['Elevation'].values[0]
        #If idx 0, initialize df; otherwise concatenate. Skip first 2 rows (metadata)
        df = pd.read_csv(url,skiprows=2)
        if (list(df.columns.values)!=['Year','Month','Day','Hour','Minute','GHI','DHI',
                                            'DNI','Wind Speed','Temperature']):
            print('Cols do not line up!')
        if idx == 0: dfAllYears = df.copy()
        else: dfAllYears = pd.concat([dfAllYears,df],ignore_index=True)
        time.sleep(5)
    dfAllYears.to_csv(os.path.join(outputDir,createNSRDBFilename(lat,lon,years)),index=False)
    return siteElev,siteLat,siteLon

#Create filename
def createNSRDBFilename(lat,lon,years):
    if 'tmy' in years: tmyFlag = 'tmy'
    else: tmyFlag = ''
    return 'nsrdb_'+str(round(lat,2))+'_'+str(round(lon,2))+tmyFlag+'.csv'

master()