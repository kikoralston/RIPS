#Michael Craig, 18 May 2017
#This script contains several functions for extracting and processing
#WIND actual power output data.
#Power output is given for 5-minute intervals from 2007 through 2012, including
#leap year days in 2008 and 2012. Timezone is UTC.

from netCDF4 import Dataset
import numpy as np
import datetime as dt
import os
from AuxFuncs import * 

#Create 1d list of datetimes for period that WIND data is for (2007-2012). Note
#taht WIND data timezone is UTC, although I don't define datetimes here as that.
def createDts():
    startYear = 2007
    endYear = 2013
    leapDays = 2
    sampleMins = 5 #sampled every 5 mins
    startDt = dt.datetime(startYear,1,1,0,0)
    sampleDts = [startDt + dt.timedelta(minutes=x) for x in range(0,(endYear-startYear)*8760*60+2*24*60,sampleMins)]
    hourlyDts = [startDt + dt.timedelta(hours=x) for x in range((endYear-startYear)*8760+2*24)]
    return sampleDts,hourlyDts

#Extracts lat, lon, and power output data from given CDF file.
def extractCDFData(file):
    data = Dataset(file,mode='r')
    print(data)
    lat,lon = data.longitude,data.latitude
    power = data.variables['power'][:]
    return lat,lon,power #float, float, np array



################ PROCESS SERC SITE POWER GENERATION DATA #######################
#This code block isolates SERC metadata and saves it to a file; calculates
#and saves hourly power output for all years for each SERC site & saves
#that data to separate files; and calculates annual CFs for each plant.
def calculateAnnualCFsSERC(states):
    metadata = importSiteMetadata()
    sercMetadata = getSitesInStates(states,metadata)
    write2dListToCSV(sercMetadata,'E:\\toolkit_sites_v7_SERC.csv')
    sampleDts,hourlyDts = createDts()
    siteCol,capacCol = sercMetadata[0].index('site_id'),sercMetadata[0].index('capacity')
    siteAnnualCfs = [['Site'] + [str(year) +'CF' for year in range(2007,2013)]]
    for site in sercMetadata[1:]:
        siteId,siteCapac = int(site[siteCol]),float(site[capacCol])
        lat,lon,power = extractCDFData(os.path.join('E:','met_data',str(siteId//500),str(siteId)+'.nc'))
        # np.savetxt('E:\\test' + str(siteId) + '.csv',power)
        powerHourly = averagePowerToHourly(power,siteId,hourlyDts)
        write2dListToCSV(powerHourly2d,os.path.join('E:','hourlyPowerSERC','powerhourly_' + str(siteId) + '.csv'))
        cfs = [val/siteCapac for val in powerHourly]
        # np.savetxt('E:\\cfstest' + str(siteId) + '.csv',np.array(cfs))
        annualCfs = averageCfsAnnually(cfs,hourlyDts)
        siteAnnualCfs.append([siteId] + annualCfs)
    write2dListToCSV(siteAnnualCfs,'sercAnnualCfs.csv')

#Import metadata for all of WIND dataset
def importSiteMetadata():
    return readCSVto2dList('E:\\toolkit_sites_v7.csv')

#Isolates metadata for given list of states
def getSitesInStates(coreSERCStates,metadata):
    heads = metadata[0]
    stateCol = heads.index('State')
    return [heads] + [row for row in metadata[1:] if row[stateCol] in coreSERCStates]

#Averages 5-min power output for a given site to hourly values, and returns
#2d list of [datetime,hourly power output] w/ headers..
def averagePowerToHourly(power,siteId,hourlyDts):
    sampleMins = 5 #power in 5 min time steps
    stepsPerHour = 60//5
    # assert(len(powerHourly)==8760*6+48)
    powerHourly = [np.average(power[idx:idx+stepsPerHour]) for idx in range(0,len(power),stepsPerHour)]
    #For writing purposes
    powerHourly2d = [['Datetime',siteId]]
    for idx in range(len(powerHourly)): powerHourly2d.append([hourlyDts[idx],powerHourly[idx]])
    return powerHourly

#Averages hourly CFs to annual level.
def averageCfsAnnually(cfs,hourlyDts):
    annualCfs = list()
    for currYear in range(2007,2013):
        yearIdxs = [idx for idx in range(len(hourlyDts)) if hourlyDts[idx].year==currYear]
        yearCfs = [cfs[idx] for idx in yearIdxs]
        annualCfs.append(sum(yearCfs)/len(yearCfs))
    return annualCfs

#Next few lines write hourly power output data for SERC states to folder.
# coreSERCStates = {'North Carolina','South Carolina','Georgia','Mississippi','Alabama',
#         'Kentucky','Tennessee'}
# calculateAnnualCFsSERC(coreSERCStates)
################################################################################

########### GET POWER OUTPUT FOR GIVEN SITE AND YEAR IN SITE'S LOCAL TIMEZONE ##
#This function takes a particular site; metadata (for just SERC) region;
#year for which WIND data is extracted; and hourly DTs over timespan of WIND data.
#Outputs: 1d list of hourly power output for year in local timezone.
def getSitePowerInActualTzForYear(siteId,metadata,dataYear,hourlyDts):
    TZOFFSETDICT = {'UTCtoCST':-6,'CSTtoCST':0,'ESTtoCST':-1,'CSTtoEST':1,'UTCtoEST':-5}
    #Get timezone & time offset
    lat,lon,state = getSiteLoc(siteId,metadata)
    siteTz = tzOfLatLon(lat,lon,state)
    hourOffset = TZOFFSETDICT['UTCto' + siteTz]
    #Change hourly DTs & get idxs in year
    hourlyDtsTz = [hourlyDt + dt.timedelta(hours=hourOffset) for hourlyDt in hourlyDts]
    yearIdxs = [idx for idx in range(len(hourlyDtsTz)) if hourlyDtsTz[idx].year == dataYear]
    #Import data and return data for given year
    allData = readCSVto2dList(os.path.join('E:','hourlyPowerSERC','powerhourly_' + str(siteId) + '.csv'))
    powCol = allData[0].index(str(siteId))
    hourlyPow = [float(row[powCol]) for row in allData[1:]]
    return [hourlyPow[idx] for idx in yearIdxs]

#Returns lat, lon, and state of given site
def getSiteLoc(siteId,metadata):
    siteCol,latCol = metadata[0].index('site_id'),metadata[0].index('latitude')
    lonCol,stateCol = metadata[0].index('longitude'),metadata[0].index('State')
    siteIds = [row[siteCol] for row in metadata]
    siteIdx = siteIds.index(siteId)
    siteRow = metadata[siteIdx]
    return siteRow[latCol],siteRow[lonCol],siteRow[stateCol]

#Returns timezone (EST or CST) of given lat, lon, and state position for SERC.
def tzOfLatLon(lat,lon,state):
    STATETIMEZONES = {'North Carolina':'EST','South Carolina':'EST','Virginia':'EST',
                  'Georgia':'EST','Mississippi':'CST','Alabama':'CST','Louisiana':'CST',
                  'Missouri':'CST','Arkansas':'CST','Illinois':'CST','Kentucky':'CSTorEST',
                  'Tennessee':'CSTorEST','Texas':'CST'}
    if STATETIMEZONES[state] == 'CST' or STATETIMEZONES[state] == 'EST':
        tz = STATETIMEZONES[state]
    else: #in KY or TN, which is half & half CST & EST
        kyLine = [(36.601261, -84.861318),(38.048865, -86.251172)]
        tnLine = [(36.601261, -85.076648),(34.997791, -85.605184)]
        if state == 'Tennessee': line = tnLine
        elif state == 'Kentucky': line = kyLine
        if siteEastOfLine(line,float(lat),float(lon)): tz = 'EST'
        else: tz = 'CST'
    return tz

#Determines whether lat,lon is in EST or CST for KY & TN, which include 2 TZs.
def siteEastOfLine(line,siteLat,siteLong):
    (deltaLat,deltaLong) = (line[0][0]-line[1][0],line[0][1]-line[1][1])
    lineSlope = deltaLat/deltaLong
    intercept = line[1][0] - lineSlope * line[1][1] #b = y - mx
    longOnLineForSiteLat = (siteLat-intercept)/lineSlope #x = (y-b)/m
    return siteLong > longOnLineForSiteLat #long decreases (more negative) west across US

#Test this codeblock
# sercMetadata = readCSVto2dList('E:\\toolkit_sites_v7_SERC.csv')
# sampleDts,hourlyDts = createDts()
# powOutput = getSitePowerInActualTzForYear('8859',sercMetadata,2011,hourlyDts)
# write2dListToCSV([powOutput],'E:\\test8859.csv')
################################################################################
