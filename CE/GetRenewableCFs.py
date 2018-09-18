# Michael Craig, 6 July 2016
# Get wind and solar capacity factors

# This script gets wind and solar capacity factors for existing and new renewables
# by mapping capacity to the best available units in each zone.

import os, copy, datetime
from statistics import mode
from AuxFuncs import *
from AssignCellsToIPMZones import locInZone, get_centroid_zone


def getRenewableCFs(genFleet, startWindCapacForCFs, startSolarCapacForCFs, desiredTz,
                    dataRoot, windGenDataYr, currZone, fipsToZones, fipsToPolys, subHour=False):
    """GET RENEWABLE CAPACITY FACTORS

    Function takes in generator fleet that has already been isolated to a single zone and maps data to CFs for
    RE sites in same zone.

    :param genFleet:
    :param startWindCapacForCFs:
    :param startSolarCapacForCFs:
    :param desiredTz:
    :param dataRoot:
    :param windGenDataYr:
    :param currZone:
    :param fipsToZones:
    :param fipsToPolys:
    :param subHour: if True reads sub hourly data for Wind Data
    :return:
    """
    # Isolate wind & solar units
    plantTypeCol = genFleet[0].index('PlantType')
    windUnits = [genFleet[0]] + [row for row in genFleet if row[plantTypeCol] == 'Wind']
    solarUnits = [genFleet[0]] + [row for row in genFleet if row[plantTypeCol] == 'Solar PV']

    # Get list of wind / solar sites in region

    windDir = os.path.join(dataRoot, 'WINDSERCData')
    solarDir = os.path.join(dataRoot, 'NRELSolarPVData', 'SERC')

    if len(windUnits) > 1:  # have wind in zone
        (windCFs, windCfsDtHr, windCfsDtSubhr, ewdIdAndCapac) = getWindCFs(windUnits, windDir, startWindCapacForCFs,
                                                                           desiredTz, windGenDataYr, fipsToZones,
                                                                           fipsToPolys, currZone, subHour=subHour)
        windCfsDtHr, windCfsDtSubhr = rotate(windCfsDtHr), rotate(windCfsDtSubhr)
    else:
        windCFs, windCfsDtHr, windCfsDtSubhr, ewdIdAndCapac = None, None, None, None

    if len(solarUnits) > 1:
        # solar data has only sub hour data
        (solarCFs, solarCfsDtHr, solarCfsDtSubhr, solarFilenameAndCapac) = getSolarCFs(solarUnits, solarDir,
                                                                                       startSolarCapacForCFs, desiredTz,
                                                                                       fipsToZones, fipsToPolys,
                                                                                       currZone)
        solarCfsDtHr, solarCfsDtSubhr = rotate(solarCfsDtHr), rotate(solarCfsDtSubhr)
    else:
        solarCFs, solarCfsDtHr, solarCfsDtSubhr, solarFilenameAndCapac = None, None, None, None

    return (windCFs, windCfsDtHr, windCfsDtSubhr, ewdIdAndCapac, solarCFs, solarCfsDtHr, solarCfsDtSubhr,
            solarFilenameAndCapac)


def getWindCFs(windUnits, windDir, startWindCapacForCFs, desiredTz, windGenDataYr, fipsToZones, fipsToPolys, currZone,
               subHour=False):
    """GET WIND Capacity FACTORS

    All CFs output in EST

    :param windUnits:
    :param windDir:
    :param startWindCapacForCFs:
    :param desiredTz:
    :param windGenDataYr:
    :param fipsToZones:
    :param fipsToPolys:
    :param currZone:
    :param subHour: if True reads wind sub hourly data
    :return:
    """
    # Get total wind capacity by zone
    capacCol = windUnits[0].index('Capacity (MW)')
    windCapacInZone = sum([float(row[capacCol]) for row in windUnits[1:]])

    # Get wind plants per state
    (ewdIdAndCapac, ewdMetadata) = getBestWindIdsInZone(windDir, windCapacInZone, startWindCapacForCFs,
                                                        fipsToZones, fipsToPolys, currZone)
    # Import CFs for each wind plant
    if subHour:
        allSiteCfsHourly, allSiteCfsSubhourly, avgFleetCfHr = [], [], []
    else:
        allSiteCfsHourly, allSiteCfsSubhourly, avgFleetCfHr = [], None, []

    idCol = ewdIdAndCapac[0].index('Id')
    datasetCapacCol = ewdIdAndCapac[0].index('DatasetCapacity')

    for site in ewdIdAndCapac[1:]:
        (siteId, datasetCapac) = (site[idCol], site[datasetCapacCol])
        if 'NoMoreSites' in siteId:
            print('No more wind sites!')
            if avgFleetCfHr == []:
                # If doing new RE CFs & existing RE > potential RE capac for CFs,
                # can end up w/ nothing in allSiteCFs. In that case, new RE CFs
                # should jsut be fleet average CF.
                if allSiteCfsHourly == []:
                    a, tempAllSiteCfsHourly, tempAllSiteCfsSubhourly, tempEwdIdAndCapac = getWindCFs(windUnits,
                                                                                                     windDir, 0,
                                                                                                     desiredTz,
                                                                                                     windGenDataYr,
                                                                                                     fipsToZones,
                                                                                                     fipsToPolys,
                                                                                                     currZone,
                                                                                                     subHour=subHour)
                    avgFleetCfHr, avgFleetCfSubhr = calcCapacWtdFleetCf(tempEwdIdAndCapac, tempAllSiteCfsHourly,
                                                                        tempAllSiteCfsSubhourly)
                else:
                    avgFleetCfHr, avgFleetCfSubhr = calcCapacWtdFleetCf(ewdIdAndCapac, allSiteCfsHourly,
                                                                        allSiteCfsSubhourly)

            siteCfsHourly, siteCfsSubhourly = copy.deepcopy(avgFleetCfHr), copy.deepcopy(avgFleetCfSubhr)

        else:
            siteCfsHourly, siteCfsSubhourly = getWindSiteCfs(windDir, siteId, datasetCapac, desiredTz, windGenDataYr,
                                                             subHour=subHour)

        addSiteCfsToAggList(siteCfsHourly, siteCfsSubhourly, allSiteCfsHourly, allSiteCfsSubhourly, siteId)

    allSiteCfsHourOfYear = [['HourOfYear'] + [val for val in range(1, 8761)]] + copy.deepcopy(allSiteCfsHourly[1:])

    return allSiteCfsHourOfYear, allSiteCfsHourly, allSiteCfsSubhourly, ewdIdAndCapac


def addSiteCfsToAggList(siteCfsHourly, siteCfsSubhourly, allSiteCfsHourly, allSiteCfsSubhourly, siteId, subHour=False):

    if allSiteCfsHourly == []:
        allSiteCfsHourly.append(siteCfsHourly[0])
        allSiteCfsHourly.append([siteId] + siteCfsHourly[1][1:])  # replace header w/ site ID
        if subHour:
            allSiteCfsSubhourly.append(siteCfsSubhourly[0])
            allSiteCfsSubhourly.append([siteId] + siteCfsSubhourly[1][1:])
    else:
        allSiteCfsHourly.append([siteId] + siteCfsHourly[1][1:])  # replace header w/ site ID
        if subHour:
            allSiteCfsSubhourly.append([siteId] + siteCfsSubhourly[1][1:])

    assert (len(siteCfsHourly[0]) == len(allSiteCfsHourly[0]))

    if subHour:
        assert (len(siteCfsSubhourly[1]) == len(allSiteCfsSubhourly[1]))


def getBestWindIdsInZone(windDir, windCapacInZone, startWindCapacForCFs, fipsToZones, fipsToPolys, currZone):
    # ewdMetadataFilename = os.path.join(windDir,'eastern_wind_dataset_site_summary.csv')

    metadata = readCSVto2dList(os.path.join(windDir, 'toolkit_sites_v7_SERC.csv'))
    ewdCfCol = metadata[0].index('capacity_factor')
    ewdCapacCol = metadata[0].index('capacity')
    ewdSiteNumberCol = metadata[0].index('site_id')

    idAndCapacs = getWindOrSolarIdsInZonesDecreasingCF(metadata, windCapacInZone, ewdCfCol, ewdCapacCol,
                                                       ewdSiteNumberCol, startWindCapacForCFs, fipsToZones,
                                                       fipsToPolys, currZone, 'Wind')
    return idAndCapacs, metadata


def getWindOrSolarIdsInZonesDecreasingCF(metadata, capacInZone, cfCol, capacCol, siteNumberOrFileCol,
                                         startRECapacForCFs, fipsToZones, fipsToPolys, currZone, windOrSolar):
    """

    Assumes 70 and 12 MW capacities for wind and solar, respectively, rather than using capacities of plants
    attached to them wind database. Do this because wind capacity in WIND database is unrealistically small.

    :param metadata:
    :param capacInZone:
    :param cfCol:
    :param capacCol:
    :param siteNumberOrFileCol:
    :param startRECapacForCFs:
    :param fipsToZones:
    :param fipsToPolys:
    :param currZone:
    :param windOrSolar:
    :return:
    """
    idAndCapacs = [['Id', 'DatasetCapacity', 'FleetCapacity']]
    (cfs, capacs, siteNumbers) = getPlantInfoInZone(metadata, startRECapacForCFs, cfCol, capacCol, siteNumberOrFileCol,
                                                    fipsToZones, fipsToPolys, currZone)
    currZoneCapac = 0
    while currZoneCapac < capacInZone:
        if len(cfs) == 0:
            fleetCapac = capacInZone - currZoneCapac
            siteName = 'NoMoreSites'
            idAndCapacs.append([siteName, fleetCapac, fleetCapac])
            currZoneCapac += fleetCapac
        else:

            # get index of site with max CF
            maxCfIdx = cfs.index(max(cfs))

            if cfs[maxCfIdx] > 1E-3:  # some solar sites have 0 CF!
                if windOrSolar == 'Wind':
                    datasetCapac = 70  # avg wind fleet capac in US (NEEDS v5.15, 860 2015)
                else:
                    datasetCapac = 12  # avg solar capac in US (NEEDS v5.15, 860 2015)

                # Trim capacity if unit capacity > spare capacity before reach state capac - capac when start saving CFs
                fleetCapac = min(datasetCapac, capacInZone - currZoneCapac)
                if currZoneCapac < startRECapacForCFs <= currZoneCapac + fleetCapac:
                    fleetCapac = startRECapacForCFs - currZoneCapac
                currZoneCapac += fleetCapac

                if currZoneCapac > startRECapacForCFs:
                    # only start appending to list after it is greater than existing capacity
                    idAndCapacs.append([siteNumbers[maxCfIdx], datasetCapac, fleetCapac])

            # remove from lists
            cfs.pop(maxCfIdx)
            capacs.pop(maxCfIdx)
            siteNumbers.pop(maxCfIdx)

    return idAndCapacs


def getPlantInfoInZone(metadata, startRECapacForCFs, cfCol, capacCol, siteNumberOrFileCol, fipsToZones, fipsToPolys,
                       currZone, *args):
    """

    Match by zone

    :param metadata:
    :param startRECapacForCFs:
    :param cfCol:
    :param capacCol:
    :param siteNumberOrFileCol:
    :param fipsToZones:
    :param fipsToPolys:
    :param currZone:
    :param args:
    :return:
    """

    if 'longitude' in metadata[0]:
        # wind has lat/long already given

        latCol, lonCol = metadata[0].index('latitude'), metadata[0].index('longitude')
        plantsInRegionOrZone = [row for row in metadata[1:] if locInZone(float(row[latCol]), float(row[lonCol]),
                                                                         currZone, fipsToZones, fipsToPolys)]
    else:
        # solar lat/lon is in filename

        nameCol = metadata[0].index('File')
        plantsInRegionOrZone = [row for row in metadata[1:] if locInZone(float(getCoordsFromFilename(row[nameCol])[0]),
                                                                         float(getCoordsFromFilename(row[nameCol])[1]),
                                                                         currZone, fipsToZones, fipsToPolys)]
    cfs = [float(row[cfCol]) for row in plantsInRegionOrZone]
    capacs = [float(row[capacCol]) for row in plantsInRegionOrZone]
    siteNumbers = [row[siteNumberOrFileCol] for row in plantsInRegionOrZone]

    return cfs, capacs, siteNumbers


def calcCapacWtdFleetCf(idAndCapac, siteCfsHr, siteCfsSubhr):
    fleetCapacCol = idAndCapac[0].index('FleetCapacity')
    idToFleetCapac = getFleetToCapacDict(idAndCapac)
    capacWtdCfsHr = calcCapacWtdFleetCfHrOrSubhr(idToFleetCapac, siteCfsHr)
    capacWtdCfsSubhr = calcCapacWtdFleetCfHrOrSubhr(idToFleetCapac, siteCfsSubhr)
    return capacWtdCfsHr, capacWtdCfsSubhr


def getFleetToCapacDict(idAndCapac):
    idCol = idAndCapac[0].index('Id')
    fleetCapacCol = idAndCapac[0].index('FleetCapacity')
    idToFleetCapac = dict()
    for row in idAndCapac[1:]:
        idToFleetCapac[row[idCol]] = row[fleetCapacCol]
    return idToFleetCapac


def calcCapacWtdFleetCfHrOrSubhr(idToFleetCapac, siteCfs):
    (totalCapac, totalGen) = (0, [])
    if siteCfs is not None:
        for row in siteCfs[1:]:
            (currId, currCfs) = (row[0], row[1:])
            currCapac = idToFleetCapac[currId]
            totalCapac += currCapac
            gens = [val * currCapac for val in currCfs]
            if totalGen == []:
                totalGen = copy.copy(gens)
            else:
                totalGen = [totalGen[idx] + gens[idx] for idx in range(len(gens))]
        capacWtdCfs = [copy.deepcopy(siteCfs[0])]
        capacWtdCfs.append(['AnnualAvgCf'] + [val / totalCapac for val in totalGen])
    else:
        capacWtdCfs = None

    return capacWtdCfs


def getWindSiteCfs(windDir, siteId, siteCapac, desiredTz, windGenDataYr, subHour=False):
    """

    :param windDir: dir w/ wind data
    :param siteId: site ID to get gen data for
    :param siteCapac: wind site capac
    :param desiredTz: desired timezone
    :param windGenDataYr: year for wind gen data
    :param subHour: if True reads sub hourly data
    :return: 2 2d lists, both have first row = datetime. 1 2d list = hourly CFs, 2nd 2d list = subhourly CFs.
             Also row labels
    """
    hourlyFile = 'powerhourly_' + siteId + '.csv'
    hourlyGen = readCSVto2dList(os.path.join(windDir, 'hourlyPowerSERC', hourlyFile))
    datetimeAndHourlyGen = convertTimeToDatetimeInTgtTz(hourlyGen, 'wind', siteId, desiredTz, 'UTC')
    datetimeAndHourlyGenInYr = [datetimeAndHourlyGen[0]] + [row for row in datetimeAndHourlyGen[1:]
                                                            if row[0].year == windGenDataYr]
    hourlyCfs = convertToCfs(datetimeAndHourlyGenInYr, siteCapac)

    if subHour:
        subhourlyFile = 'powersubhourly_' + siteId + '.csv'
        subhourlyGen = readCSVto2dList(os.path.join(windDir, 'subhourlyPowerSERC', subhourlyFile))
        datetimeAndSubhourlyGen = convertTimeToDatetimeInTgtTz(subhourlyGen, 'wind', siteId, desiredTz, 'UTC')
        datetimeAndSubhourlyGenInYr = [datetimeAndSubhourlyGen[0]] + [row for row in datetimeAndSubhourlyGen[1:]
                                                                      if row[0].year == windGenDataYr]
        subhourlyCfs = convertToCfs(datetimeAndSubhourlyGenInYr, siteCapac)
    else:
        subhourlyCfs = None

    return hourlyCfs, subhourlyCfs


def convertTimeToDatetimeInTgtTz(genData, windOrSolar, siteOrFilename, tgtTz, siteTz):
    """Converts datetimes to target timezone

    :param genData: gen data (2d list w/ datetime in col 1 and gen data in col 2)
    :param windOrSolar: (string) whether processing wind or solar gen data
    :param siteOrFilename: site ID or filename
    :param tgtTz:
    :param siteTz:
    :return: 2d list w/ gen data (datetime in col 1, gen data in col 2)
    """
    datetimeAndGen = [['datetime' + tgtTz, 'power(MW)' + siteOrFilename]]
    tzOffsetDict = {'UTCtoCST': -6, 'CSTtoCST': 0, 'ESTtoCST': -1, 'CSTtoEST': 1, 'UTCtoEST': -5,
                    'ESTtoEST': 0}
    timezoneOffset = tzOffsetDict[siteTz + 'to' + tgtTz]

    if windOrSolar == 'wind':
        dateAndTimeCol = genData[0].index('Datetime')
        genCol = genData[0].index(siteOrFilename)
    elif windOrSolar == 'solar':
        dateAndTimeCol = genData[0].index('LocalTime')
        genCol = genData[0].index('Power(MW)')
    for row in genData[1:]:
        if windOrSolar == 'wind':
            year, month, day, hour, minute = divideWindDt(row[dateAndTimeCol])
        else:
            year, month, day, hour, minute = divideSolarDatetime(row[dateAndTimeCol])
        rowDatetime = datetime.datetime(year, month, day, hour, minute)
        rowDatetimeCST = rowDatetime + datetime.timedelta(hours=timezoneOffset)
        datetimeAndGen.append([rowDatetimeCST, float(row[genCol])])
    return datetimeAndGen


# Return year,month,day,hour,minute from datetime in wind gen data
# Given as YYYY-MM-DD HH:MM:SS
def divideWindDt(windDt):
    windDate, windTime = windDt.split(' ')
    dateSplit = windDate.split('-')
    timeSplit = windTime.split(':')
    return (int(dateSplit[0]), int(dateSplit[1]), int(dateSplit[2]),
            int(timeSplit[0]), int(timeSplit[1]))


# Return year,month,day,hour,minute from datetime in solar gen data
def divideSolarDatetime(solarDatetime):
    baseYear = 2000  # solar year is given as '06', so add 2000
    solarDate, solarTime = solarDatetime.split(' ')
    solarDateSplit = solarDate.split('/')
    solarTimeSplit = solarTime.split(':')
    return (baseYear + int(solarDateSplit[2]), int(solarDateSplit[0]), int(solarDateSplit[1]),
            int(solarTimeSplit[0]), int(solarTimeSplit[1]))


# Converts subhourly to hourly gen by counting # of power output entires for each
# hour and generator, then averaging them together.
# Inputs: subhourly (10 or 5 min) power output (2d list, col 1 = datetime CST,
# col 2 = gen)
# Outputs: average hourly gen *2d list, col 1 = datetime CST,
# col 2 = average hourly gen for each gen).
def convertGenToHourly(genSubhourly, tgtTz):
    datetimeCol = genSubhourly[0].index('datetime' + tgtTz)
    hourlyGen = [copy.deepcopy(genSubhourly[0])]
    countGen = [copy.deepcopy(genSubhourly[0])]
    hourlyGenAverage = [copy.deepcopy(genSubhourly[0])]
    dtHourToRowDict = dict()
    lastRowDtToHour = datetime.datetime(1980, 1, 1, 1, 1)  # random datetime
    for row in genSubhourly[1:]:
        rowDt = row[datetimeCol]
        rowDtToHour = datetime.datetime(rowDt.year, rowDt.month, rowDt.day, rowDt.hour, 0)
        if rowDtToHour == lastRowDtToHour:
            hourlyGen[-1][1] += row[1]
            countGen[-1][1] += 1
        else:
            hourlyGen.append([rowDtToHour] + [row[1]])
            countGen.append([rowDtToHour] + [1])
        lastRowDtToHour = rowDtToHour
    for idx in range(1, len(hourlyGen)):
        hourlyGenAverage.append([hourlyGen[idx][datetimeCol], hourlyGen[idx][1] / countGen[idx][1]])
    return hourlyGenAverage


# Inputs: 2d list (datetime 1st col, gen 2nd col, w/ headers), capacity of curr wind gen
# Outputs: 2d list (datetime first row, gen 2nd row, w/ labels)
def convertToCfs(datetimeAndGen, siteCapac):
    dateCol, genCol = 0, 1
    dateInfoHoriz = [row[dateCol] for row in datetimeAndGen]
    cfsHoriz = [datetimeAndGen[0][genCol]] + [float(row[genCol]) / siteCapac for row in datetimeAndGen[1:]]
    return [dateInfoHoriz, cfsHoriz]


##### SOLAR CFS #####
# All CFs output in CST
def getSolarCFs(solarUnits, solarDir, startSolarCapacForCFs, desiredTz, fipsToZones, fipsToPolys, currZone):
    # Get total wind capacity by state
    capacCol = solarUnits[0].index('Capacity (MW)')
    solarCapacInZone = sum([float(row[capacCol]) for row in solarUnits[1:]])

    # Get solar plants in NREL dataset per state until capacity met
    (solarFilenameAndCapac, solarFilenameAndCapacAndTz,
     solarMetadata) = getBestSolarIdsInStates(solarDir, solarCapacInZone, startSolarCapacForCFs,
                                              fipsToZones, fipsToPolys, currZone)

    # Import CFs for each wind plant
    idCol = solarFilenameAndCapacAndTz[0].index('Id')
    datasetCapacCol = solarFilenameAndCapacAndTz[0].index('DatasetCapacity')
    tzCol = solarFilenameAndCapacAndTz[0].index('Timezone')

    allSiteCfsHourly, allSiteCfsSubhourly, avgFleetCfHr = [], [], []

    for site in solarFilenameAndCapacAndTz[1:]:
        (siteFilename, datasetSiteCapac, siteTz) = (site[idCol], site[datasetCapacCol], site[tzCol])
        if 'NoMoreSites' in siteFilename:
            print('No more solar sites!')
            if avgFleetCfHr == []:
                # If doing new RE CFs & existing RE > potential RE capac for CFs, can end up w/ nothing in allSiteCFs.
                # In that case, new RE CFs should jsut be fleet average CF.
                if allSiteCfsHourly == []:
                    a, tempAllSiteCfsHr, tempAllSiteCfsSubhr, tempFileAndCapac = getSolarCFs(solarUnits, solarDir, 0,
                                                                                             desiredTz, fipsToZones,
                                                                                             fipsToPolys, currZone)
                    avgFleetCfHr, avgFleetCfSubhr = calcCapacWtdFleetCf(tempFileAndCapac, tempAllSiteCfsHr,
                                                                        tempAllSiteCfsSubhr)
                else:
                    avgFleetCfHr, avgFleetCfSubhr = calcCapacWtdFleetCf(solarFilenameAndCapac, allSiteCfsHourly,
                                                                        allSiteCfsSubhourly)

            siteCfsHourly, siteCfsSubhourly = copy.deepcopy(avgFleetCfHr), copy.deepcopy(avgFleetCfSubhr)
        else:
            siteCfsHourly, siteCfsSubhourly = getSolarSiteCfs(solarDir, siteFilename, datasetSiteCapac, siteTz,
                                                              desiredTz)
        addSiteCfsToAggList(siteCfsHourly, siteCfsSubhourly, allSiteCfsHourly, allSiteCfsSubhourly, siteFilename,
                            subHour=True)

    allSiteCfsHourOfYear = [['HourOfYear'] + [val for val in range(1, 8761)]] + copy.deepcopy(allSiteCfsHourly[1:])

    return (allSiteCfsHourOfYear, allSiteCfsHourly, allSiteCfsSubhourly, solarFilenameAndCapac)


def getBestSolarIdsInStates(solarDir, solarCapacByState, startSolarCapacForCFs,
                            fipsToZones, fipsToPolys, currZone):
    solarMetadataFilename = os.path.join(solarDir, 'SolarCapacityFactorsNRELSERC.csv')
    solarMetadata = readCSVto2dList(solarMetadataFilename)
    solarCfCol = solarMetadata[0].index('CF')
    solarCapacCol = solarMetadata[0].index('PlantSize')
    solarFilenameCol = solarMetadata[0].index('File')
    idAndCapacs = getWindOrSolarIdsInZonesDecreasingCF(solarMetadata, solarCapacByState,
                                                       solarCfCol, solarCapacCol, solarFilenameCol,
                                                       startSolarCapacForCFs, fipsToZones,
                                                       fipsToPolys, currZone, 'Solar')
    idCol = idAndCapacs[0].index('Id')
    solarIdAndCapacAndTz = [idAndCapacs[0] + ['Timezone']]
    for idAndCapac in idAndCapacs[1:]:
        solarIdAndCapacAndTz.append(idAndCapac + [timezoneOfSolarSite(idAndCapac[idCol], currZone)])

    return (idAndCapacs, solarIdAndCapacAndTz, solarMetadata)


def timezoneOfSolarSite(solarFilename, currZone):
    # in KY or TN, which is half & half CST & EST
    kyLine = [(36.601261, -84.861318), (38.048865, -86.251172)]
    tnLine = [(36.601261, -85.076648), (34.997791, -85.605184)]
    gaAlBorder = [(34.9894, -85.644), (29.893, -84.854)]
    if currZone == 'S_C_TVA' or currZone == 'S_C_KY' or currZone == 'S_SOU':
        if currZone == 'S_C_TVA':
            line = tnLine
        elif currZone == 'S_C_KY':
            line = kyLine
        elif currZone == 'S_SOU':
            line = gaAlBorder

        if solarFilename == 'NoMoreSites':
            print('No More solar sites at zone {}'.format(currZone))
            (siteLat, siteLong) = get_centroid_zone(genparam, currZone)
        else:
            (siteLat, siteLong) = getCoordsFromFilename(solarFilename)

        if siteEastOfLine(line, float(siteLat), float(siteLong)):
            tz = 'EST'
        else:
            tz = 'CST'

    elif currZone == 'S_VACA':
        tz = 'EST'

    return tz


def getCoordsFromFilename(solarFilename):
    latStart = solarFilename.index('_') + 1
    latEnd = solarFilename[latStart:].index('_')
    lat = solarFilename[latStart:(latStart + latEnd)]
    longStart = solarFilename.index('-')
    longEnd = solarFilename[longStart:].index('_')
    longitude = solarFilename[longStart:(longStart + longEnd)]
    return (lat, longitude)


# Long = x coord, lat = y coord
def siteEastOfLine(line, siteLat, siteLong):
    (deltaLat, deltaLong) = (line[0][0] - line[1][0], line[0][1] - line[1][1])
    lineSlope = deltaLat / deltaLong
    intercept = line[1][0] - lineSlope * line[1][1]  # b = y - mx
    longOnLineForSiteLat = (siteLat - intercept) / lineSlope  # x = (y-b)/m
    return siteLong > longOnLineForSiteLat  # long decreases (more negative) west across US


# Solar data can be CST or EST, but move all to desiredTz. Output from this function
# is generation in desiredTz for entire year.
def getSolarSiteCfs(solarDir, siteFilename, siteCapac, siteTz, desiredTz):
    genData = readCSVto2dList(os.path.join(solarDir, 'AllSERC', siteFilename))
    datetimeAndGenSubhourly = convertTimeToDatetimeInTgtTz(genData, 'solar', siteFilename, desiredTz, siteTz)
    datetimeAndGenHourly = convertGenToHourly(datetimeAndGenSubhourly, desiredTz)
    datetimeAndGenSuhHourlyInYr = slimSolarGenToYear(datetimeAndGenSubhourly, desiredTz, siteFilename)
    datetimeAndGenHourlyInYr = slimSolarGenToYear(datetimeAndGenHourly, desiredTz, siteFilename)
    subhourlyCfs = convertToCfs(datetimeAndGenSuhHourlyInYr, siteCapac)
    hourlyCfs = convertToCfs(datetimeAndGenHourlyInYr, siteCapac)
    return hourlyCfs, subhourlyCfs


# Slims solar generation data to full year.
def slimSolarGenToYear(gen, desiredTz, siteFilename):
    dtCol, genCol = gen[0].index('datetime' + desiredTz), gen[0].index('power(MW)' + siteFilename)
    currYear = mode([row[dtCol].year for row in gen[1:]])
    firstYrDt = gen[[row[dtCol].year for row in gen[1:]].index(currYear) + 1][dtCol]  # add 1 since skip 1st row
    timeStepMins = (gen[2][dtCol] - gen[1][dtCol]).seconds // 60  # convert seconds to minutes
    missingVals = firstYrDt.hour * 60 // timeStepMins
    if missingVals > 0:
        genInYr = [gen[0]] + [[datetime.datetime(currYear, 1, 1, 0, 0) + datetime.timedelta(minutes=timeStepMins * ctr),
                               gen[1][genCol]] for ctr in range(missingVals)] + [row for row in gen[1:] if
                                                                                 row[dtCol].year == currYear]
    else:
        genInYr = copy.deepcopy(gen)
    return genInYr


################################################################################
################################################################################
################################################################################


################################################################################
####### MISC. DATA #############################################################
################################################################################
# STATETIMEZONES = {'North Carolina':'EST','South Carolina':'EST','Virginia':'EST',
#                   'Georgia':'EST','Mississippi':'CST','Alabama':'CST','Louisiana':'CST',
#                   'Missouri':'CST','Arkansas':'CST','Illinois':'CST','Kentucky':'CSTorEST',
#                   'Tennessee':'CSTorEST','Texas':'CST'}

# STATEABBREVS = {'North Carolina':'NC','South Carolina':'SC','Virginia':'VA',
#                   'Georgia':'GA','Mississippi':'MS','Alabama':'AL','Louisiana':'LA',
#                   'Missouri':'MO','Arkansas':'AR','Illinois':'IL','Kentucky':'KY',
#                   'Tennessee':'TN','Texas':'TX'}
################################################################################
################################################################################
################################################################################

################################################################################
####### TEST FUNCTIONS #########################################################
################################################################################
def testTimezoneAssignment():
    print('Testing timezone assignment')
    assert (timezoneOfSolarSite('Arkansas', 'hello.csv') == 'CST')
    assert (timezoneOfSolarSite('South Carolina', 'hello.csv') == 'EST')
    assert (getCoordsFromFilename('Actual_29.15_-90.15_2006_UPV_140MW_5_Min.csv') == ('29.15', '-90.15'))
    kyLine = [(36.601261, -84.861318), (38.048865, -86.251172)]
    tnLine = [(36.601261, -85.076648), (34.997791, -85.605184)]
    assert (siteEastOfLine(kyLine, 35, -500) == False)
    assert (siteEastOfLine(kyLine, 30, -1) == True)
    assert (siteEastOfLine(kyLine, 36.55, -88.15) == False)
    assert (siteEastOfLine(kyLine, 37.493, -82.679) == True)
    assert (siteEastOfLine(kyLine, 37.525, -84.557) == True)
    assert (siteEastOfLine(tnLine, 35, -500) == False)
    assert (siteEastOfLine(tnLine, 30, -1) == True)
    assert (siteEastOfLine(tnLine, 35.066, -85.009) == True)
    assert (siteEastOfLine(tnLine, 36.541, -86.101) == False)
    assert (timezoneOfSolarSite('Kentucky', 'Actual_38.05_-84.55_2006_DPV_34MW_5_Min.csv') == 'EST')
    assert (timezoneOfSolarSite('Kentucky', 'Actual_38.85_-84.75_2006_DPV_31MW_5_Min.csv') == 'EST')
    assert (timezoneOfSolarSite('Kentucky', 'Actual_36.55_-88.15_2006_UPV_29MW_5_Min.csv') == 'CST')
    assert (timezoneOfSolarSite('Tennessee', 'Actual_36.65_-87.35_2006_DPV_35MW_5_Min.csv') == 'CST')
    assert (timezoneOfSolarSite('Tennessee', 'Actual_34.95_-85.25_2006_DPV_38MW_5_Min.csv') == 'EST')
    print('Passed')

################################################################################
################################################################################
################################################################################
