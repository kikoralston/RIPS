# Michael Craig
# August 11, 2016
# Script for modifying generator capacity with water T from RBM.

# 1) Gets list of all grid cells from .spat.
# 2) Saves T data for each segment in each grid cell into folders for each cell.
# 3) Averages T data across all segments for each cell and saves that data in folder for each cell.
# 4) Creates dictionaries mapping generators to cells and vice versa.
# 5) Using dictionary, pairs generator with water T time series and modifies capacity.

import os, csv, copy
import numpy as np
import pandas as pd
import datetime as dt
from AuxFuncs import *
from GAMSAuxFuncs import *
from CurtailmentRegressions import (calcCurtailmentForGenOrTech, loadRegCoeffs, getKeyCurtailParams,
                                    getCoeffsForGenOrTech)


################################################################################
####### MASTER FUNCTIONS #######################################################
################################################################################
def processRBMDataIntoIndividualCellFiles(rbmDataDir, tempAndSpatFilename,
                                          rbmOutputDir, nsegFilename, locPrecision, outputHeaders, numCellsToProcess):
    gridCellLatLongs = getAllGridCellLatLongsInSpatFile(rbmDataDir, tempAndSpatFilename)
    saveAndAverageTDataForAllCells(gridCellLatLongs, rbmDataDir, rbmOutputDir, tempAndSpatFilename,
                                   nsegFilename, locPrecision, outputHeaders, numCellsToProcess)


# Outputs: dictionaryof ORIS+UNITID to 2d vertical list of [datetime,curtailment(mw)],
# list of (ORIS+UNIT,gen lat long, cell lat long), and list of hourly curtailments for existing gens.
def determineHrlyCurtailmentsForExistingGens(locPrecision, genFleet, rbmOutputDir, curtailmentYear,
                                             modelName, ptCurtailed, ptCurtailedRegs, incCurtailments, incRegs,
                                             envRegMaxT, dataRoot,
                                             coolDesignT, resultsDir):
    (genToCellLatLongsDict, cellLatLongToGenDict,
     genToCellLatLongsList) = getGenToCellAndCellToGenDictionaries(genFleet)
    hrlyCurtailmentsAllGensInTgtYr, hrlyCurtailmentsList = calculateGeneratorCurtailments(cellLatLongToGenDict,
                                                                                          rbmOutputDir, locPrecision,
                                                                                          curtailmentYear, genFleet,
                                                                                          modelName, ptCurtailed,
                                                                                          ptCurtailedRegs,
                                                                                          incCurtailments, incRegs,
                                                                                          envRegMaxT, dataRoot,
                                                                                          coolDesignT, resultsDir)
    return hrlyCurtailmentsAllGensInTgtYr, genToCellLatLongsList, hrlyCurtailmentsList


################################################################################
################################################################################
################################################################################

################################################################################
####### PROCESS RBM DATA INTO INDIVIDUAL CELL FILES ############################
################################################################################
# Get all grid cell lat/longs in spat file
def getAllGridCellLatLongsInSpatFile(rbmDataDir, tempAndSpatFilename):
    spatFile = os.path.join(rbmDataDir, tempAndSpatFilename + '.Spat')
    f = open(spatFile, 'r')
    lineCtr = 0
    gridCellLatLongs = set()
    while 1:
        line = f.readline().rstrip("\n")
        lineCtr += 1  # current line number
        if line == "": break
        lineSplit = line.split()
        (cellLat, cellLong) = (float(lineSplit[4]), float(lineSplit[5]))
        gridCellLatLongs.add((cellLat, cellLong))
    f.close()
    return gridCellLatLongs


def saveAndAverageTDataForAllCells(gridCellLatLongs, rbmDataDir, rbmOutputDir, tempAndSpatFilename,
                                   nsegFilename, locPrecision, outputHeaders, numCellsToProcess):
    """Average temperature data for all segments in each cell, then save all data

    :param gridCellLatLongs:
    :param rbmDataDir:
    :param rbmOutputDir:
    :param tempAndSpatFilename:
    :param nsegFilename:
    :param locPrecision:
    :param outputHeaders:
    :param numCellsToProcess:
    """
    print('Num grid cells:', len(gridCellLatLongs))
    (cellCtr, totalCellCtr, tgtCellLatLongs) = (0, 0, [])
    for (cellLat, cellLong) in gridCellLatLongs:
        cellCtr += 1
        totalCellCtr += 1
        tgtCellLatLongs.append((cellLat, cellLong))
        if cellCtr == numCellsToProcess or totalCellCtr == len(gridCellLatLongs):
            saveAndAverageTDataForGridCells(rbmDataDir, rbmOutputDir, tempAndSpatFilename,
                                            nsegFilename, locPrecision, tgtCellLatLongs, outputHeaders)
            (cellCtr, tgtCellLatLongs) = (0, [])
            print('Cell count:', totalCellCtr)


def saveAndAverageTDataForGridCells(rbmDataDir, rbmOutputDir, tempAndSpatFilename,
                                    nsegFilename, locPrecision, tgtCellLatLongs, outputHeaders):
    outputDirs = createOutputDirs(rbmOutputDir, tempAndSpatFilename, locPrecision, tgtCellLatLongs)
    (numTotalSegments, numDays) = getNumSegmentsAndDays(rbmDataDir, nsegFilename)
    (cellInfo, mapSegmentNumToLatLong, mapLatLongToSegmentNum) = getLinesInSpatFileForCell(rbmDataDir,
                                                                                           tempAndSpatFilename,
                                                                                           tgtCellLatLongs,
                                                                                           numTotalSegments)
    (allSegmentData, waterTAllSegments) = readCellTemperatureData(rbmDataDir, tempAndSpatFilename,
                                                                  cellInfo, numTotalSegments, outputHeaders)
    writeCellTemperatureData(allSegmentData, cellInfo, outputDirs, locPrecision,
                             numDays, outputHeaders, mapSegmentNumToLatLong)
    averageAndSaveWaterTOverAllSegmentsInCells(waterTAllSegments, outputDirs, locPrecision,
                                               mapLatLongToSegmentNum, mapSegmentNumToLatLong)


def createOutputDirs(baseOutputDir, tempAndSpatFilename, locPrecision, tgtCellLatLongs):
    """CREATE OUTPUT DIRECTORY FOR CELL T DATA

    :param baseOutputDir:
    :param tempAndSpatFilename:
    :param locPrecision:
    :param tgtCellLatLongs:
    :return:
    """
    outputDirs = dict()
    for (cellLat, cellLong) in tgtCellLatLongs:
        outputDir = os.path.join(baseOutputDir,
                                 '%.*f_%.*f' % (locPrecision, cellLat, locPrecision, cellLong))
        if not os.path.exists(outputDir): os.makedirs(outputDir)
        outputDirs[(cellLat, cellLong)] = outputDir
    return outputDirs


def getNumSegmentsAndDays(dataDir, nsegFilename):
    """EXTRACT TOTAL NUM SEGMENTS AND DAYS IN MODEL RUN

    :param dataDir:
    :param nsegFilename:
    :return:
    """
    nseg_nday = np.loadtxt(os.path.join(dataDir, nsegFilename))
    (numTotalSegments, numDays) = (int(nseg_nday[0]), int(nseg_nday[1]))
    return (numTotalSegments, numDays)


def getLinesInSpatFileForCell(dataDir, tempAndSpatFilename, tgtCellLatLongs, numTotalSegments):
    """GET ALL LINES IN SPATIAL FILE THAT CORRESPOND TO CELL

    Gets line numbers in spatial file that match cell. Includes all reaches & segments
    EXCEPT segment 2 in outlet cell (= max segment #) b/c no data for that segment.

    :param dataDir:
    :param tempAndSpatFilename:
    :param tgtCellLatLongs:
    :param numTotalSegments:
    :return:
    """
    (mapSegmentNumToLatLong, mapLatLongToSegmentNum) = (dict(), dict())
    spatFile = os.path.join(dataDir, tempAndSpatFilename + '.Spat')  # .Temp file
    f = open(spatFile, 'r')
    totalSegmentNum = 0
    cellInfo = dict()  # line num of target grid cell : [reach index, segment index] for that line
    while 1:  # loop over each line in the .Spat file
        line = f.readline().rstrip("\n")
        totalSegmentNum = totalSegmentNum + 1  # current line number
        if line == "": break
        lineSplit = line.split()
        (lat, lon) = (float(lineSplit[4]), float(lineSplit[5]))
        if (lat, lon) in tgtCellLatLongs:  # line is for one of target grid cells
            if totalSegmentNum != numTotalSegments:  # if not segment 2 in outlet cell (no data)
                (reachIndex, segmentIndex) = (int(lineSplit[0]), int(lineSplit[6]))  # seg = 1 or 2; reach varies
                cellInfo[totalSegmentNum] = [reachIndex, segmentIndex]
                if totalSegmentNum not in mapSegmentNumToLatLong:
                    mapSegmentNumToLatLong[totalSegmentNum] = (lat, lon)
                    if (lat, lon) in mapLatLongToSegmentNum:
                        mapLatLongToSegmentNum[(lat, lon)].append(totalSegmentNum)
                    else:
                        mapLatLongToSegmentNum[(lat, lon)] = [totalSegmentNum]
    f.close()
    return (cellInfo, mapSegmentNumToLatLong, mapLatLongToSegmentNum)


def readCellTemperatureData(dataDir, tempAndSpatFilename, cellInfo, numTotalSegments, outputHeaders):
    """READ CELL TEMPERATURE DATA

    Temperature data: hour 1, all cells + segments + reaches, then hour 2, all cells + segments + reaches, etc.

    :param dataDir:
    :param tempAndSpatFilename:
    :param cellInfo:
    :param numTotalSegments:
    :param outputHeaders:
    :return:
    """
    tempFile = os.path.join(dataDir, tempAndSpatFilename + '.Temp')  # .Temp file

    # total segment number (=line #in .spat file): np array of year, month, day, flow(cfs), streamTemp(degC).
    allSegmentData = {}

    for totalSegmentNum in cellInfo:
        allSegmentData[totalSegmentNum] = []  # initialize b/c adding to lists in dict later

    waterTAllSegments = {}  # store just water T for each segment
    for totalSegmentNum in cellInfo: waterTAllSegments[totalSegmentNum] = dict()
    f = open(tempFile, 'r')
    lineCtr = 0
    while 1:  # loop over each line in the .Temp file
        tempLineRaw = f.readline().rstrip("\n")
        lineCtr = lineCtr + 1
        if tempLineRaw == "": break
        totalSegmentNum = lineCtr % numTotalSegments  # corresponding line number in the .Spat file
        if totalSegmentNum in cellInfo:
            (tLineData, year, month, day, waterT) = processRawTemperatureLineData(tempLineRaw, outputHeaders)
            # Save entire line of segment data
            allSegmentData[totalSegmentNum].append(tLineData)
            # Save just water T data of segment into nested dictionary - total segment number: date : waterT
            waterTAllSegments[totalSegmentNum][createDateLabel(year, month, day)] = waterT
    f.close()
    for i in allSegmentData: allSegmentData[i] = np.asarray(allSegmentData[i])
    return (allSegmentData, waterTAllSegments)


def processRawTemperatureLineData(tempLineRaw, outputHeaders):
    (tDataDict, year, month, day, waterT) = readRawTemperatureLineData(tempLineRaw)
    tLineData = processTempDataDictIntoRow(tDataDict, outputHeaders)
    return (tLineData, year, month, day, waterT)


def readRawTemperatureLineData(tempLineRaw):
    lineSplit = tempLineRaw.split()
    decimal_year = lineSplit[0]
    year = int(decimal_year.split('.')[0])
    day_of_year = int(lineSplit[1])
    if day_of_year > 360 and float('0.' + decimal_year.split('.')[1]) <= 0.005:  # correct bad decimal year integer part
        year = year - 1
    date = dt.datetime(year, 1, 1) + dt.timedelta(days=day_of_year - 1)  # convert day of year to date
    flow = float(lineSplit[8])
    streamT = float(lineSplit[5])
    headwaterT = float(lineSplit[6])
    airT = float(lineSplit[7])
    tDataDict = {'year': year, 'month': date.month, 'day': date.day, 'flow': flow, 'streamT': streamT,
                 'headwaterT': headwaterT, 'airT': airT}
    return (tDataDict, year, date.month, date.day, streamT)


def processTempDataDictIntoRow(tempDataDict, outputHeaders):
    tLineData = [''] * len(outputHeaders)
    (yearCol, monthCol, dayCol, flowCol, streamTCol, headTCol, airTCol) = getDataColNums(outputHeaders)
    tLineData[yearCol] = tempDataDict['year']
    tLineData[monthCol] = tempDataDict['month']
    tLineData[dayCol] = tempDataDict['day']
    tLineData[flowCol] = tempDataDict['flow']
    tLineData[streamTCol] = tempDataDict['streamT']
    tLineData[headTCol] = tempDataDict['headwaterT']
    tLineData[airTCol] = tempDataDict['airT']
    return tLineData


def getDataColNums(cols):
    (yearCol, monthCol, dayCol, flowCol, streamTCol, headTCol, airTCol) = (cols.index('Year'),
                                                                           cols.index('Month'), cols.index('Day'),
                                                                           cols.index('Streamflow(cfs)'),
                                                                           cols.index('StreamT(degC)'),
                                                                           cols.index('HeadwaterT(degC)'),
                                                                           cols.index('AirT(degC)'))
    return (yearCol, monthCol, dayCol, flowCol, streamTCol, headTCol, airTCol)


def createDateLabel(year, month, day):
    return '%s-%s-%s' % (year, month, day)


###### WRITE CELL TEMPERATURE DATA ##############################################
def writeCellTemperatureData(allSegmentData, cellInfo, outputDirs, locPrecision, numDays,
                             outputHeaders, mapSegmentNumToLatLong):
    for totalSegmentNum in cellInfo:
        (cellLat, cellLong) = mapSegmentNumToLatLong[totalSegmentNum]
        outputDir = outputDirs[(cellLat, cellLong)]
        baseFilename = createBaseFilenameToReadOrWrite(locPrecision, cellLat, cellLong)
        fullFilename = baseFilename + '_reach%d_seg%d' % (cellInfo[totalSegmentNum][0], cellInfo[totalSegmentNum][1])
        fullFilePath = os.path.join(outputDir, fullFilename)
        f = open(fullFilePath, 'w')
        f.write(createHeaderStr(outputHeaders))
        dataCurr = allSegmentData[totalSegmentNum]
        for i in range(numDays):
            f.write('%d %d %d %.1f %.2f %.2f %.2f\n' % (dataCurr[i, 0], dataCurr[i, 1], dataCurr[i, 2],
                                                        dataCurr[i, 3], dataCurr[i, 4], dataCurr[i, 5], dataCurr[i, 6]))
        f.close()


def createBaseFilenameToReadOrWrite(locPrecision, inputLat, inputLong):
    return '%.*f_%.*f' % (locPrecision, inputLat, locPrecision, inputLong)


def createHeaderStr(outputHeaders):
    headStr = ''
    for idx in range(len(outputHeaders)):
        headStr += outputHeaders[idx]
        if outputHeaders[idx] != outputHeaders[-1]: headStr += ' '  # not last header
    return headStr + '\n'


###### AVERAGE AND SAVE CELL TEMPERATURE DATA ##################################
def averageAndSaveWaterTOverAllSegmentsInCells(waterTAllSegments, outputDirs, locPrecision,
                                               mapLatLongToSegmentNum, mapSegmentNumToLatLong):
    for (cellLat, cellLon) in mapLatLongToSegmentNum:
        currCellSegmentNums = mapLatLongToSegmentNum[(cellLat, cellLon)]
        waterTAllSegmentsInCell = isolateWaterTDataForCell(waterTAllSegments, currCellSegmentNums)
        waterTAvgInCell = averageWaterTOverAllSegmentsInCell(waterTAllSegmentsInCell)
        sortedAvgWaterT2dList = convertAvgWaterTTo2dList(waterTAvgInCell)
        verticalSortedAvgWaterT2dList = flip2dList(sortedAvgWaterT2dList)
        outputDir = outputDirs[(cellLat, cellLon)]
        saveAverageWaterT(verticalSortedAvgWaterT2dList, outputDir, locPrecision, cellLat, cellLon)


def isolateWaterTDataForCell(waterTAllSegments, currCellSegmentNums):
    waterTAllSegmentsInCell = dict()
    for segmentNum in currCellSegmentNums:
        waterTAllSegmentsInCell[segmentNum] = waterTAllSegments[segmentNum]
    return waterTAllSegmentsInCell


def averageWaterTOverAllSegmentsInCell(waterTAllSegments):
    waterTSums = dict()
    for segment in waterTAllSegments:
        waterTSegment = waterTAllSegments[segment]
        for date in waterTSegment:
            waterTSums[date] = waterTSums.get(date, 0) + waterTSegment[date]
    waterTAvgInCell = dict()
    numSegments = len(waterTAllSegments)
    for date in waterTSums: waterTAvgInCell[date] = waterTSums[date] / numSegments
    return waterTAvgInCell


def convertAvgWaterTTo2dList(waterTAvgInCell):
    sortedAvgWaterT2dList = []
    sortedDatesList = createDatesList(waterTAvgInCell)
    sortedAvgWaterT2dList.append(['Datetime'] + sortedDatesList)
    waterTAvgList = ['AverageWaterT(degC)'] + convertWaterTDictToSortedList(sortedDatesList, waterTAvgInCell)
    sortedAvgWaterT2dList.append(waterTAvgList)
    return sortedAvgWaterT2dList


# Outputs sorted list of datetimes
def createDatesList(waterTAvgInCell):
    datesList = []
    for listDate in waterTAvgInCell:
        (year, month, day) = getElementsOfDate(listDate)
        dateAsDatetime = dt.date(year, month, day)
        datesList.append(dateAsDatetime)
    return sorted(datesList)


def getElementsOfDate(listDate):
    year = int(listDate[:4])
    restOfDate = listDate[5:]
    monthEndIdx = restOfDate.index('-')
    month = int(restOfDate[:monthEndIdx])
    day = int(restOfDate[monthEndIdx + 1:])
    return (year, month, day)


def convertWaterTDictToSortedList(sortedDatesList, waterTDict):
    waterTDates = [''] * len(sortedDatesList)
    for listDate in waterTDict:
        (year, month, day) = getElementsOfDate(listDate)
        dateAsDatetime = dt.date(year, month, day)
        waterTDates[sortedDatesList.index(dateAsDatetime)] = waterTDict[listDate]
    return waterTDates


def flip2dList(list2d):
    flippedList = []
    for colIdx in range(len(list2d[0])):
        newList = []
        for rowIdx in range(len(list2d)):
            newList.append(list2d[rowIdx][colIdx])
        flippedList.append(newList)
    return flippedList


def saveAverageWaterT(waterTAvgs2dList, outputDir, locPrecision, inputLat, inputLong):
    filenameToSave = createAverageTFilename(locPrecision, inputLat, inputLong)
    write2dListToCSV(waterTAvgs2dList, os.path.join(outputDir, filenameToSave))


def createAverageTFilename(locPrecision, inputLat, inputLong):
    baseFilename = createBaseFilenameToReadOrWrite(locPrecision, inputLat, inputLong)
    return baseFilename + 'Average.csv'


################################################################################
################################################################################
################################################################################

def getGenToCellAndCellToGenDictionaries(genFleet):
    """MAP GENERATORS TO AND FROM RBM GRID CELLS

    Returns 2 dicts for all gens in fleet that have lat/long coords. 1)gen:[gen,(gen lat, gen long),
    (cell lat, cell long)]. 2) (cell lat, clell long):genID. & list of [[genID,(genlat,genlong),(celllat,celllong)]]

    :param genFleet:
    :return:
    """
    (fleetLatCol, fleetLongCol) = (genFleet[0].index('Latitude'), genFleet[0].index('Longitude'))
    genToCellLatLongsList = [['ORIS+UnitID', 'GenLat,GenLong', 'CellLat,CellLong']]
    genToCellLatLongsDict = dict()
    cellLatLongToGenDict = dict()
    for row in genFleet[1:]:
        genID = createGenSymbol(row, genFleet[0])
        (genLat, genLong) = (row[fleetLatCol], row[fleetLongCol])
        cellLoc = 'NA'
        if genLat != 'NA' and genLat != '':
            cellLoc = find125GridMaurerLatLong(float(genLat), float(genLong))
            if cellLoc in cellLatLongToGenDict:
                cellLatLongToGenDict[cellLoc].append(genID)
            else:
                cellLatLongToGenDict[cellLoc] = [genID]
        genToCellLatLongsDict[genID] = [genID, (genLat, genLong), cellLoc]
        genToCellLatLongsList.append([genID, (genLat, genLong), cellLoc])
    return (genToCellLatLongsDict, cellLatLongToGenDict,
            genToCellLatLongsList)


# Get (lat,lon) of 1/8 grid cell that a (lat, lon) point falls in
def find125GridMaurerLatLong(lat, lon):
    lat_grid = np.around(8.0 * lat - 0.5) / 8.0 + 0.0625
    lon_grid = np.around(8.0 * lon - 0.5) / 8.0 + 0.0625
    return (lat_grid, lon_grid)


################################################################################
################################################################################
def calculateGeneratorCurtailments(cellLatLongToGenDict, rbmOutputDir, locPrecision,
                                   curtailmentYear, genFleet, modelName, ptCurtailed, ptCurtailedRegs, incCurtailments,
                                   incRegs, envRegMaxT, dataRoot, coolDesignT, resultsDir):
    """CURTAIL GENERATOR CAPACITY WITH WATER TEMPERATURES

    Calculates generator curtailments. If generator isn't curtailed (not right plant type, cell data not avaialble,
    etc.), ignored here. CalculateHourlyCapacs... script handles those cases (assume curtailment = 0).
    cellLatLongToGenDict = dict of (cell lat, cell long):gen ID for all gens w/ lat/lon coords in genFleet,
    including gens that should nto be curtailed.

    :param cellLatLongToGenDict:
    :param rbmOutputDir:
    :param locPrecision:
    :param curtailmentYear:
    :param genFleet:
    :param modelName:
    :param ptCurtailed:
    :param ptCurtailedRegs:
    :param incCurtailments:
    :param incRegs:
    :param envRegMaxT:
    :param dataRoot:
    :param coolDesignT:
    :param resultsDir:
    :return:
    """
    hrlyCurtailmentsAllGensInTgtYr, hrlyCurtailmentsList = dict(), list()
    regCoeffs = loadRegCoeffs(dataRoot, 'capacity.json')  # dict of plant type: cooling type: cool design T: param: coeffs
    allCellFolders = os.listdir(rbmOutputDir)

    for (cellLat, cellLong) in cellLatLongToGenDict:  # this maps gen lat/lon to gen IDs; cell data may not exist
        cellFoldername = createBaseFilenameToReadOrWrite(locPrecision, cellLat, cellLong)
        if cellFoldername in allCellFolders:
            metAndWaterData = loadWaterAndMetData(cellLat, cellLong, rbmOutputDir, locPrecision,
                                                  curtailmentYear, cellFoldername, resultsDir)
            gensInCell = cellLatLongToGenDict[(cellLat, cellLong)]  # list of ORIS-UNITID in cell
            for gen in gensInCell:
                (plantType, hr, fuelAndCoalType, coolType, fgdType, state, capac) = getKeyCurtailParams(gen, genFleet)
                coeffs = getCoeffsForGenOrTech(plantType, hr, fuelAndCoalType, coolType,
                                               fgdType, ptCurtailed, regCoeffs, coolDesignT)
                if (coeffs is not None) and (plantType in ptCurtailedRegs):  # skip gens that aren't curtailed
                    hrlyCurtailmentsGen = calcCurtailmentForGenOrTech(plantType, hr, fuelAndCoalType,
                                                                      coolType, fgdType, state, capac, ptCurtailed,
                                                                      ptCurtailedRegs, metAndWaterData,
                                                                      incCurtailments, incRegs, envRegMaxT, coolDesignT,
                                                                      coeffs)
                    hrlyCurtailmentsAllGensInTgtYr[gen] = hrlyCurtailmentsGen
                    hrlyCurtailmentsList.append([gen] + hrlyCurtailmentsGen)
        else:
            print('Cell not in folders!:', cellFoldername)
    return hrlyCurtailmentsAllGensInTgtYr, hrlyCurtailmentsList


def loadWaterAndMetData(cellLat, cellLong, rbmOutputDir, locPrecision, curtailmentYear,
                        cellFoldername, resultsDir):
    """LOAD WATER AND MET DATA BY CELL ON HOURLY BASIS

    Get DF w/ met (air T and rh) and water T data for cell

    :param cellLat:
    :param cellLong:
    :param rbmOutputDir:
    :param locPrecision:
    :param curtailmentYear:
    :param cellFoldername:
    :param resultsDir:
    :return:
    """
    # Load pandas DF w/ air T and rh and datetime just for current year
    metData = loadMetData(cellLat, cellLong, rbmOutputDir, locPrecision, curtailmentYear)
    # Load cell's avg water T data (2d list of Y-M-D, water T (vertical))
    hourlyWaterT = loadWaterData(cellLat, cellLong, rbmOutputDir, locPrecision, curtailmentYear, cellFoldername)
    assert (len(hourlyWaterT) == metData.shape[0])
    metData['waterC'] = pd.Series(hourlyWaterT)
    metData['waterF'] = metData.loc[:, 'waterC'] * 9 / 5 + 32
    metData['airF'] = metData.loc[:, 'tC'] * 9 / 5 + 32
    metData.to_csv(os.path.join(resultsDir, 'metAndWater' + str(cellLat) + '_' + str(cellLong) + '.csv'))
    return metData


def loadMetData(cellLat, cellLong, rbmOutputDir, locPrecision, curtailmentYear):
    """LOAD WATER T AND METEOROLOGICAL VARIABLES

    This should load cell-specific met data, but don't have that currently, so just load regional data.

    :param cellLat:
    :param cellLong:
    :param rbmOutputDir:
    :param locPrecision:
    :param curtailmentYear:
    :return:
    """
    metData = pd.read_csv('/Users/kiko/Documents/CE/AreaTVACellsallCurtailEnvRegsCnoneSnor//' +
                          'demandAndMetDfS_C_TVA' + str(curtailmentYear) + '.csv')
    return metData


###### LOAD WATER T ############################################################
def loadWaterData(cellLat, cellLong, rbmOutputDir, locPrecision, curtailmentYear, cellFoldername):
    cellTemperature = loadCellAvgWaterT(cellLat, cellLong, rbmOutputDir, locPrecision, cellFoldername)
    waterTInCurtailYear = getCellWaterTsInCurtailYear(cellTemperature, curtailmentYear)
    hourlyWaterT = list(np.array([[val] * 24 for val in waterTInCurtailYear]).flatten())
    return hourlyWaterT


def loadCellAvgWaterT(cellLat, cellLong, rbmOutputDir, locPrecision, cellFoldername):
    cellFolder = os.path.join(rbmOutputDir, cellFoldername)
    averageTFilename = createAverageTFilename(locPrecision, cellLat, cellLong)
    cellTemperature = readCSVto2dList(os.path.join(cellFolder, averageTFilename))
    return cellTemperature


# Looking up gen by gen lat/lon. May not be cell @ corresponding location.
# def loadDummyWaterTData(rbmOutputDir,allCellFolders,locPrecision):
#     fakeWaterT = 2
#     dummyFolder = allCellFolders[0]
#     (cellLat,cellLong) = getCellLatAndLongFromFolderName(dummyFolder)
#     averageTFilename = createAverageTFilename(locPrecision,cellLat,cellLong)
#     dummyCellT = readCSVto2dList(os.path.join(rbmOutputDir,dummyFolder,averageTFilename))
#     (dateCol,waterTCol) = (dummyCellT[0].index('Datetime'),dummyCellT[0].index('AverageWaterT(degC)'))
#     for row in dummyCellT[1:]: row[waterTCol] = fakeWaterT
#     return dummyCellT

def getCellLatAndLongFromFolderName(dummyFolder):
    dividerIdx = dummyFolder.index('_')
    return (float(dummyFolder[:dividerIdx]), float(dummyFolder[dividerIdx + 1:]))


###### ISOLATE DATA FOR YEAR OF ANALYSIS #######################################
# Return 2 1d lists of dates and water Ts in year of analysis. Daily basis.
def getCellWaterTsInCurtailYear(cellTemperature, curtailmentYear):
    (dateCol, waterTCol) = (cellTemperature[0].index('Datetime'), cellTemperature[0].index('AverageWaterT(degC)'))
    rowsInCurtailmentYear = [row for row in cellTemperature[1:]
                             if int(getElementsOfDate(row[dateCol])[0]) == curtailmentYear]  # skip header row
    datesInCurtailmentYear = [row[dateCol] for row in rowsInCurtailmentYear]
    temperaturesInCurtailmentYear = [float(row[waterTCol]) for row in rowsInCurtailmentYear]
    return temperaturesInCurtailmentYear