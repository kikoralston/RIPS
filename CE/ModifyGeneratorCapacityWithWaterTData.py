"""
Michael Craig
August 11, 2016
Script for modifying generator capacity with water T from RBM.

1) Gets list of all grid cells from .spat.
2) Saves T data for each segment in each grid cell into folders for each cell.
3) Averages T data across all segments for each cell and saves that data in folder for each cell.
4) Creates dictionaries mapping generators to cells and vice versa.
5) Using dictionary, pairs generator with water T time series and modifies capacity.
"""

import os, csv, copy
import numpy as np
import pandas as pd
import datetime as dt
import time
from AuxFuncs import *
from GAMSAuxFuncs import *
from CurtailmentRegressions import (calcCurtailmentForGenOrTech, loadRegCoeffs, getKeyCurtailParams,
                                    getCoeffsForGenOrTech)


################################################################################
####### MASTER FUNCTIONS #######################################################
################################################################################
def processRBMDataIntoIndividualCellFiles(curtailparam):
    """Main routine to pre process RBM Raw data to individual cell files

    This routine reads the raw output files from a RBM simulations and saves the relevant data into separate CSV
    files which will be used by the CE/UC models.

    :param curtailparam: object of class Curtailmentparameters
    """
    gridCellLatLongs = getAllGridCellLatLongsInSpatFile(curtailparam.rbmDataDir, curtailparam.tempAndSpatFilename)
    saveAndAverageTDataForAllCells(gridCellLatLongs, curtailparam)


def determineHrlyCurtailmentsForExistingGens(genFleet, currYear, modelName, genparam, curtailparam, resultsDir):
    """Computes times series of hourly curtailments for EXISTING generators

    :param genFleet:
    :param currYear:
    :param modelName:
    :param genparam:
    :param curtailparam:
    :param resultsDir:
    :return: dictionaryof ORIS+UNITID to 2d vertical list of [datetime,curtailment(mw)],
             list of (ORIS+UNIT,gen lat long, cell lat long), and
             list of hourly curtailments for existing gens.
    """
    (genToCellLatLongsDict, cellLatLongToGenDict,
     genToCellLatLongsList) = getGenToCellAndCellToGenDictionaries(genFleet)
    hrlyCurtailmentsAllGensInTgtYr, hrlyCurtailmentsList = \
        calculateGeneratorCurtailments(cellLatLongToGenDict, currYear, genFleet, modelName, genparam,
                                       curtailparam, resultsDir)
    return hrlyCurtailmentsAllGensInTgtYr, genToCellLatLongsList, hrlyCurtailmentsList


#####################################################################################################
# ----------- ROUTINES USED IN THE PRE-PROCESSING OF RBM DATA INTO INDIVIDUAL CELL FILES -----------
#####################################################################################################

def getAllGridCellLatLongsInSpatFile(rbmDataDir, tempAndSpatFilename):
    """Get all grid cell lat/longs in spat file

    SPAT FILE COLUMNS:
    reach, cell, row #, column #, lat, lon, segment # (within the cell).
    lat/lon are in decimal degree.
    Note each cell has two segments - so it's the same lat/lon for reach segment (i know - kind of redundant).

    :param rbmDataDir:
    :param tempAndSpatFilename:
    :return:
    """
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


def saveAndAverageTDataForAllCells(gridCellLatLongs, curtailparam):
    """Average temperature data for all segments in each cell, then save all data

    :param gridCellLatLongs:
    :param curtailparam:
    """

    start_time = time.time()

    print('TOTAL number of grid cells in data set: {0:8d}\n'.format(len(gridCellLatLongs)))
    (cellCtr, totalCellCtr, tgtCellLatLongs) = (0, 0, [])
    for (cellLat, cellLong) in gridCellLatLongs:
        batch_start = time.time()
        cellCtr += 1
        totalCellCtr += 1

        # first compile list of cells to process
        tgtCellLatLongs.append((cellLat, cellLong))
        if cellCtr == curtailparam.numCellsToProcess or totalCellCtr == len(gridCellLatLongs):
            # create files with processed data from list of cells
            print('Number of grid cells being processed on this batch: {0:8d}'.format(len(tgtCellLatLongs)))
            saveAndAverageTDataForGridCells(tgtCellLatLongs, curtailparam)

            (cellCtr, tgtCellLatLongs) = (0, [])

            elapsed_time_batch = time.time() - batch_start
            elapsed_time_tot = time.time() - start_time
            h1 = int(elapsed_time_batch / (60 * 60))
            m1 = int((elapsed_time_batch % (60 * 60)) / 60)
            s1 = int(elapsed_time_batch % 60.)

            ht = int(elapsed_time_tot / (60 * 60))
            mt = int((elapsed_time_tot % (60 * 60)) / 60)
            st = int(elapsed_time_tot % 60.)

            print('\nElapsed time: Batch = {0:>02d}m:{1:>02d}s; Total = {2:03d}h:{3:>02d}m:{4:>02d}s'.format(m1, s1, ht,
                                                                                                             mt, st))
            print('--- Total Cell count: {} ---\n'.format(totalCellCtr))


def saveAndAverageTDataForGridCells(tgtCellLatLongs, curtailparam, initialyear=2000):
    """

    :param tgtCellLatLongs:
    :param curtailparam:
    :param initialyear: (integer) Data from previous years will not be saved in the files.

    """
    outputDirs = createOutputDirs(tgtCellLatLongs, curtailparam)
    (numTotalSegments, numDays) = getNumSegmentsAndDays(curtailparam.rbmDataDir, curtailparam.nsegFilename)
    (cellInfo, mapSegmentNumToLatLong, mapLatLongToSegmentNum) = getLinesInSpatFileForCell(tgtCellLatLongs,
                                                                                           numTotalSegments,
                                                                                           curtailparam)
    (allSegmentData, waterTAllSegments) = readCellTemperatureData(cellInfo, numTotalSegments, curtailparam)

    writeCellTemperatureData(allSegmentData, cellInfo, outputDirs, curtailparam.locPrecision,
                             numDays, curtailparam.outputHeaders, mapSegmentNumToLatLong)

    averageAndSaveWaterTOverAllSegmentsInCells(waterTAllSegments, outputDirs, curtailparam.locPrecision,
                                               mapLatLongToSegmentNum, initialyear)


def createOutputDirs(tgtCellLatLongs, curtailparam):
    """CREATES OUTPUT DIRECTORIES FOR CELL T DATA BEING PROCESSED

    :param tgtCellLatLongs: list with (lat, lon) tuples for cells being processed
    :param curtailparam: object of class Curtailmentparameters
    :return: dictionary with (cellLat, cellLong): 'path/to/folder'
    """

    print('Creating folders.... ', end="")
    outputDirs = dict()
    for (cellLat, cellLong) in tgtCellLatLongs:
        outputDir = os.path.join(curtailparam.rbmOutputDir,
                                 '%.*f_%.*f' % (curtailparam.locPrecision, cellLat,
                                                curtailparam.locPrecision, cellLong))
        if not os.path.exists(outputDir): os.makedirs(outputDir)
        outputDirs[(cellLat, cellLong)] = outputDir
    print('Done!')

    return outputDirs


def getNumSegmentsAndDays(dataDir, nsegFilename):
    """EXTRACT TOTAL NUM SEGMENTS AND DAYS IN MODEL RUN

    This function reads the total number of segments and days in a model run. This information is in the a text file
    'nsegFilename'. The name of this file is defined in Curtailment parameters

    :param dataDir: string with location of nsegFilename.
    :param nsegFilename: string with name of nsegFilename. The name of this file is defined in Curtailment parameters.
    :return: tuple with number of segments and number of days
    """
    nseg_nday = np.loadtxt(os.path.join(dataDir, nsegFilename))
    (numTotalSegments, numDays) = (int(nseg_nday[0]), int(nseg_nday[1]))
    return (numTotalSegments, numDays)


def getLinesInSpatFileForCell(tgtCellLatLongs, numTotalSegments, curtailparam):
    """GET ALL LINES IN SPATIAL FILE THAT CORRESPOND TO CELL

    Gets line numbers in spatial file that match cell. Includes all reaches & segments
    EXCEPT segment 2 in outlet cell (= max segment #) b/c no data for that segment.

    :param tgtCellLatLongs:
    :param numTotalSegments:
    :param curtailparam:
    :return:
    """
    (mapSegmentNumToLatLong, mapLatLongToSegmentNum) = (dict(), dict())
    spatFile = os.path.join(curtailparam.rbmDataDir, curtailparam.tempAndSpatFilename + '.Spat')  # .Spat file
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
    return cellInfo, mapSegmentNumToLatLong, mapLatLongToSegmentNum


def readCellTemperatureData(cellInfo, numTotalSegments, curtailparam):
    """READ CELL TEMPERATURE DATA

    Temperature data: hour 1, all cells + segments + reaches, then hour 2, all cells + segments + reaches, etc.

    :param cellInfo:
    :param numTotalSegments:
    :param curtailparam:
    :return:
    """
    print('Reading .Temp file... ', end='', flush=True)

    tempFile = os.path.join(curtailparam.rbmDataDir, curtailparam.tempAndSpatFilename + '.Temp')  # .Temp file

    # dictionary. total segment number (=line #in .spat file): np array of year, month, day, flow(cfs), streamTemp(degC).
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
            (tDataDict, year, month, day, waterT) = readRawTemperatureLineData(tempLineRaw)
            tLineData = processTempDataDictIntoRow(tDataDict, curtailparam.outputHeaders)

            # Save entire line of segment data
            allSegmentData[totalSegmentNum].append(tLineData)

            # Save data of segment into nested dictionary - total segment number: date : (waterT, airT, flow)
            waterTAllSegments[totalSegmentNum][createDateLabel(year, month, day)] = (tDataDict['streamT'],
                                                                                     tDataDict['airT'],
                                                                                     tDataDict['flow'])

    f.close()

    for i in allSegmentData: allSegmentData[i] = np.asarray(allSegmentData[i])

    print('Done!')

    return allSegmentData, waterTAllSegments


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

    return tDataDict, year, date.month, date.day, streamT


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
    """Auxiliary function that maps header names to indexes of columns

    This is an auxiliary function used in 'processTempDataDictIntoRow' that maps the header names to indexes of columns

    :param cols: list of strings with names of variables in .Temp file. This list is defined in Curtailment parameters
    :return: tuple with the indexes of the variables in a specific arbitrary order (defined by the variables
             returned)
    """
    (yearCol, monthCol, dayCol, flowCol, streamTCol, headTCol, airTCol) = (cols.index('Year'),
                                                                           cols.index('Month'), cols.index('Day'),
                                                                           cols.index('Streamflow(cfs)'),
                                                                           cols.index('StreamT(degC)'),
                                                                           cols.index('HeadwaterT(degC)'),
                                                                           cols.index('AirT(degC)'))
    return (yearCol, monthCol, dayCol, flowCol, streamTCol, headTCol, airTCol)


def createDateLabel(year, month, day):
    return '%s-%s-%s' % (year, month, day)


def writeCellTemperatureData(allSegmentData, cellInfo, outputDirs, locPrecision, numDays,
                             outputHeaders, mapSegmentNumToLatLong):
    """WRITE FILES WITH CELL TEMPERATURE DATA

    :param allSegmentData:
    :param cellInfo:
    :param outputDirs: (string) path to base output folder
    :param locPrecision: (integer) number of decimal digits in lat and long
    :param numDays:
    :param outputHeaders: (string) header of output file
    :param mapSegmentNumToLatLong:
    """
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
    """Creates string with name of folder with cell data

    This function creates a string with the formatted name of the folder that contains the data
    for the respective grid cell

    :param locPrecision:
    :param inputLat: latitude of grid cell
    :param inputLong: longitude of grid cell
    :return: string with name of folder (e.g. 34.4375_-86.4375)
    """
    return '%.*f_%.*f' % (locPrecision, inputLat, locPrecision, inputLong)


def createHeaderStr(outputHeaders):
    headStr = ''
    for idx in range(len(outputHeaders)):
        headStr += outputHeaders[idx]
        if outputHeaders[idx] != outputHeaders[-1]: headStr += ' '  # not last header
    return headStr + '\n'


def averageAndSaveWaterTOverAllSegmentsInCells(waterTAllSegments, outputDirs, locPrecision,
                                               mapLatLongToSegmentNum, initialyear=2000):
    """AVERAGE AND SAVE CELL TEMPERATURE DATA TO FILE

    :param waterTAllSegments: dictionary with water data for all segments
    :param outputDirs: (string) base path to output folder
    :param locPrecision: (integer) number of decimal digits in lat and long values
    :param mapLatLongToSegmentNum: (dictionary) mapping of (lat, lon) to segments
    :param initialyear: (integer) data from previous years will not be saved to files
    """
    for (cellLat, cellLon) in mapLatLongToSegmentNum:
        print('Processing and saving data for cell ({lat:.{p}f},{lon:.{p}f})... '.format(lat=cellLat, lon=cellLon,
                                                                                         p=locPrecision), end='',
              flush=True)
        currCellSegmentNums = mapLatLongToSegmentNum[(cellLat, cellLon)]

        # create subdictionary by isolating data for segments in currCellSegmentNums
        waterTAllSegmentsInCell = {k: waterTAllSegments[k] for k in currCellSegmentNums}

        # average data over all segment in cell
        waterTAvgInCell = averageWaterTOverAllSegmentsInCell(waterTAllSegmentsInCell)

        # convert to PD Data Frame
        waterTAvgInCell = pd.DataFrame.from_dict(waterTAvgInCell, orient='index')
        waterTAvgInCell.rename({0: 'waterT', 1: 'AirT', 2: 'flow'}, inplace=True, axis='columns')
        waterTAvgInCell.index.name = 'date'
        waterTAvgInCell.reset_index(inplace=True)
        waterTAvgInCell['date'] = pd.to_datetime(waterTAvgInCell['date'])
        # sort according to date
        waterTAvgInCell = waterTAvgInCell.sort_values(by=['date'])
        waterTAvgInCell = waterTAvgInCell.reset_index(drop=True)
        # filter years
        waterTAvgInCell = waterTAvgInCell[(waterTAvgInCell['date'].dt.year >= initialyear)]
        waterTAvgInCell['date'] = waterTAvgInCell['date'].dt.strftime('%Y-%m-%d')

        # convert to numpy array for saving to file
        namecols = np.array(waterTAvgInCell.keys()).reshape(1, 4)
        waterTAvgInCell = np.array(waterTAvgInCell)

        waterTAvgInCell = np.concatenate((namecols, waterTAvgInCell), axis=0)

        outputDir = outputDirs[(cellLat, cellLon)]
        saveAverageWaterT(waterTAvgInCell, outputDir, locPrecision, cellLat, cellLon)

        print('Done!')


def averageWaterTOverAllSegmentsInCell(waterTAllSegments):
    waterTSums = dict()
    for segment in waterTAllSegments:
        waterTSegment = waterTAllSegments[segment]
        for date in waterTSegment:
            waterTSums[date] = np.add(waterTSums.get(date, 0), waterTSegment[date])
    numSegments = len(waterTAllSegments)
    for date in waterTSums: waterTSums[date] = waterTSums[date] / numSegments
    return waterTSums


def getElementsOfDate(listDate):
    year = int(listDate[:4])
    restOfDate = listDate[5:]
    monthEndIdx = restOfDate.index('-')
    month = int(restOfDate[:monthEndIdx])
    day = int(restOfDate[monthEndIdx + 1:])
    return (year, month, day)


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

    Takes the 2d list with generator fleet data and creates dictionaries mapping individual plants to grid cells
    lat and long

    :param genFleet: 2d list with generator fleet data
    :return: genToCellLatLongsDict: dictionary gen:[gen,(gen lat, gen long), (cell lat, cell long)]
             cellLatLongToGenDict:  dictionary (cell lat, clell long):genID.
             genToCellLatLongsList: 2d list [[genID,(genlat,genlong),(celllat,celllong)]]

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

    return genToCellLatLongsDict, cellLatLongToGenDict, genToCellLatLongsList


def find125GridMaurerLatLong(lat, lon):
    """Get (lat,lon) of 1/8 grid cell that a (lat, lon) point falls in

    :param lat: latitude value
    :param lon:  longitude value
    :return: tuple (lat_grid, long_grid) of grid cell that contains the point (lat,lon)
    """

    lat_grid = np.around(8.0 * lat - 0.5) / 8.0 + 0.0625
    lon_grid = np.around(8.0 * lon - 0.5) / 8.0 + 0.0625
    return lat_grid, lon_grid


################################################################################
################################################################################
def loadWaterAndMetData(cellLat, cellLong, rbmOutputDir, locPrecision, curtailmentYear,
                        cellFoldername, resultsDir):
    """LOAD WATER AND MET DATA BY CELL ON HOURLY BASIS

    Read text files with data for met (air T and rh) and water T data for one specific cell. The data files were
    created by the pre-processing function 'processRBMDataIntoIndividualCellFiles' and saved in a folder 'rbmOutputDir'.

    :param cellLat: latitude of grid cell
    :param cellLong: longitude of grid cell
    :param rbmOutputDir: string with folder of pre-processed RBM (Stream temperature, Stream flow, etc...) data. This
                         pre-processed data is created by the function 'processRBMDataIntoIndividualCellFiles'
    :param locPrecision: integer with number of decimal digits in latitude and longitude
    :param curtailmentYear: integer with year that curtailment is being simulated
    :param cellFoldername:
    :param resultsDir: string with folder of output data of CE model
    :return: panda data frame
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


def calculateGeneratorCurtailments(cellLatLongToGenDict, curtailmentYear, genFleet, modelName, genparam,
                                   curtailparam, resultsDir):
    """CURTAIL GENERATOR CAPACITY WITH WATER TEMPERATURES

    Calculates generator curtailments. If generator isn't curtailed (not right plant type, cell data not avaialble,
    etc.), ignored here. CalculateHourlyCapacs... script handles those cases (assume curtailment = 0).
    cellLatLongToGenDict = dict of (cell lat, cell long):gen ID for all gens w/ lat/lon coords in genFleet,
    including gens that should nto be curtailed.

    :param cellLatLongToGenDict:
    :param curtailmentYear:
    :param genFleet:
    :param modelName:
    :param genparam:
    :param curtailparam:
    :param resultsDir:
    :return:
    """
    hrlyCurtailmentsAllGensInTgtYr, hrlyCurtailmentsList = dict(), list()
    regCoeffs = loadRegCoeffs(genparam.dataRoot, 'capacity.json')  # dict of planttype:coolingtype:cooldesignT:param:coeffs
    allCellFolders = os.listdir(curtailparam.rbmOutputDir)

    # this maps gen lat/lon to gen IDs; cell data may not exist
    for (cellLat, cellLong) in cellLatLongToGenDict:
        cellFoldername = createBaseFilenameToReadOrWrite(curtailparam.locPrecision, cellLat, cellLong)
        if cellFoldername in allCellFolders:
            metAndWaterData = loadWaterAndMetData(cellLat, cellLong, curtailparam.rbmOutputDir,
                                                  curtailparam.locPrecision, curtailmentYear, cellFoldername,
                                                  resultsDir)
            gensInCell = cellLatLongToGenDict[(cellLat, cellLong)]  # list of ORIS-UNITID in cell
            for gen in gensInCell:
                (plantType, hr, fuelAndCoalType, coolType, fgdType, state, capac) = getKeyCurtailParams(gen, genFleet)
                coeffs = getCoeffsForGenOrTech(plantType, coolType, genparam.ptCurtailed, regCoeffs,
                                               genparam.coolDesignT)
                if (coeffs is not None) and (plantType in genparam.ptCurtailedRegs):  # skip gens that aren't curtailed
                    hrlyCurtailmentsGen = calcCurtailmentForGenOrTech(plantType, hr, fuelAndCoalType, coolType, fgdType,
                                                                      state, capac, metAndWaterData, coeffs,
                                                                      genparam, curtailparam)
                    hrlyCurtailmentsAllGensInTgtYr[gen] = hrlyCurtailmentsGen
                    hrlyCurtailmentsList.append([gen] + hrlyCurtailmentsGen)
        else:
            print('Cell not in folders!:', cellFoldername)
    return hrlyCurtailmentsAllGensInTgtYr, hrlyCurtailmentsList


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


def loadWaterData(cellLat, cellLong, rbmOutputDir, locPrecision, curtailmentYear, cellFoldername):
    """LOAD WATER TEMPERATURE DATA

    This function reads file with water data simulated by RBM and returns a list with the data

    :param cellLat:
    :param cellLong:
    :param rbmOutputDir:
    :param locPrecision:
    :param curtailmentYear:
    :param cellFoldername:
    :return:
    """
    cellTemperature = loadCellAvgWaterT(cellLat, cellLong, rbmOutputDir, locPrecision, cellFoldername)
    waterTInCurtailYear = getCellWaterTsInCurtailYear(cellTemperature, curtailmentYear)
    hourlyWaterT = list(np.array([[val] * 24 for val in waterTInCurtailYear]).flatten())
    return hourlyWaterT


def loadCellAvgWaterT(cellLat, cellLong, rbmOutputDir, locPrecision, cellFoldername):
    """

    :param cellLat: latitide of grid cell
    :param cellLong:  longitude of grid cell
    :param rbmOutputDir: string with path to folder with processed RBM data
    :param locPrecision:  integer with number of decimal digits in logintude and latitude values
    :param cellFoldername:
    :return:
    """
    cellFolder = os.path.join(rbmOutputDir, cellFoldername)
    averageTFilename = createAverageTFilename(locPrecision, cellLat, cellLong)
    cellTemperature = readCSVto2dList(os.path.join(cellFolder, averageTFilename))
    return cellTemperature


def getCellLatAndLongFromFolderName(dummyFolder):
    """

    :param dummyFolder:
    :return:
    """
    dividerIdx = dummyFolder.index('_')
    return (float(dummyFolder[:dividerIdx]), float(dummyFolder[dividerIdx + 1:]))


def getCellWaterTsInCurtailYear(cellTemperature, curtailmentYear):
    """ISOLATE DATA FOR YEAR OF ANALYSIS

    Return 2 1d lists of dates and water Ts in year of analysis. Daily basis.

    :param cellTemperature:
    :param curtailmentYear: year being analyzed (integer)
    :return: 2 1d lists of dates and water Ts in year of analysis. Daily basis.
    """
    (dateCol, waterTCol) = (cellTemperature[0].index('Datetime'), cellTemperature[0].index('AverageWaterT(degC)'))
    rowsInCurtailmentYear = [row for row in cellTemperature[1:]
                             if int(getElementsOfDate(row[dateCol])[0]) == curtailmentYear]  # skip header row
    datesInCurtailmentYear = [row[dateCol] for row in rowsInCurtailmentYear]
    temperaturesInCurtailmentYear = [float(row[waterTCol]) for row in rowsInCurtailmentYear]
    return temperaturesInCurtailmentYear


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
