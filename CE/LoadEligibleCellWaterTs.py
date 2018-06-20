""""
Michael Craig, 17 May 2017
Functions load cell water T data from cells eligible for new generators by first setting
the folder names with future data for cells (each folder = 1 cell), then loading water temperatures
from those folders.
"""

import os, copy
from AuxFuncs import *
from AssignCellsToIPMZones import mapCellToIPMZone
from ModifyGeneratorCapacityWithWaterTData import (getGenToCellAndCellToGenDictionaries,
                                        createBaseFilenameToReadOrWrite)


def loadEligibleCellWaterTs(cellsEligibleForNewPlants, genFleet, rbmOutputDir, locPrecision, zones,
                            fipsToZones, fipsToPolys, currYear):
    """GET DAILY AVERAGE WATER TS FOR CELLS FOR ALL YEARS

    Current options for which cells can hsot new plants: 'All' (get all cell data)
    or 'WithGen' (get data for cells that have generator inside at beginning of CE runs).

    :param cellsEligibleForNewPlants:
    :param genFleet:
    :param rbmOutputDir:
    :param locPrecision:
    :param zones:
    :param fipsToZones:
    :param fipsToPolys:
    :param currYear:
    :return:
    """
    allCellFoldersInZone, eligibleCellFolders = setCellFolders(rbmOutputDir, cellsEligibleForNewPlants, genFleet,
                                                              locPrecision, zones, fipsToZones, fipsToPolys)
    return loadCellWaterTs(eligibleCellFolders, allCellFoldersInZone, rbmOutputDir, currYear)


def setCellFolders(rbmOutputDir, cellsEligibleForNewPlants, genFleet, locPrecision, zones, fipsToZones,
                   fipsToPolys):
    """Get list of folders with future data for 1 cell each

    :param rbmOutputDir:
    :param cellsEligibleForNewPlants:
    :param genFleet:
    :param locPrecision:
    :param zones:
    :param fipsToZones:
    :param fipsToPolys:
    :return:
    """
    #allCellFolders = os.listdir(rbmOutputDir)
    allCellFolders = [name for name in os.listdir(rbmOutputDir) if os.path.isdir(os.path.join(rbmOutputDir, name))]
    allCellFoldersInZone = isolateCellsInZones(allCellFolders, zones, fipsToZones, fipsToPolys)
    if cellsEligibleForNewPlants == 'all':
        eligibleCellFolders = copy.deepcopy(allCellFoldersInZone)
    elif cellsEligibleForNewPlants == 'withGens':
        # (genToCellLatLongsDictValues,genToCellLatLongsDict,cellLatLongToGenDict,
        #             genToCellLatLongsList) = getGenToCellAndCellToGenDictionaries(genFleet)
        (genToCellLatLongsDict,cellLatLongToGenDict,
                    genToCellLatLongsList) = getGenToCellAndCellToGenDictionaries(genFleet)
        eligibleCellFolders = []
        for (cellLat,cellLong) in cellLatLongToGenDict:
            eligibleCellFolders.append(createBaseFilenameToReadOrWrite(locPrecision,cellLat,cellLong))
    return allCellFoldersInZone, eligibleCellFolders


def isolateCellsInZones(allCellFolders, zones, fipsToZones, fipsToPolys):
    """Isolates cells in zones of analysis and returns list of cell folders in zones

    :param allCellFolders:
    :param zones:
    :param fipsToZones:
    :param fipsToPolys:
    :return:
    """
    allCellFoldersInZone = list()
    for cell in allCellFolders:
        cellZone = mapCellToIPMZone(cell, fipsToZones, fipsToPolys)
        if cellZone in zones: allCellFoldersInZone.append(cell)
    return allCellFoldersInZone


def loadCellWaterTs(eligibleCellFolders, allCellFoldersInZone, rbmOutputDir, currYear=None):
    """Returns dict of cell folder name : [[Datetime],[AverageWaterT(degC)]]

    :param eligibleCellFolders:
    :param allCellFoldersInZone:
    :param rbmOutputDir:
    :param currYear: (integer) year of data being imported. If None, imports all years
    :return:
    """
    eligibleCellWaterTs = dict()
    for cellFolder in eligibleCellFolders:
        if cellFolder in allCellFoldersInZone:

            cellFolderAvgTFilename = cellFolder + 'Average.csv'
            avgTFilePath = os.path.join(rbmOutputDir, cellFolder, cellFolderAvgTFilename)
            cellAvgWaterT = readCSVto2dList(avgTFilePath)

            # filter current year
            if currYear is not None:
                dateCol = cellAvgWaterT[0].index('date')
                waterTCol = cellAvgWaterT[0].index('waterT')
                airTCol = cellAvgWaterT[0].index('AirT')
                flowCol = cellAvgWaterT[0].index('flow')

                # convert to numeric
                cellAvgWaterTCurrYear = [cellAvgWaterT[0]]

                for row in cellAvgWaterT[1:]:
                    if str(currYear) in row[dateCol]:
                        cellAvgWaterTCurrYear = cellAvgWaterTCurrYear + \
                                                [[row[dateCol], float(row[waterTCol]), float(row[airTCol]),
                                                  float(row[flowCol])]]
            else:
                cellAvgWaterTCurrYear = cellAvgWaterT

            eligibleCellWaterTs[cellFolder] = cellAvgWaterTCurrYear

    return eligibleCellWaterTs
