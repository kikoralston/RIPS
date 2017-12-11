#Michael Craig, 17 May 2017
#Functions load cell water T data from cells eligible for new generators
#by first setting  the folder names with future data for cells (each folder = 1 cell),
#then loading water temperatures from those folders.

import os, copy    
from AuxFuncs import *
from AssignCellsToIPMZones import mapCellToIPMZone
from ModifyGeneratorCapacityWithWaterTData import (getGenToCellAndCellToGenDictionaries,
                                        createBaseFilenameToReadOrWrite)


########### GET DAILY AVERAGE WATER TS FOR CELLS FOR ALL YEARS #################
#Current options for which cells can hsot new plants: 'All' (get all cell data) 
#or 'WithGen' (get data for cells that have generator inside at beginning of CE runs).
def loadEligibleCellWaterTs(cellsEligibleForNewPlants,genFleet,rbmOutputDir,
            locPrecision,zones,fipsToZones,fipsToPolys):
    allCellFoldersInZone,eligibleCellFolders = setCellFolders(rbmOutputDir,
        cellsEligibleForNewPlants,genFleet,locPrecision,zones,fipsToZones,fipsToPolys)
    return loadCellWaterTs(eligibleCellFolders,allCellFoldersInZone,rbmOutputDir)

#Get list of folders with future data for 1 cell each
def setCellFolders(rbmOutputDir,cellsEligibleForNewPlants,genFleet,locPrecision,
                    zones,fipsToZones,fipsToPolys):
    allCellFolders = os.listdir(rbmOutputDir)
    allCellFoldersInZone = isolateCellsInZones(allCellFolders,zones,fipsToZones,fipsToPolys)
    if cellsEligibleForNewPlants=='all':
        eligibleCellFolders = copy.deepcopy(allCellFoldersInZone)
    elif cellsEligibleForNewPlants=='withGens':
        # (genToCellLatLongsDictValues,genToCellLatLongsDict,cellLatLongToGenDict,
        #             genToCellLatLongsList) = getGenToCellAndCellToGenDictionaries(genFleet)
        (genToCellLatLongsDict,cellLatLongToGenDict,
                    genToCellLatLongsList) = getGenToCellAndCellToGenDictionaries(genFleet)
        eligibleCellFolders = []
        for (cellLat,cellLong) in cellLatLongToGenDict:
            eligibleCellFolders.append(createBaseFilenameToReadOrWrite(locPrecision,cellLat,cellLong))
    return allCellFoldersInZone,eligibleCellFolders

#Isolates cells in zones of analysis and returns list of cell folders in zones
def isolateCellsInZones(allCellFolders,zones,fipsToZones,fipsToPolys):
    allCellFoldersInZone = list()
    for cell in allCellFolders:
        cellZone = mapCellToIPMZone(cell,fipsToZones,fipsToPolys)
        if cellZone in zones: allCellFoldersInZone.append(cell)
    return allCellFoldersInZone

#Returns dict of cell folder name : [[Datetime],[AverageWaterT(degC)]]
def loadCellWaterTs(eligibleCellFolders,allCellFoldersInZone,rbmOutputDir): 
    eligibleCellWaterTs = dict()
    for cellFolder in eligibleCellFolders:
        if cellFolder in allCellFoldersInZone:
            cellFolderAvgTFilename = cellFolder + 'Average.csv'
            avgTFilePath = os.path.join(rbmOutputDir,cellFolder,cellFolderAvgTFilename)
            cellAvgWaterT = readCSVto2dList(avgTFilePath)
            eligibleCellWaterTs[cellFolder] = cellAvgWaterT #[[Datetime],[AverageWaterT(degC)]]
    return eligibleCellWaterTs
