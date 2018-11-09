""""
Michael Craig, 17 May 2017
Functions load cell water T data from cells eligible for new generators by first setting
the folder names with future data for cells (each folder = 1 cell), then loading water temperatures
from those folders.

August 2018
UW changed the format of the files with water data to NETCDF. I added a flag 'netcdf' to some functions to read them.
"""

import os, copy, sys
import numpy as np
from AuxFuncs import *
import pickle as pk
from AssignCellsToIPMZones import mapCellToIPMZone
from ModifyGeneratorCapacityWithWaterTData import getGenToCellAndCellToGenDictionaries
from AuxCurtailmentFuncs import get_all_cells_from_netcdf, order_cells_by_flow, get_all_cells_in_zone, loadCellWaterTs
from PreProcessRBM import createBaseFilenameToReadOrWrite

# genparam.cellsEligibleForNewPlants, genFleet, curtailparam.rbmOutputDir, curtailparam.locPrecision, genparam.ipmZones, genparam.fipsToZones, genparam.fipsToPolys, currYear
# cellsEligibleForNewPlants, genFleet, rbmOutputDir, locPrecision, zones, fipsToZones, fipsToPolys, currYear


def loadEligibleCellWaterTs(genFleet, currYear, genparam, curtailparam):
    """GET DAILY AVERAGE WATER TS FOR CELLS FOR ALL YEARS

    Current options for which cells can hsot new plants: 'All' (get all cell data)
    or 'WithGen' (get data for cells that have generator inside at beginning of CE runs).

    :param genFleet:
    :param currYear:
    :param genparam:
    :param curtailparam:
    :return: dict of {cell folder name : [[Datetime],[AverageWaterT(degC)], [AirT], [flow]]}
    """
    allCellFoldersInZone, eligibleCellFolders = setCellFolders(genFleet, currYear, genparam, curtailparam)

    dict_out = dict()
    for gcm in curtailparam.listgcms:
        dict_out[gcm] = loadCellWaterTs(eligibleCellFolders, allCellFoldersInZone, curtailparam, gcm, currYear)

    return dict_out


# rbmOutputDir, cellsEligibleForNewPlants, genFleet, locPrecision, zones, fipsToZones, fipsToPolys
# rbmOutputDir, cellsEligibleForNewPlants, genFleet, locPrecision, zones, fipsToZones, fipsToPolys, netcdf=True

def setCellFolders(genFleet, currYear, genparam, curtailparam, netcdf=True):
    """Get list of folders with future data for 1 cell each

    :param genFleet: 2d list with generator fleet
    :param netcdf: (boolean) if true reads water data from NETCDF files (new format from August 2018)
    :return:
    """

    if netcdf:
        # get cells from first gcm in list
        gcm = curtailparam.listgcms[0]

        allCellFolders = get_all_cells_from_netcdf(curtailparam, gcm)
        allCellFoldersInZone = get_all_cells_in_zone(allCellFolders, genparam, curtailparam)

    else:
        allCellFolders = [name for name in os.listdir(curtailparam.rbmOutputDir)
                          if os.path.isdir(os.path.join(curtailparam.rbmOutputDir, name))]

        allCellFoldersInZone = isolateCellsInZones(allCellFolders, genparam)

    if genparam.cellsEligibleForNewPlants == 'all':
        eligibleCellFolders = copy.deepcopy(allCellFoldersInZone)

    elif genparam.cellsEligibleForNewPlants == 'maxflow':
        # get all cells in each zone ordered by flow (descending order)
        best_cells_zone_dict = order_cells_by_flow(genparam, curtailparam, currYear, sys.maxsize)

        # get list of cells with existing generators in it (including build decisions from previous CE runs)
        cellLatLongToGenDict = getGenToCellAndCellToGenDictionaries(genFleet)[1]
        cells_list = list(cellLatLongToGenDict.keys())
        cells_list_str = list(map(lambda x: createBaseFilenameToReadOrWrite(curtailparam.locPrecision, x[0], x[1]),
                                  cells_list))

        # for each zone remove cells that already have a generator in it and keep only 100 best cells
        for z in best_cells_zone_dict:
            auxList = [x for x in best_cells_zone_dict[z] if x not in cells_list_str]
            auxList = auxList[:100]
            best_cells_zone_dict[z] = auxList

        # combine into one single list with all eligible cells
        eligibleCellFolders = []
        for z in best_cells_zone_dict:
            eligibleCellFolders = eligibleCellFolders + best_cells_zone_dict[z]

    elif genparam.cellsEligibleForNewPlants == 'withGens':
        # (genToCellLatLongsDictValues,genToCellLatLongsDict,cellLatLongToGenDict,
        #             genToCellLatLongsList) = getGenToCellAndCellToGenDictionaries(genFleet)
        (genToCellLatLongsDict, cellLatLongToGenDict, genToCellLatLongsList) = \
            getGenToCellAndCellToGenDictionaries(genFleet)

        eligibleCellFolders = []
        for (cellLat, cellLong) in cellLatLongToGenDict:
            eligibleCellFolders.append(createBaseFilenameToReadOrWrite(curtailparam.locPrecision, cellLat, cellLong))

    else:
        print('Parameter \'cellsEligibleForNewPlants\' must be either \'all\' or \'maxflow\' or \'withGens\'....')
        print('Ending simulation!')
        sys.exit()

    return allCellFoldersInZone, eligibleCellFolders


def isolateCellsInZones(allCellFolders, genparam):
    """Isolates cells in zones of analysis and returns list of cell folders in zones

    :param allCellFolders:
    :param zones:
    :param fipsToZones:
    :param fipsToPolys:
    :return:
    """
    allCellFoldersInZone = list()
    for cell in allCellFolders:
        cellZone = mapCellToIPMZone(cell, genparam.fipsToZones, genparam.fipsToPolys)
        if cellZone in genparam.impZones: allCellFoldersInZone.append(cell)
    return allCellFoldersInZone