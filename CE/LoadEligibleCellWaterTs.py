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
from ModifyGeneratorCapacityWithWaterTData import (getGenToCellAndCellToGenDictionaries, read_waterdata_netcdf,
                                                   get_all_cells_from_netcdf, order_cells_by_flow,
                                                   getCellLatAndLongFromFolderName)
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

    return loadCellWaterTs(eligibleCellFolders, allCellFoldersInZone, curtailparam, currYear)


# rbmOutputDir, cellsEligibleForNewPlants, genFleet, locPrecision, zones, fipsToZones, fipsToPolys
# rbmOutputDir, cellsEligibleForNewPlants, genFleet, locPrecision, zones, fipsToZones, fipsToPolys, netcdf=True

def setCellFolders(genFleet, currYear, genparam, curtailparam, netcdf=True):
    """Get list of folders with future data for 1 cell each

    :param genFleet: 2d list with generator fleet
    :param netcdf: (boolean) if true reads water data from NETCDF files (new format from August 2018)
    :return:
    """

    if netcdf:

        allCellFolders = get_all_cells_from_netcdf(curtailparam)

        # read file with dictionary mapping cells to IPM zones (this file must be in the VICS/RBM folder)
        with open(os.path.join(curtailparam.rbmRootDir, 'cells2zones.pk'), 'rb') as f:
            cells2zones = pk.load(f)

        allCellFoldersInZone = [c for c in allCellFolders if cells2zones[c] in genparam.ipmZones]
    else:
        allCellFolders = [name for name in os.listdir(curtailparam.rbmOutputDir)
                          if os.path.isdir(os.path.join(curtailparam.rbmOutputDir, name))]

        allCellFoldersInZone = isolateCellsInZones(allCellFolders, genparam)

    if genparam.cellsEligibleForNewPlants == 'all':
        eligibleCellFolders = copy.deepcopy(allCellFoldersInZone)

    elif genparam.cellsEligibleForNewPlants == 'maxflow':
        best_cells_zone_dict = order_cells_by_flow(genparam, curtailparam, currYear)

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


def loadCellWaterTs(eligibleCellFolders, allCellFoldersInZone, curtailparam, currYear=None, netcdf=True):
    """Returns dict of {cell folder name : [[Datetime], [AverageWaterT(degC)], [AirT], [flow]]}

    :param eligibleCellFolders:
    :param allCellFoldersInZone:
    :param curtailparam: object of class Curtailparameters
    :param currYear: (integer) year of data being imported. If None, imports all years
    :param netcdf: (boolean) if true reads water data from NETCDF files (new format from August 2018)

    :return: dict of {cell folder name : [[Datetime],[AverageWaterT(degC)], [AirT], [flow]]}
    """
    eligibleCellWaterTs = dict()

    if netcdf:
        # TODO: find a way to combine different GCMs
        gcm = curtailparam.listgcms[0]

        # reading flow/streamT file. substitute '_' to '.'
        gcm = gcm.replace('_', '.')
        gcm = gcm.replace('rcp', 'RCP')

        fname = curtailparam.basenamestreamT
        fname = os.path.join(curtailparam.rbmDataDir, fname.format(gcm))
        waterT, date, lons, lats = read_waterdata_netcdf(fname, 'T_stream', currYear)

        fname = curtailparam.basenameflow
        fname = os.path.join(curtailparam.rbmDataDir, fname.format(gcm))
        streamflow = read_waterdata_netcdf(fname, 'streamflow', currYear)[0]

        for cellFolder in eligibleCellFolders:
            if cellFolder in allCellFoldersInZone:

                cellLat, cellLon = getCellLatAndLongFromFolderName(cellFolder)

                ix = np.argwhere(lats == cellLat).flatten()[0]
                iy = np.argwhere(lons == cellLon).flatten()[0]

                cellAvgWaterTCurrYear = [['date', 'waterT', 'flow']]

                for i, d in enumerate(date):
                    cellAvgWaterTCurrYear = cellAvgWaterTCurrYear + \
                                            [[date[i], waterT[i, ix, iy], streamflow[i, ix, iy]]]

                eligibleCellWaterTs[cellFolder] = cellAvgWaterTCurrYear

    else:
        for cellFolder in eligibleCellFolders:
            if cellFolder in allCellFoldersInZone:

                cellFolderAvgTFilename = cellFolder + 'Average.csv'
                avgTFilePath = os.path.join(curtailparam.rbmOutputDir, cellFolder, cellFolderAvgTFilename)
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
