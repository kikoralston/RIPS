"""
Michael Craig
August 11, 2016
Script for modifying generator capacity with water T from RBM.

1) Gets list of all grid cells from .spat.
2) Saves T data for each segment in each grid cell into folders for each cell.
3) Averages T data across all segments for each cell and saves that data in folder for each cell.
4) Creates dictionaries mapping generators to cells and vice versa.
5) Using dictionary, pairs generator with water T time series and modifies capacity.

August 2018
UW changed the format of the files with water data to NETCDF. I added a flag 'netcdf' to some functions to read them.
"""

import os, csv, copy
import numpy as np
import pandas as pd
import datetime as dt
import netCDF4 as nc
import pickle as pk
import time
from AuxFuncs import *
import multiprocessing as mp
from GAMSAuxFuncs import *
from AuxCurtailmentFuncs import *
from CurtailmentRegressions import (calcCurtailmentForGenOrTech, loadRegCoeffs, getKeyCurtailParams,
                                    getCoeffsForGenOrTech)
import progressbar


def determineHrlyCurtailmentsForExistingGens(genFleet, currYear, genparam, curtailparam):
    """Computes times series of hourly curtailments for EXISTING generators

    :param genFleet:
    :param currYear:
    :param genparam:
    :param curtailparam:
    :return: dictionary of ORIS+UNITID to 2d vertical list of [datetime,curtailment(mw)],
             list of (ORIS+UNIT,gen lat long, cell lat long)

    """
    (genToCellLatLongsDict, cellLatLongToGenDict,
     genToCellLatLongsList) = getGenToCellAndCellToGenDictionaries(genFleet)

    # compile list of arguments to call multiprocessing
    args_list = [[cellLatLongToGenDict, currYear, genFleet, genparam, curtailparam, gcm]
                 for gcm in curtailparam.listgcms]

    ncores = min(len(curtailparam.listgcms), genparam.ncores_py)

    if ncores == 1:
        # do normal map instead of pool.map to avoid overhead
        list_curtailments = list(map(worker_generator_curtailments, args_list))
    else:
        # if ncores > 1, use multiprocessing
        with mp.Pool(processes=ncores) as pool:
            list_curtailments = pool.map(worker_generator_curtailments, args_list)

    # create nested dictionary
    hrlyCurtailmentsAllGensInTgtYr = dict(zip(curtailparam.listgcms, list_curtailments))

    return hrlyCurtailmentsAllGensInTgtYr, genToCellLatLongsList


def worker_generator_curtailments(list_args):
    """ Worker function for using multiprocessing with function calculateGeneratorCurtailments

    :param list_args: list with arguments for function
    :return:
    """

    curtailments_out = calculateGeneratorCurtailments(cellLatLongToGenDict=list_args[0], currYear=list_args[1],
                                                      genFleet=list_args[2], genparam=list_args[3],
                                                      curtailparam=list_args[4], gcm=list_args[5])

    return curtailments_out


def getGenToCellAndCellToGenDictionaries(genFleet):
    """MAP GENERATORS TO AND FROM RBM GRID CELLS

    Takes the 2d list with generator fleet data and creates dictionaries mapping individual plants to grid cells
    lat and long

    :param genFleet: 2d list with generator fleet data
    :return: genToCellLatLongsDict: dictionary gen:[gen,(gen lat, gen long), (cell lat, cell long)]
             cellLatLongToGenDict:  dictionary (cell lat, cell long):genID.
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


def calculateGeneratorCurtailments(cellLatLongToGenDict, currYear, genFleet, genparam, curtailparam, gcm,
                                   netcdf=True, pbar=True):

    """CURTAIL GENERATOR CAPACITY WITH WATER TEMPERATURES

    Calculates generator curtailments. If generator isn't curtailed (not right plant type, cell data not avaialble,
    etc.), ignored here. CalculateHourlyCapacs... script handles those cases (assume curtailment = 0).
    cellLatLongToGenDict = dict of (cell lat, cell long):gen ID for all gens w/ lat/lon coords in genFleet,
    including gens that should nto be curtailed.

    :param cellLatLongToGenDict:
    :param currYear:
    :param genFleet:
    :param genparam:
    :param curtailparam:
    :param netcdf: True if data comes in netcdf format
    :param pbar: True to show progress bar
    :return: hrlyCurtailmentsAllGensInTgtYr: nested dictionary with {gcm:{orisId: [1d np array with hourly curtailments]}}
             hrlyCurtailmentsList: 2d list [[orisId] + [hourly curtailments]]
    """

    hrlyCurtailmentsAllGensInTgtYr, hrlyCurtailmentsList = dict(), list()

    if netcdf:
        allCellFolders = get_all_cells_from_netcdf(curtailparam, gcm)
    else:
        allCellFolders = os.listdir(curtailparam.rbmOutputDir)

    # dict of planttype:coolingtype:cooldesignT:param:coeffs
    regCoeffs = loadRegCoeffs(genparam.dataRoot, 'capacity.json')

    if pbar:
        ProgressBar = progressbar.ProgressBar
    else:
        ProgressBar = progressbar.NullBar

    bar = ProgressBar()

    allCellFoldersInZone = get_all_cells_in_zone(allCellFolders, genparam, curtailparam)

    cellWaterTsForGens = loadCellWaterTs(allCellFolders, allCellFoldersInZone, curtailparam, gcm, currYear)

    fname_meteo = curtailparam.basenamemeteo
    fname_meteo = os.path.join(curtailparam.rbmDataDir, fname_meteo.format(gcm, currYear))
    meteodata = read_netcdf_full(currYear, fname_meteo, curtailparam)

    # this maps gen lat/lon to gen IDs; cell data may not exist
    for (cellLat, cellLong) in bar(cellLatLongToGenDict):

        cellFoldername = createBaseFilenameToReadOrWrite(curtailparam.locPrecision, cellLat, cellLong)

        if cellFoldername in allCellFoldersInZone:

            metAndWaterData = loadWaterAndMetData(currYear, cellLat, cellLong, genparam, curtailparam,
                                                  metdatatot=meteodata, waterDatatot=cellWaterTsForGens)

            gensInCell = cellLatLongToGenDict[(cellLat, cellLong)]  # list of ORIS-UNITID in cell

            for gen in gensInCell:

                (plantType, hr, fuelAndCoalType, coolType,
                 fgdType, state, capac, coolDesignT) = getKeyCurtailParams(gen, genFleet)

                coeffs = getCoeffsForGenOrTech(plantType, coolType, genparam.ptCurtailed, regCoeffs, coolDesignT)

                if (coeffs is not None) and (plantType in genparam.ptCurtailedRegs):  # skip gens that aren't curtailed
                    hrlyCurtailmentsGen = calcCurtailmentForGenOrTech(plantType, fuelAndCoalType, coolType, state,
                                                                      capac, coolDesignT, metAndWaterData, coeffs,
                                                                      genparam, curtailparam)

                    # add result to dictionary
                    hrlyCurtailmentsAllGensInTgtYr[gen] = hrlyCurtailmentsGen

                    # concatenate lists and add to the master list
                    hrlyCurtailmentsList.append([gen] + list(hrlyCurtailmentsGen))

        else:
            print('Cell not in folders!:', cellFoldername)

    return hrlyCurtailmentsAllGensInTgtYr


def find125GridMaurerLatLong(lat, lon):
    """Get (lat,lon) of 1/8 grid cell that a (lat, lon) point falls in

    :param lat: latitude value
    :param lon:  longitude value
    :return: tuple (lat_grid, long_grid) of grid cell that contains the point (lat,lon)
    """

    lat_grid = np.around(8.0 * lat - 0.5) / 8.0 + 0.0625
    lon_grid = np.around(8.0 * lon - 0.5) / 8.0 + 0.0625
    return lat_grid, lon_grid


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


