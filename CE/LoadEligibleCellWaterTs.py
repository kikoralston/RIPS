""""
Michael Craig, 17 May 2017
Functions load cell water T data from cells eligible for new generators by first setting
the folder names with future data for cells (each folder = 1 cell), then loading water temperatures
from those folders.

August 2018
UW changed the format of the files with water data to NETCDF. I added a flag 'netcdf' to some functions to read them.
"""

import os, copy
import netCDF4 as nc
import numpy as np
import datetime as dt
import pandas as pd
import pickle as pk
from AuxFuncs import *
import pickle as pk
from AssignCellsToIPMZones import mapCellToIPMZone
from ModifyGeneratorCapacityWithWaterTData import (getGenToCellAndCellToGenDictionaries,
                                                   createBaseFilenameToReadOrWrite)


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
    :return:
    """
    allCellFoldersInZone, eligibleCellFolders = setCellFolders(genFleet, genparam, curtailparam)

    return loadCellWaterTs(eligibleCellFolders, allCellFoldersInZone, curtailparam, currYear)


# rbmOutputDir, cellsEligibleForNewPlants, genFleet, locPrecision, zones, fipsToZones, fipsToPolys

# rbmOutputDir, cellsEligibleForNewPlants, genFleet, locPrecision, zones, fipsToZones, fipsToPolys, netcdf=True

def setCellFolders(genFleet, genparam, curtailparam, netcdf=True):
    """Get list of folders with future data for 1 cell each

    :param genFleet: 2d list with generator fleet
    :param netcdf: (boolean) if true reads water data from NETCDF files (new format from August 2018)
    :return:
    """

    if netcdf:
        # TODO: fname is the name of the netcdf file with water/meteo data. need to decide how to define it
        dataset1 = nc.Dataset(os.path.join(curtailparam.rbmDataDir, fname))

        lats = dataset1.variables['lat'][:]
        lons = dataset1.variables['lon'][:]

        allCellFolders = ['{}_{}'.format(la, lo) for la in lats for lo in lons]

        # read file with dictionary mapping cells to IPM zones (this file must be in the VICS/RBM folder)
        with open(os.path.join(curtailparam.rbmRootDir, 'cells2zones.pk'), 'rb') as f:
            cells2zones = pk.load(f)

        allCellFoldersInZone = [c for c in allCellFolders if cells2zones[c] in zones]
    else:
        allCellFolders = [name for name in os.listdir(curtailparam.rbmOutputDir)
                          if os.path.isdir(os.path.join(curtailparam.rbmOutputDir, name))]

        allCellFoldersInZone = isolateCellsInZones(allCellFolders, genparam)

    if genparam.cellsEligibleForNewPlants == 'all':
        eligibleCellFolders = copy.deepcopy(allCellFoldersInZone)

    elif genparam.cellsEligibleForNewPlants == 'withGens':
        # (genToCellLatLongsDictValues,genToCellLatLongsDict,cellLatLongToGenDict,
        #             genToCellLatLongsList) = getGenToCellAndCellToGenDictionaries(genFleet)
        (genToCellLatLongsDict,cellLatLongToGenDict,
                    genToCellLatLongsList) = getGenToCellAndCellToGenDictionaries(genFleet)
        eligibleCellFolders = []
        for (cellLat, cellLong) in cellLatLongToGenDict:
            eligibleCellFolders.append(createBaseFilenameToReadOrWrite(curtailparam.locPrecision, cellLat, cellLong))

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
        waterT, date, lons, lats = read_waterdata_netcdf(curtailparam.rbmDataDir, 'serc.NorESM1-M.RCP85.stream_T.nc', 'T_stream', currYear)
        streamflow = read_waterdata_netcdf(curtailparam.rbmDataDir, 'serc.NorESM1-M.RCP85.KW.flow.MACA.regulated.nc', 'streamflow', currYear)[0]

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


def read_waterdata_netcdf(path_data, fname, varname, curryear):
    """ Reads netcdf file with water data and returns arrays with data, lats and longs

    :param path_data:
    :param fname:
    :param varname:
    :param curryear:
    :return:
            data_values: array with data values
            date: 1d array with dates
            lons: 1d array with longitude values
            lats: 1d array with latitude values
    """
    dataset1 = nc.Dataset(os.path.join(path_data, fname))

    lats = dataset1.variables['lat'][:]
    lons = dataset1.variables['lon'][:]

    # Variables to read:
    #       - 'T_stream'
    #       - 'streamflow'

    # get reference date of data and number of days of data
    strdate = dataset1.variables['time'].units.replace('days since', '').strip()

    ref_date = dt.datetime.strptime(strdate, '%Y-%m-%d %H:%M:%S').date()
    n_days = dataset1.variables['time'].shape[0]

    # create array with dates
    date_array = pd.date_range(start=ref_date, periods=n_days, freq='D')

    if curryear is not None:
        # get indexes of days for year == 'curryear'

        year = curryear
        idx_year = np.where(date_array.year == year)[0]

        # slice data for this year
        if varname == 'T_stream':
            data_values = dataset1.variables[varname][idx_year, :, :, :]
        else:
            data_values = dataset1.variables[varname][idx_year, :, :]
    else:
        data_values = dataset1.variables[varname][:]

    date_array = date_array[idx_year]

    if varname == 'T_stream':
        # for T_stream average data of multiple segments for each day
        pos_no_seg = dataset1.variables['T_stream'].dimensions.index('no_seg')
        data_values = np.mean(data_values, axis=pos_no_seg)

    return data_values, date_array, lons, lats


def read_waterdata_cell(path_data, currYear, cellLat, cellLon):
    """ Reads netcdf file with water data and returns panda data frame for data for given cell

    :param path_data:
    :param fname:
    :param varname:
    :param curryear:
    :return:
            data_values: array with data values
            date: 1d array with dates
            lons: 1d array with longitude values
            lats: 1d array with latitude values
    """
    waterT, date, lons, lats = read_waterdata_netcdf(path_data, 'serc.NorESM1-M.RCP85.stream_T.nc', 'T_stream',
                                                     currYear)
    streamflow = read_waterdata_netcdf(path_data, 'serc.NorESM1-M.RCP85.KW.flow.MACA.regulated.nc', 'streamflow',
                                       currYear)[0]

    ix = np.argwhere(lats == cellLat).flatten()[0]
    iy = np.argwhere(lons == cellLon).flatten()[0]

    # get data for cell
    waterT = waterT.data[:, ix, iy]
    streamflow = streamflow.data[:, ix, iy]

    df_water = pd.DataFrame(data={'date': date, 'waterT': waterT, 'flow': streamflow})

    return df_water

