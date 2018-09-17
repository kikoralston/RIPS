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
from GAMSAuxFuncs import *
from CurtailmentRegressions import (calcCurtailmentForGenOrTech, loadRegCoeffs, getKeyCurtailParams,
                                    getCoeffsForGenOrTech)
from PreProcessRBM import createBaseFilenameToReadOrWrite
import progressbar


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
        calculateGeneratorCurtailments(cellLatLongToGenDict, currYear, genFleet, modelName, genparam, curtailparam,
                                       resultsDir)

    return hrlyCurtailmentsAllGensInTgtYr, genToCellLatLongsList, hrlyCurtailmentsList


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


def loadWaterAndMetData(curtailmentYear, cellLat, cellLong, genparam, curtailparam, metdatatot=None, waterDatatot=None,
                        netcdf=True):
    """LOAD WATER AND MET DATA BY CELL ON HOURLY BASIS

    Read text files with data for met (air T and rh) and water T data for one specific cell. The data files were
    created by the pre-processing function 'processRBMDataIntoIndividualCellFiles' and saved in a folder 'rbmOutputDir'.

    :param curtailmentYear: integer with year that curtailment is being simulated
    :param cellLat:
    :param cellLong:
    :param genparam:
    :param curtailparam:
    :param metdatatot:
    :param waterDatatot:
    :param netcdf: (boolean) if true reads water data from NETCDF files (new format from August 2018)
    :return: panda data frame
    """

    cellFoldername = createBaseFilenameToReadOrWrite(curtailparam.locPrecision, cellLat, cellLong)

    if metdatatot is None:
        # if no object with pre loaded meteo data is given, read netcdf file

        # Load netcdf files w/ air T, rh, air pressure and datetime just for current year
        if genparam.analysisArea == 'test':
            # for now read with dara from some random
            metData = read_netcdf(34.3125, -86.6875, curtailmentYear, curtailparam)
        else:
            metData = read_netcdf(cellLat, cellLong, curtailmentYear, curtailparam)
    else:
        # get cell data from pre loaded meteo data
        metData = read_dict_necdf(cellLat, cellLong, metdatatot)

    if waterDatatot is None:
        if netcdf:
            waterData = read_waterdata_cell(curtailparam, curtailmentYear, cellLat, cellLong)

        else:
            # load csv files with water temp and stream flows and filter current year
            cellFolder = os.path.join(curtailparam.rbmOutputDir, cellFoldername)
            averageTFilename = cellFoldername + 'Average.csv'
            waterData = pd.read_csv(filepath_or_buffer=os.path.join(cellFolder, averageTFilename))
            waterData = waterData[['date', 'waterT', 'flow']]

            # filter current year
            y = [int(date[0:4]) == int(curtailmentYear) for date in waterData['date']]
            waterData = waterData[y]
            waterData = waterData.reset_index(drop=True)
    else:
        waterData = copy.deepcopy(waterDatatot[cellFoldername])
        # convert from 2d list to pandas DF
        headers = waterData.pop(0)
        waterData = pd.DataFrame(waterData, columns=headers)

    # expand water data to hourly
    waterData = expand_df_hourly(waterData)

    # merge water and meteo data for this cell
    totalData = waterData.merge(metData, on='date')

    # create Farenheit
    totalData['waterF'] = totalData.loc[:, 'waterT'] * 9 / 5 + 32
    totalData['airF'] = totalData.loc[:, 'airT'] * 9 / 5 + 32

    return totalData


def calculateGeneratorCurtailments(cellLatLongToGenDict, curtailmentYear, genFleet, modelName, genparam,
                                   curtailparam, resultsDir, netcdf=True):
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
    :param netcdf:
    :return:
    """

    if netcdf:
        allCellFolders = get_all_cells_from_netcdf(curtailparam)
    else:
        allCellFolders = os.listdir(curtailparam.rbmOutputDir)

    hrlyCurtailmentsAllGensInTgtYr, hrlyCurtailmentsList = dict(), list()

    # dict of planttype:coolingtype:cooldesignT:param:coeffs
    regCoeffs = loadRegCoeffs(genparam.dataRoot, 'capacity.json')

    # this maps gen lat/lon to gen IDs; cell data may not exist
    for (cellLat, cellLong) in progressbar.progressbar(cellLatLongToGenDict):
    # for (cellLat, cellLong) in cellLatLongToGenDict:
        #print('cell: {}, {}'.format(cellLat, cellLong))

        cellFoldername = createBaseFilenameToReadOrWrite(curtailparam.locPrecision, cellLat, cellLong)

        if cellFoldername in allCellFolders:
            metAndWaterData = loadWaterAndMetData(curtailmentYear, cellLat, cellLong, genparam, curtailparam,
                                                  netcdf=netcdf)
            gensInCell = cellLatLongToGenDict[(cellLat, cellLong)]  # list of ORIS-UNITID in cell

            for gen in gensInCell:
                #print('genId: {}'.format(gen))
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

    return hrlyCurtailmentsAllGensInTgtYr, hrlyCurtailmentsList


def getCellLatAndLongFromFolderName(dummyFolder):
    """

    :param dummyFolder:
    :return:
    """
    cellLat, cellLon = dummyFolder.split('_')
    return float(cellLat), float(cellLon)


def expand_df_hourly(df):
    """utility function to expand a daily data frame to hourly

    :param df: daily data frame (check names of columns)
    :return: data frame
    """
    a = pd.to_datetime(df['date'], format='%Y-%m-%d')
    start_day, end_day = min(a), max(a)
    end_day = end_day + dt.timedelta(hours=23)

    date2 = pd.date_range(start_day, end_day, freq='h')
    daystr = date2.strftime('%Y-%m-%d')

    df2 = pd.DataFrame(data={'date_hour': date2, 'date': daystr})
    df['date'] = df['date'].dt.strftime('%Y-%m-%d')

    df2 = df2.merge(df, on='date', how='left')

    del df2['date']
    df2 = df2.rename(columns={'date_hour': 'date'})

    return df2


def read_netcdf_full(currYear, fname, curtailparam):
    """Reads NETCDF file with meteo data for current year

    :param currYear: (integer) current year
    :param fname: (string) full path to netcdf file
    :param curtailparam: object of class Curtailmentparameters
    :return: dictionary with arrays of data {'var name': array}
    """

    if os.path.isfile(fname):

        dataset = nc.Dataset(fname)

        # Extract data from NetCDF file
        lats = dataset.variables['lat'][:]
        lons = dataset.variables['lon'][:]
        temp = dataset.variables['temp'][:]
        air_pressure = dataset.variables['air_pressure'][:]
        rel_humid = dataset.variables['rel_humid'][:]

        # create array with hourly time stamps for year
        start = dt.datetime(currYear, 1, 1)
        end = dt.datetime(currYear, 12, 31, 23, 00, 00)
        date_array = pd.date_range(start=start, end=end, freq='H')

        dict_out = {'date': date_array, 'lat': lats, 'lon': lons, 'airT': temp, 'P': air_pressure, 'rh': rel_humid}
    else:
        print('{} not found!'.format(fname))
        dict_out = None

    return dict_out


def read_netcdf(cellLat, cellLon, currYear, curtailparam):
    """Reads NETCDF file with meteo data for current year and gets data for specified cell

    :param cellLat: (numeric) latitude of cell
    :param cellLon:  (numeric) longitude of cell
    :param currYear: (integer) current year
    :param curtailparam: object of class Curtailmentparameters
    :return: pandas data frame with hourly meteo data
    """
    fname = curtailparam.basenamemeteo

    # TODO: find a way to combine different GCMs
    gcm = curtailparam.listgcms[0]

    if os.path.isfile(os.path.join(curtailparam.rbmDataDir, fname.format(gcm, currYear))):

        fname = os.path.join(curtailparam.rbmDataDir, fname.format(gcm, currYear))
        dict_out = read_netcdf_full(currYear, fname, curtailparam)

        lats = dict_out['lat']
        lons = dict_out['lon']
        temp = dict_out['airT']
        air_pressure = dict_out['P']
        rel_humid = dict_out['rh']

        ix = np.argwhere(lats == cellLat).flatten()[0]
        iy = np.argwhere(lons == cellLon).flatten()[0]

        # create array with hourly time stamps for year
        start = dt.datetime(currYear, 1, 1)
        end = dt.datetime(currYear, 12, 31, 23, 00, 00)
        date_array = pd.date_range(start=start, end=end, freq='H')

        df_out = pd.DataFrame({'date': date_array, 'airT': temp[:, ix, iy], 'P': air_pressure[:, ix, iy],
                               'rh': rel_humid[:, ix, iy]},
                              columns=['date', 'airT', 'rh', 'P'])
    else:
        print('No METEO data (NETCDF files) on folder {0:} for year {1:4d}!'.format(curtailparam.rbmDataDir, currYear))
        df_out = None

    return df_out


def read_dict_necdf(cellLat, cellLon, meteoData):
    """Get meteorological data for cell stored in a dictionary that loaded this data from the netcdf file

    see function 'read_netcdf_full'

    :param cellLat: (numeric) latitude of cell
    :param cellLon: (numeric) longitude of cell
    :param meteoData: (dictionary) dictionary with complete meteo data from the netcdf file
    :return: pandas data frame with hourly meteo data
    """

    lats = meteoData['lat']
    lons = meteoData['lon']
    temp = meteoData['airT']
    air_pressure = meteoData['P']
    rel_humid = meteoData['rh']
    date_array = meteoData['date']

    ix = np.argwhere(lats == cellLat).flatten()[0]
    iy = np.argwhere(lons == cellLon).flatten()[0]

    df_out = pd.DataFrame({'date': date_array, 'airT': temp[:, ix, iy], 'P': air_pressure[:, ix, iy],
                           'rh': rel_humid[:, ix, iy]},
                          columns=['date', 'airT', 'rh', 'P'])

    return df_out


def convert_2dList_netcdf(listcurtail, curtailparam, fnameuw, fnameout='~/test.nc', ):
    """
    This function converts a 2d list with new generator curtailments to a netcdf file


    :param listcurtail: 2d list with curtailment data
    :param curtailparam: object of class curtailmentparameters
    :param fnameuw: (string) name of netcdf file with original meteo data from UW (no path)
    :param fnameout: (string) name of resulting netcdf file that will be created (full path)
    :return: nothing
    """

    prec = curtailparam.locPrecision

    # read netcdf file with meteo data form UW to get spatial limits
    dataset = nc.Dataset(os.path.join(curtailparam.rbmDataDir, fnameuw))
    time = dataset.variables['time'][:]
    lats = dataset.variables['lat'][:]
    lons = dataset.variables['lon'][:]
    temp = dataset.variables['temp'][:]

    # initialize new netcdf file

    datasetout = nc.Dataset(os.path.expanduser(fnameout), 'w')

    latdim = datasetout.createDimension('lat', len(lats))
    londim = datasetout.createDimension('lon', len(lons))
    timedim = datasetout.createDimension('time', len(time))

    # Create coordinate variables for 4-dimensions
    timevar = datasetout.createVariable('time', np.int32, ('time',))
    latsvar = datasetout.createVariable('latitude', np.float32, ('lat',))
    lonsvar = datasetout.createVariable('longitude', np.float32, ('lon',))

    timevar[:] = time
    latsvar[:] = lats
    lonsvar[:] = lons

    # get set of curtailed techs
    curtailed_techs_all = [l[0][0] for l in listcurtail]
    curtailed_techs = list(set(curtailed_techs_all))

    # loop through curtailed techs and extract all spatial results for each
    for tech in curtailed_techs:

        print('Saving spatial data of plant {} to netcdf file...'.format(tech))

        idx = [i for i, y in enumerate(curtailed_techs_all) if y == tech]

        curt_data_tech = [listcurtail[i] for i in idx]
        list_keys = [row[0] for row in curt_data_tech]

        arr_data = np.zeros(temp.shape)

        for la in lats:
            for lo in lons:

                ix = np.argwhere(lats == la).flatten()[0]
                iy = np.argwhere(lons == lo).flatten()[0]

                # tuple
                label_tech = (tech, "{0:.{2}f}_{1:.{2}f}".format(la, lo, prec))

                if label_tech in list_keys:
                    idx_tech = list_keys.index(label_tech)

                    # convert to number and to numpy array and store into 3d array
                    z = np.array([float(zz) for zz in curt_data_tech[idx_tech][1:]])
                    arr_data[:, ix, iy] = z

        # convert to masked array
        arr_data = np.ma.array(arr_data, mask=temp.mask)

        data_netcdf = datasetout.createVariable(tech, np.float32, ('time', 'lat', 'lon'))

        data_netcdf[:] = arr_data

    datasetout.close()


def read_bulk_water_meteo(currYear, curtailparam):

    t0 = time.time()
    basenamemeteo = curtailparam.basenamemeteo

    dataset = nc.Dataset(os.path.join(curtailparam.rbmDataDir, basenamemeteo.format(currYear)))

    # Extract data from NetCDF file
    lats = dataset.variables['lat'][:]
    lons = dataset.variables['lon'][:]
    temp = dataset.variables['temp'][:]
    air_pressure = dataset.variables['air_pressure'][:]
    rel_humid = dataset.variables['rel_humid'][:]

    # create array with hourly time stamps for year
    start = dt.datetime(currYear, 1, 1)
    end = dt.datetime(currYear, 12, 31, 23, 00, 00)
    date_array = pd.date_range(start=start, end=end, freq='H')

    waterT = np.zeros(temp.shape)
    flow = np.zeros(temp.shape)

    i = 1
    for la in lats:
        for lo in lons:
            cellFoldername = createBaseFilenameToReadOrWrite(curtailparam.locPrecision, la, lo)

            ix = np.argwhere(lats == la).flatten()[0]
            iy = np.argwhere(lons == lo).flatten()[0]

            cellFolder = os.path.join(curtailparam.rbmOutputDir, cellFoldername)
            averageTFilename = cellFoldername + 'Average.csv'

            if os.path.isfile(os.path.join(cellFolder, averageTFilename)):

                print('{0:d}. Reading {1}'.format(i, cellFoldername))

                waterData = pd.read_csv(filepath_or_buffer=os.path.join(cellFolder, averageTFilename))
                waterData = waterData[['date', 'waterT', 'flow']]

                # filter current year
                y = [int(date[0:4]) == int(currYear) for date in waterData['date']]
                waterData = waterData[y]
                waterData = waterData.reset_index(drop=True)

                waterData = expand_df_hourly(waterData)

                waterT[:, ix, iy] = np.array(waterData['waterT'])
                flow[:, ix, iy] = np.array(waterData['flow'])

                i = i + 1

    # convert to masked arrays
    waterT = np.ma.array(waterT, mask=temp.mask)
    flow = np.ma.array(flow, mask=temp.mask)

    print(str_elapsedtime(t0))


def read_waterdata_netcdf(fname, varname, curryear):
    """ Reads netcdf file with water data and returns arrays with data, lats and longs

    :param fname:
    :param varname:
    :param curryear:
    :return:
            data_values: array with data values
            date: 1d array with dates
            lons: 1d array with longitude values
            lats: 1d array with latitude values
    """
    dataset1 = nc.Dataset(fname)

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


def read_waterdata_cell(curtailparam, currYear, cellLat, cellLon):
    """ Reads netcdf file with water data and returns panda data frame for data for given cell

    :param curtailparam
    :param currYear:
    :return:
            data_values: array with data values
            date: 1d array with dates
            lons: 1d array with longitude values
            lats: 1d array with latitude values
    """

    # TODO: find a way to combine different GCMs
    gcm = curtailparam.listgcms[0]

    # reading flow/streamT file. substitute '_' to '.'
    gcm = gcm.replace('_', '.')

    fname = curtailparam.basenamestreamT
    fname = os.path.join(curtailparam.rbmDataDir, fname.format(gcm))
    waterT, date, lons, lats = read_waterdata_netcdf(fname, 'T_stream', currYear)

    fname = curtailparam.basenameflow
    fname = os.path.join(curtailparam.rbmDataDir, fname.format(gcm))
    streamflow = read_waterdata_netcdf(fname, 'streamflow', currYear)[0]

    ix = np.argwhere(lats == cellLat).flatten()[0]
    iy = np.argwhere(lons == cellLon).flatten()[0]

    # get data for cell
    waterT = waterT.data[:, ix, iy]
    streamflow = streamflow.data[:, ix, iy]

    df_water = pd.DataFrame(data={'date': date, 'waterT': waterT, 'flow': streamflow})

    return df_water


def get_all_cells_from_netcdf(curtailparam):

    fname = curtailparam.basenameflow
    gcm = curtailparam.listgcms[0]

    # reading flow file. substitute '_' to '.'
    gcm = gcm.replace('_', '.')

    if os.path.isfile(os.path.join(curtailparam.rbmDataDir, fname.format(gcm))):

        dataset1 = nc.Dataset(os.path.join(curtailparam.rbmDataDir, fname.format(gcm)))

        lats = dataset1.variables['lat'][:]
        lons = dataset1.variables['lon'][:]
        mask = dataset1.variables['streamflow'][0, :, :].mask

        # compile list of cell grids. Ignore those that are masked.
        allCellFolders = ['{}_{}'.format(la, lo) for (i, la) in enumerate(lats)
                          for (j, lo) in enumerate(lons) if not mask[i, j]]

    return allCellFolders


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


def check_water_data(genparam, curtailparam, currYear):

    allCellFolders = get_all_cells_from_netcdf(curtailparam)

    for cellname in allCellFolders:

        cellLat, cellLong = tuple(map(float, cellname.split('_')))

        metAndWaterData = loadWaterAndMetData(currYear, cellLat, cellLong, genparam, curtailparam, netcdf=netcdf)

        print('cell: {}_{}'.format(cellLat, cellLong))
        print('  range water T: {} - {}'.format(metAndWaterData['waterT'].min(), metAndWaterData['waterT'].max()))
        print('  # na: {}'.format(sum(metAndWaterData['waterT'].isna())))
        print('  range flow: {} - {}'.format(metAndWaterData['flow'].min(), metAndWaterData['flow'].max()))
        print('  # na: {}'.format(sum(metAndWaterData['flow'].isna())))
        print()


def order_cells_by_flow(genparam, curtailparam, currYear, n=100):

    allCellFolders = get_all_cells_from_netcdf(curtailparam)

    # TODO: find a way to combine different GCMs
    gcm = curtailparam.listgcms[0]

    # reading flow/streamT file. substitute '_' to '.'
    gcm = gcm.replace('_', '.')

    fname = curtailparam.basenameflow
    fname = os.path.join(curtailparam.rbmDataDir, fname.format(gcm))
    streamflow, date, lons, lats = read_waterdata_netcdf(fname, 'streamflow', currYear)

    avg_stream_flow = np.mean(streamflow.data[:, :, :], axis=0)

    # read file with dictionary mapping cells to IPM zones (this file must be in the VICS/RBM folder)
    with open(os.path.join(curtailparam.rbmRootDir, 'cells2zones.pk'), 'rb') as f:
        cells2zones = pk.load(f)

    best_cells_by_zone = {}

    for z in genparam.ipmZones:
        cellFoldersInZone = [c for c in allCellFolders if cells2zones[c] == z]

        avg_stream_flow_in_zone = []

        for c in cellFoldersInZone:
            cellLat, cellLon = tuple(map(float, c.split('_')))

            ix = np.argwhere(lats == cellLat).flatten()[0]
            iy = np.argwhere(lons == cellLon).flatten()[0]

            avg_stream_flow_in_zone = avg_stream_flow_in_zone + [avg_stream_flow[ix, iy]]

        # add to pandas data frame to sort

        df = pd.DataFrame({'cell': cellFoldersInZone, 'flow': avg_stream_flow_in_zone})

        df = df.sort_values(by=['flow'], ascending=False)

        df = df.head(n)

        best_cells = list(df['cell'])

        best_cells_by_zone[z] = best_cells

    return best_cells_by_zone




