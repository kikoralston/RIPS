# October 2018
# This module has auxiliary functions for reading and processing meteorological and water data from the different GCMs
# These functions are used by modules ModifyGeneratorCapacityWithWaterTDAta.py, ModifyNewTechCapacityWithWaterTDAta.py, LoadEligibleCellWaterTs.py
# Previously, most of these functions were spread over these modules, but these was creating some circular references.


import os
import copy
import netCDF4 as nc
import pandas as pd
import numpy as np
import pickle as pk
import datetime as dt


def loadWaterAndMetData(curtailmentYear, cellLat, cellLong, genparam, curtailparam, gcm, metdatatot=None,
                        waterDatatot=None, netcdf=True):
    """Load water and met data of a single cell (and a single gcm) on hourly basis for given year

    This function can take dictionaries with previously loaded full gridded datasets of weather (`metdatatot`) and
    water data (`waterDatatot`) and it filters the data of a single grid cell. If these dictionaries have not been
    loaded previously, the function reads from the netcdf files.

    :param curtailmentYear: (int) year of simulation
    :param cellLat: (float) latitude of grid cell
    :param cellLong: (float) longitude of grid cell
    :param genparam: object of type :mod:`Generalparameters`
    :param curtailparam: object of type :mod:`Curtailmentparameters`
    :param gcm: (string) GCM name
    :param metdatatot: dict with pre-loaded complete meteo data from the netcdf file (see :func:`.read_netcdf_full`).
                       If this is `None` the function loads the data directly from the netcdf file (:func:`.read_netcdf`)
    :param waterDatatot: dict of {cell folder name : [[Datetime],[AverageWaterT(degC)], [AirT], [flow]]}. If value is
                         `None`, it loads data from netcdf file (see :func:`.read_waterdata_cell`)
    :param netcdf: (boolean) if true reads water data from NETCDF files (Old format does not work anymore)
    :return: panda data frame with hourly water and weather data for given year and given cell
    """

    cellFoldername = createBaseFilenameToReadOrWrite(curtailparam.locPrecision, cellLat, cellLong)

    if metdatatot is None:
        # if no object with pre loaded meteo data is given, read netcdf file
        metData = read_netcdf(cellLat, cellLong, curtailmentYear, curtailparam, gcm)
    else:
        # get cell data from pre loaded meteo data
        metData = read_dict_netcdf(cellLat, cellLong, metdatatot)

    if waterDatatot is None:
        if netcdf:
            waterData = read_waterdata_cell(curtailparam=curtailparam, currYear=curtailmentYear,  cellLat=cellLat,
                                            cellLon=cellLong, gcm=gcm)
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


def loadCellWaterTs(eligibleCellFolders, allCellFoldersInZone, curtailparam, gcm, currYear=None, netcdf=True):
    """Load water and weather data for all eligible grid cells in study

    Returns dict of {cell folder name : [[Datetime, AverageWaterT(degC), flow]]}. Each element in the dictionary is a
    2d list.

    :param eligibleCellFolders: (list) list with names of eligible cells for study
    :param allCellFoldersInZone: (list) list with names of cells inside zones included in study
    :param curtailparam: object of type :mod:`Curtailmentparameters`
    :param gcm: (string) name of gcm being imported
    :param currYear: (integer) year of data being imported. If None, imports all years
    :param netcdf: (boolean) if true reads water data from NETCDF files (new format from August 2018)

    :return: dict of {cell folder name : [[Datetime, AverageWaterT(degC), flow]]}
    """
    eligibleCellWaterTs = dict()

    if netcdf:
        name_gcm = gcm

        # reading flow/streamT file. substitute '_' to '.'
        name_gcm = name_gcm.replace('_', '.')
        name_gcm = name_gcm.replace('rcp', 'RCP')

        fname = curtailparam.basenamestreamT
        fname = os.path.join(curtailparam.rbmRootDir, fname.format(name_gcm))
        waterT, date, lons, lats = read_waterdata_netcdf(fname, 'T_stream', currYear)

        fname = curtailparam.basenameflow
        fname = os.path.join(curtailparam.rbmRootDir, fname.format(name_gcm))
        streamflow = read_waterdata_netcdf(fname, 'streamflow', currYear)[0]

        for cellFolder in eligibleCellFolders:
            if cellFolder in allCellFoldersInZone:

                if waterT is not None and streamflow is not None:

                    cellLat, cellLon = getCellLatAndLongFromFolderName(cellFolder)

                    ix = np.argwhere(lats == cellLat).flatten()[0]
                    iy = np.argwhere(lons == cellLon).flatten()[0]

                    cellAvgWaterTCurrYear = [['date', 'waterT', 'flow']]

                    for i, d in enumerate(date):
                        cellAvgWaterTCurrYear = cellAvgWaterTCurrYear + \
                                                [[date[i], waterT[i, ix, iy], streamflow[i, ix, iy]]]

                    eligibleCellWaterTs[cellFolder] = cellAvgWaterTCurrYear

                else:
                    eligibleCellWaterTs[cellFolder] = None

    else:
        # old format
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


def read_netcdf_full(currYear, fname, curtailparam):
    """Reads NETCDF file with meteo data for current year

    This function reads the netcdf file and returns a dictionary with 1-d arrays of the dimension variables
    (time, lat, lon) and 3-d arrays of weather and water data (air temp, air pressure, relative humidity).

    It reads data for all cells. Netcdf files have data for a single year.

    :param currYear: (integer) current year
    :param fname: (string) full path to netcdf file
    :param curtailparam: object of type :mod:`Curtailmentparameters`
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
        dataset.close()
    else:
        print('{} not found!'.format(fname))
        dict_out = None

    return dict_out


def read_netcdf(cellLat, cellLon, currYear, curtailparam, gcm):
    """Reads NETCDF file with meteo data for current year and gets data for specified cell

    :param cellLat: (numeric) latitude of cell
    :param cellLon:  (numeric) longitude of cell
    :param currYear: (integer) current year
    :param curtailparam: object of type :mod:`Curtailmentparameters`
    :return: pandas data frame with hourly meteo data
    """
    fname = curtailparam.basenamemeteo

    if os.path.isfile(os.path.join(curtailparam.rbmRootDir, fname.format(gcm, currYear))):

        fname = os.path.join(curtailparam.rbmRootDir, fname.format(gcm, currYear))

        dataset = nc.Dataset(fname)

        # Extract lat and lon from NetCDF file
        lats = dataset.variables['lat'][:]
        lons = dataset.variables['lon'][:]

        # get indexes of position of cell
        ix = np.argwhere(lats == cellLat).flatten()[0]
        iy = np.argwhere(lons == cellLon).flatten()[0]

        # get meteodata for cell
        temp = dataset.variables['temp'][:, ix, iy]
        air_pressure = dataset.variables['air_pressure'][:, ix, iy]
        rel_humid = dataset.variables['rel_humid'][:, ix, iy]

        dataset.close()

        # create array with hourly time stamps for year
        start = dt.datetime(currYear, 1, 1)
        end = dt.datetime(currYear, 12, 31, 23, 00, 00)
        date_array = pd.date_range(start=start, end=end, freq='H')

        # create data frame
        df_out = pd.DataFrame({'date': date_array, 'airT': temp, 'P': air_pressure, 'rh': rel_humid},
                              columns=['date', 'airT', 'rh', 'P'])
    else:
        print('No METEO data (NETCDF files) on folder {0:} for year {1:4d}!'.format(curtailparam.rbmRootDir, currYear))
        df_out = None

    return df_out


def read_dict_netcdf(cellLat, cellLon, meteoData):
    """Get meteorological data for cell stored in a dictionary with full dataset from the netcdf file

    see function :func:`.read_netcdf_full`

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


def get_all_cells_from_netcdf(curtailparam, gcm):
    """Get all cells from netcdf file

    This function reads a netcdf file for the given gcm and returns a list with cells lat and longs as strings

    :param curtailparam: object of type :mod:`Curtailmentparameters`
    :param gcm: (string) name of gcm
    :return: list with cells lat and longs as strings
    """
    fname = curtailparam.basenameflow
    name_gcm = gcm

    # reading flow file. substitute '_' to '.' and change 'RCP' to small caps
    name_gcm = name_gcm.replace('_', '.')
    name_gcm = name_gcm.replace('rcp', 'RCP')

    allCellFolders = None

    if os.path.isfile(os.path.join(curtailparam.rbmRootDir, fname.format(name_gcm))):

        dataset1 = nc.Dataset(os.path.join(curtailparam.rbmRootDir, fname.format(name_gcm)))

        lats = dataset1.variables['lat'][:]
        lons = dataset1.variables['lon'][:]
        mask = dataset1.variables['streamflow'][0, :, :].mask

        # compile list of cell grids. Ignore those that are masked.
        allCellFolders = ['{}_{}'.format(la, lo) for (i, la) in enumerate(lats)
                          for (j, lo) in enumerate(lons) if not mask[i, j]]
        dataset1.close()
    else:
        print('File {} not found!'.format(os.path.join(curtailparam.rbmRootDir, fname.format(name_gcm))))

    return allCellFolders


def get_all_cells_in_zone(allCellFolders, genparam, curtailparam):
    """ Compiles list of cells that are inside the ipm zones

    :param allCellFolders: list with all cells
    :param curtailparam: object of type :mod:`Generalparameters`
    :param curtailparam: object of type :mod:`Curtailmentparameters`
    :return: allCellFoldersInZone: list with all cells that are inside the ipm zones
    """
    # read file with dictionary mapping cells to IPM zones (this file must be in the VICS/RBM folder)
    with open(os.path.join(curtailparam.rbmRootDir, 'cells2zones.pk'), 'rb') as f:
        cells2zones = pk.load(f)

    if allCellFolders is not None:
        allCellFoldersInZone = [c for c in allCellFolders if cells2zones[c] in genparam.ipmZones]
    else:
        allCellFoldersInZone = [c for c in cells2zones if cells2zones[c] in genparam.ipmZones]

    return allCellFoldersInZone


def read_waterdata_cell(curtailparam, currYear, cellLat, cellLon, gcm):
    """ Reads netcdf file with water data and returns panda data frame for data for given cell

    :param curtailparam: object of type :mod:`Curtailmentparameters`
    :param currYear: (int) current year
    :return: data frame with following columns
             * date: 1d array with dates
             * lons: 1d array with longitude values
             * lats: 1d array with latitude values
    """

    # reading flow/streamT file. substitute '_' to '.'
    gcm = gcm.replace('_', '.')
    gcm = gcm.replace('rcp', 'RCP')

    fname = curtailparam.basenamestreamT
    fname = os.path.join(curtailparam.rbmRootDir, fname.format(gcm))
    waterT, date, lons, lats = read_waterdata_netcdf(fname, 'T_stream', currYear)

    fname = curtailparam.basenameflow
    fname = os.path.join(curtailparam.rbmRootDir, fname.format(gcm))
    streamflow = read_waterdata_netcdf(fname, 'streamflow', currYear)[0]

    ix = np.argwhere(lats == cellLat).flatten()[0]
    iy = np.argwhere(lons == cellLon).flatten()[0]

    # get data for cell
    waterT = waterT.data[:, ix, iy]
    streamflow = streamflow.data[:, ix, iy]

    df_water = pd.DataFrame(data={'date': date, 'waterT': waterT, 'flow': streamflow})

    return df_water


def read_waterdata_netcdf(fname, varname, curryear):
    """ Reads netcdf file with water data and returns arrays with data, lats and longs for given year

    :param fname: (string) complete path fo netcdf file
    :param varname: (string) name of variable to load from netcdf file
    :param curryear: (int) current year
    :return: tuple with following elements
            * data_values: array with data values
            * date: 1d array with dates
            * lons: 1d array with longitude values
            * lats: 1d array with latitude values
    """

    if os.path.isfile(fname):
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
        idx_year = None

        if curryear is not None:
            # get indexes of days for year == 'curryear'

            year = curryear
            idx_year = np.where(date_array.year == year)[0]

            # slice data for this year
            if varname == 'T_stream':
                data_values = dataset1.variables[varname][idx_year, :, :, :]
            else:
                data_values = dataset1.variables[varname][idx_year, :, :]

            date_array = date_array[idx_year]
        else:
            data_values = dataset1.variables[varname][:]

        if varname == 'T_stream':
            # for T_stream average data of multiple segments for each day
            pos_no_seg = dataset1.variables['T_stream'].dimensions.index('no_seg')
            data_values = np.mean(data_values, axis=pos_no_seg)

        dataset1.close()
    else:
        print('File {} not found!'.format(fname))
        data_values, date_array, lons, lats = None, None, None, None

    return data_values, date_array, lons, lats


def order_cells_by_flow(genparam, curtailparam, currYear, n=100, output_list=True):
    """Order cells according to average annual water flow across all gcms

    This function order cells in each zone according to annual average river flow. It returns the `n` cells with largest
    flow as an ordered list.

    :param genparam: object of type :mod:`Generalparameters`
    :param curtailparam: object of type :mod:`Curtailmentparameters`
    :param currYear: (int) current year in simulation
    :param n: (int) number of cells to consider
    :param output_list: (boolean) Type of output by zone. If True output is list of cells. If false, output is dataframe with cells and annual flows
    :return: dictionary with a list/dataframe for each zone (depending on output_list)
    """
    gcm = curtailparam.listgcms[0]
    allCellFolders = get_all_cells_from_netcdf(curtailparam, gcm)

    df_final = None

    # read file with dictionary mapping cells to IPM zones (this file must be in the VICS/RBM folder)
    with open(os.path.join(curtailparam.rbmRootDir, 'cells2zones.pk'), 'rb') as f:
        cells2zones = pk.load(f)

    for gcm in curtailparam.listgcms:

        if genparam.referenceCase:
            # in reference case just use flows of an arbitrary gcm
            gcm = 'bcc-csm1-1.RCP45'
        else:
            # reading flow/streamT file. substitute '_' to '.'
            gcm = gcm.replace('_', '.')
            gcm = gcm.replace('rcp', 'RCP')

        fname = curtailparam.basenameflow
        fname = os.path.join(curtailparam.rbmRootDir, fname.format(gcm))
        streamflow, date, lons, lats = read_waterdata_netcdf(fname, 'streamflow', currYear)

        avg_stream_flow = np.mean(streamflow.data[:, :, :], axis=0)

        for z in genparam.ipmZones:
            if allCellFolders is not None:
                cellFoldersInZone = [c for c in allCellFolders if cells2zones[c] == z]
            else:
                cellFoldersInZone = [c for c in cells2zones if cells2zones[c] == z]

            avg_stream_flow_in_zone = []

            for c in cellFoldersInZone:
                cellLat, cellLon = tuple(map(float, c.split('_')))

                ix = np.argwhere(lats == cellLat).flatten()[0]
                iy = np.argwhere(lons == cellLon).flatten()[0]

                avg_stream_flow_in_zone = avg_stream_flow_in_zone + [avg_stream_flow[ix, iy]]

            # add to pandas data frame to sort
            if df_final is None:
                df_final = pd.DataFrame({'zone': z, 'gcm': gcm, 'cell': cellFoldersInZone, 'flow': avg_stream_flow_in_zone})
            else:
                df_final = pd.concat([df_final,
                                      pd.DataFrame({'zone': z, 'gcm': gcm, 'cell': cellFoldersInZone,
                                                    'flow': avg_stream_flow_in_zone})])

    # compute mean for each cell over all gcms
    a = df_final.groupby(['cell'], as_index=False).agg({'zone': 'first', 'flow': 'mean'})

    best_cells_by_zone = {}
    for z in genparam.ipmZones:
        df = a.groupby(['zone'], as_index=False).get_group(z).sort_values(by=['flow'], ascending=False)

        df = df.head(n)

        if output_list:
            best_cells = list(df['cell'])
        else:
            best_cells = df

        best_cells_by_zone[z] = best_cells

    return best_cells_by_zone


def convert_2dList_netcdf(listcurtail, curtailparam, fnameuw, fnameout='~/test.nc', ):
    """This function converts a 2d list with new generator curtailments to a netcdf file


    :param listcurtail: (2d list) with curtailment data
    :param curtailparam: object of type :mod:`Curtailmentparameters`
    :param fnameuw: (string) name of netcdf file with original meteo data from UW (no path)
    :param fnameout: (string) name of resulting netcdf file that will be created (full path)
    """

    prec = curtailparam.locPrecision

    # read netcdf file with meteo data form UW to get spatial limits
    dataset = nc.Dataset(os.path.join(curtailparam.rbmRootDir, fnameuw))
    time = dataset.variables['time'][:]
    lats = dataset.variables['lat'][:]
    lons = dataset.variables['lon'][:]
    temp = dataset.variables['temp'][:]

    dataset.close()

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


def expand_df_hourly(df):
    """utility function to expand a daily data frame to hourly

    :param df: (pandas data frame) data frame with daily data (check names of columns)
    :return: data frame with hourly data
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


def createBaseFilenameToReadOrWrite(locPrecision, inputLat, inputLong):
    """Creates string with name of folder with cell data

    This function creates a string with the formatted name of the folder that contains the data
    for the respective grid cell

    :param locPrecision: number of decimal digits in lat and long values
    :param inputLat: latitude of grid cell
    :param inputLong: longitude of grid cell
    :return: string with name of folder (e.g. 34.4375\_-86.4375)
    """
    return '%.*f_%.*f' % (locPrecision, inputLat, locPrecision, inputLong)


def getCellLatAndLongFromFolderName(dummyFolder):
    """Get cell lat and long values from formatted string

    This function splits a string with lat and long (format: '{lat}_{lon}') into numeric lat and long values

    See :func:`.createBaseFilenameToReadOrWrite`

    :param dummyFolder: (string) string with lat an long values
    :return: tuple with numeric values of lat and long
    """
    cellLat, cellLon = dummyFolder.split('_')
    return float(cellLat), float(cellLon)


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


def read_bulk_water_meteo(currYear, curtailparam):

    t0 = time.time()
    basenamemeteo = curtailparam.basenamemeteo

    dataset = nc.Dataset(os.path.join(curtailparam.rbmRootDir, basenamemeteo.format(currYear)))

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

    dataset.close()

    # convert to masked arrays
    waterT = np.ma.array(waterT, mask=temp.mask)
    flow = np.ma.array(flow, mask=temp.mask)

    print(str_elapsedtime(t0))