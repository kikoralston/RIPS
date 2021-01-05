# Michael Craig, 6 July 2016
# Get wind and solar capacity factors
#
# This script gets wind and solar capacity factors for existing and new renewables
# by mapping capacity to the best available units in each zone.

import os
from AuxFuncs import *
from ipmzones.AssignCellsToIPMZones import (locInZone)
import pandas as pd
import datetime as dt
import multiprocessing as mp
import gc


def wrapperwind(list_args):
    """ wrapper function to read and process wind data in parallel

    :param list_args: a list with arguments for function getWindSiteCfs_2
    :return: list with:
                    - 1d list with hourly CFs for wind site
                    - 1d list with subhourly CFs for wind site
                    - string with siteId
                    - segment/block that this site is allocated
                    - original capacity of site in Wind data base (MW)
                    - capacity of site used for fleet computation (MW)
    """
    (windDir, siteId, datasetCapac, fleetCapac, segment, desiredTz, windGenDataYr, subHour) = list_args

    siteCfsHourly, siteCfsSubhourly = getWindSiteCfs_2(windDir, siteId, datasetCapac, desiredTz, windGenDataYr,
                                                       subHour=subHour)

    out = [siteCfsHourly, siteCfsSubhourly, siteId, segment, datasetCapac, fleetCapac]

    return out


def wrappersolar(site):
    """ wrapper function to read and process solar data in parallel

    :param site: a list with arguments for function getSolarSiteCfs_2
    :return: list with:
                    - 1d list with hourly CFs for wind site
                    - 1d list with subhourly CFs for wind site
                    - string with site file name
                    - segment/block that this site is allocated
                    - original capacity of site in solar data base (MW)
                    - capacity of site used for fleet computation (MW)
    """
    (solarDir, siteFilename, datasetSiteCapac, fleetCapac, segment, siteTz, desiredTz, subHour) = site

    siteCfsHourly, siteCfsSubhourly = getSolarSiteCfs_2(solarDir=solarDir, siteFilename=siteFilename,
                                                        datasetSiteCapac=datasetSiteCapac, siteTz=siteTz,
                                                        desiredTz=desiredTz, subHour=subHour, listout=True)

    out = [siteCfsHourly, siteCfsSubhourly, siteFilename, segment, datasetSiteCapac, fleetCapac]

    return out


def getRenewableCFData(currZone, genparam, sizeSegment=1000, fleetCap=70, capacInCurrFleet=0, type='wind',
                       existing=False, subHour = False):
    """Compile wind/solar potential data

    This is the main function to compile wind/solar potential generation data for both existing and new solar/wind
    generators. It reads data from the databases and prepares a dataframe with hourly generation potential

    :param currZone: (string) name of ipm zone being simulated
    :param genparam: object of type :mod:`Generalparameters`
    :param sizeSegment: (int) size of segment (in MW)
    :param fleetCap: (float) capacity in MW of each site
    :param capacInCurrFleet: (float) existing total installed capacity in MW of wind or solar
    :param type: (string) 'wind' or 'solar'
    :param existing: (boolean) If True process data for existing power plants. If False process data for new
                     candidate plants
    :return: two data frame with hourly data and subhourly data

    """
    t_start = time.time()
    t_step = time.time()

    windGenDataYr = genparam.windGenDataYr
    desiredTz = genparam.tzAnalysis

    if type.lower() == 'wind':
        renewableDir = os.path.join(genparam.dataRoot, 'WINDSERCData')
        metadata = pd.read_csv(os.path.join(renewableDir, 'toolkit_sites_v7_SERC_zones.csv'),
                               dtype={'sitenumbers': np.int})
        metadata.rename(columns={'capacity_factor': 'cfs', 'capacity': 'capacs', 'sitenumbers': 'site_id'},
                        inplace=True)
    else:
        renewableDir = os.path.join(genparam.dataRoot, 'NRELSolarPVData')
        metadata = pd.read_csv(os.path.join(renewableDir, 'SolarCapacityFactorsNRELSERC_zones.csv'))
        metadata.rename(columns={'CF': 'cfs', 'PlantSize': 'capacs', 'File': 'site_id'},
                        inplace=True)

    # read metada data of all sites in this zone and keep only relevant columns
    df_cfs = metadata[metadata['zone'] == currZone]
    df_cfs = df_cfs[['site_id', 'capacs', 'cfs']]

    # add column with size of site to be considered for the fleet computation
    df_cfs.loc[:, 'fleetCapac'] = fleetCap

    # sort sites according to CF and compute cumulative sum of capacity
    df_cfs = df_cfs.sort_values(by='cfs', ascending=False).reset_index(drop=True)
    df_cfs.loc[:, 'totalCapac'] = df_cfs['fleetCapac'].cumsum()

    # first define last block
    df_cfs_inv = df_cfs.sort_values(by=['cfs'], ascending=True)
    df_cfs_inv.loc[:, 'totalCapac'] = df_cfs_inv['fleetCapac'].cumsum()
    df_cfs_inv = df_cfs_inv[df_cfs_inv['totalCapac'] < sizeSegment]
    df_cfs_inv = df_cfs_inv.sort_values(by=['cfs'], ascending=False)

    #
    # If existing gen and current capacity == 0, just return data frame with zeros
    if existing and capacInCurrFleet == 0:

        # get array with times and dates from first file (and cap it at 8760 hours)
        if type.lower() == 'wind':
            dfauxhour, dfauxsub = getWindSiteCfs_2(windDir=renewableDir, siteId=str(df_cfs_inv['site_id'].iloc[0]),
                                                   siteCapac=df_cfs_inv['capacs'].iloc[0], desiredTz=desiredTz,
                                                   windGenDataYr=windGenDataYr, subHour=subHour, listout=False)
        else:
            sitetz = timezoneOfSolarSite(df_cfs_inv['site_id'].iloc[0], currZone)
            dfauxhour, dfauxsub = getSolarSiteCfs_2(solarDir=renewableDir,
                                                    siteFilename=str(df_cfs_inv['site_id'].iloc[0]),
                                                    datasetSiteCapac=df_cfs_inv['capacs'].iloc[0], siteTz=sitetz,
                                                    desiredTz=desiredTz, subHour=subHour, listout=False)

        dtvalues_hour = dfauxhour.iloc[:8760, 0].tolist()
        df_hour = pd.DataFrame({'datetime': dtvalues_hour, 'segment': 0, 'cfs': 0, 'gen': 0})
        df_hour = df_hour[['datetime', 'segment', 'cfs', 'gen']]

        if subHour:
            dtvalues_sub = dfauxsub.iloc[:, 0].tolist()
            df_subhour = pd.DataFrame({'datetime': dtvalues_sub, 'segment': 0, 'cfs': 0, 'gen': 0})
            df_subhour = df_subhour[['datetime', 'segment', 'cfs', 'gen']]
        else:
            df_subhour = None

        print('    FINISHED: ' + str_elapsedtime(t_start))

        # ****** EXIT FUNCTION ********
        return df_hour, df_subhour

    if existing:

        # Subtract existing capacity from cumulative capacity in data set
        df_cfs.loc[:, 'totalCapac'] = df_cfs['totalCapac'] - capacInCurrFleet

        if df_cfs.iloc[-1, df_cfs.columns.get_loc('totalCapac')] < 0:
            # if total capacity of existing sites is greater than capacity in dataset
            # fill up this delta value with CFs from last block of sites
            # (the idea is that additional sites will be built using CFs form this last block)

            deltaCapac = -1*df_cfs.iloc[-1, df_cfs.columns.get_loc('totalCapac')]
            df_cfs_inv['fleetCapac'] = deltaCapac/df_cfs_inv.shape[0]

            df_cfs = pd.concat([df_cfs, df_cfs_inv])
            df_cfs.loc[:, 'totalCapac'] = df_cfs['fleetCapac'].cumsum()

        else:
            # only keep sites that belong to existing plants
            df_cfs = df_cfs[df_cfs['totalCapac'] < fleetCap]

            # adjust capacity of last site
            df_cfs.iloc[-1, df_cfs.columns.get_loc('fleetCapac')] = fleetCap - df_cfs['totalCapac'].iloc[-1]

        # assign all sites to single block/segment
        df_cfs['segment'] = 0

    else:
        # remove sites that belong to existing plants
        df_cfs.loc[:, 'totalCapac'] = df_cfs['totalCapac'] - capacInCurrFleet
        df_cfs = df_cfs[df_cfs['totalCapac'] > 0]

        if df_cfs.shape[0] > 0:
            # adjust capacity of first site
            df_cfs.iloc[0, df_cfs.columns.get_loc('fleetCapac')] = df_cfs['totalCapac'].iloc[0]

            # assign sites to blocks/segments
            df_cfs['segment'] = df_cfs['totalCapac'].apply(lambda x: x // sizeSegment)
        else:
            # no more new sites, just use capacity factors of last block
            df_cfs = df_cfs_inv
            df_cfs['segment'] = 0

    print('    Finished reading site CF data: ' + str_elapsedtime(t_step))

    t_step = time.time()

    ncores = genparam.ncores_py
    gc.collect()

    if type.lower() == 'wind':
        list_args = [[renewableDir, str(int(site[1]['site_id'])), site[1]['capacs'], site[1]['fleetCapac'], site[1]['segment'],
                      desiredTz, windGenDataYr, subHour] for site in df_cfs.iterrows()]
        if ncores == 1:
            list_cfs = list(map(wrapperwind, list_args))
        else:
            with mp.Pool(processes=ncores) as pool:
                list_cfs = pool.map(wrapperwind, list_args)

    else:
        list_args = [[renewableDir, site[1]['site_id'], site[1]['capacs'], site[1]['fleetCapac'], site[1]['segment'],
                      timezoneOfSolarSite(site[1]['site_id'], currZone), desiredTz, subHour]
                     for site in df_cfs.iterrows()]
        if ncores == 1:
            list_cfs = list(map(wrappersolar, list_args))
        else:
            with mp.Pool(processes=ncores) as pool:
                list_cfs = pool.map(wrappersolar, list_args)

    gc.collect()

    print('    Finished reading all sites: ' + str_elapsedtime(t_step))

    t_step = time.time()

    # first compute average CFs for hourly gen
    df_hour = None
    for s in df_cfs['segment'].unique():

        list_cfs_segment = [siteCfs for siteCfs in list_cfs if siteCfs[3] == s]

        # convert to data frame
        for i, x in enumerate(list_cfs_segment):

            if df_hour is None:
                df_hour = pd.DataFrame({'datetime': x[0][0][1:], 'Id': x[2], 'capac': x[5],'cfs': x[0][1][1:], 'segment': s})
            else:
                df_hour = df_hour.append(pd.DataFrame({'datetime': x[0][0][1:], 'Id': x[2], 'capac': x[5], 'cfs': x[0][1][1:],
                                                       'segment': s}), ignore_index=True)

    # compute generation in MWh
    df_hour['gen'] = df_hour['capac'] * df_hour['cfs']

    # compute weighted average CF using vectorization and pandas
    x = df_hour.groupby(['segment', 'datetime'])['capac'].agg(np.sum).reset_index()
    df_aux = pd.merge(df_hour, x, how='left', on=['segment', 'datetime'])
    df_aux['cfsw'] = df_aux['cfs']*df_aux['capac_x']/df_aux['capac_y']

    df_hour = df_aux.groupby(['segment', 'datetime'])['cfsw', 'gen'].agg(np.sum).reset_index().rename(columns={'cfsw': 'cf'})
    # just to make sure order of hours is correct is correct
    df_hour = df_hour.sort_values(by=['segment', 'datetime'])

    del x
    del df_aux
    gc.collect()

    # then compute average CFs for sub-hourly gen
    df_subhour = None
    if subHour:
        for s in df_cfs['segment'].unique():

            list_cfs_segment = [siteCfs for siteCfs in list_cfs if siteCfs[3] == s]

            # convert to data frame
            for i, x in enumerate(list_cfs_segment):

                if df_subhour is None:
                    df_subhour = pd.DataFrame({'datetime': x[1][0][1:], 'Id': x[2], 'capac': x[5],'cfs': x[1][1][1:], 'segment': s})
                else:
                    df_subhour = df_subhour.append(pd.DataFrame({'datetime': x[1][0][1:], 'Id': x[2], 'capac': x[5],
                                                                 'cfs': x[1][1][1:], 'segment': s}), ignore_index=True)

        # compute generation in MWh
        df_subhour['gen'] = df_subhour['capac'] * df_subhour['cfs']

        # compute weighted average CF using vectorization and pandas
        x = df_subhour.groupby(['segment', 'datetime'])['capac'].agg(np.sum).reset_index()
        df_aux = pd.merge(df_subhour, x, how='left', on=['segment', 'datetime'])
        df_aux['cfsw'] = df_aux['cfs']*df_aux['capac_x']/df_aux['capac_y']

        df_subhour = df_aux.groupby(['segment', 'datetime'])['cfsw', 'gen'].agg(np.sum).reset_index().rename(columns={'cfsw': 'cf'})
        # just to make sure order of hours is correct is correct
        df_subhour = df_subhour.sort_values(by=['segment', 'datetime'])

    print('    Finished computing average CF using vectorization: ' + str_elapsedtime(t_step))

    print('    FINISHED: ' + str_elapsedtime(t_start))

    # return both data frames (let the routine that is calling it parse it to the desired format)
    return df_hour, df_subhour


def getPlantInfoInZone(metadata, cfCol, capacCol, siteNumberOrFileCol, fipsToZones, fipsToPolys, currZone,
                       return_df=False):
    """Get all meta data of renewable sites in current zone

    Match by zone

    :param metadata: 2-d list with metadata of solar or wind sites
    :param cfCol: (int) index of column with CF data
    :param capacCol: (int) index of column with capacity data (in MW)
    :param siteNumberOrFileCol: (int) index of column with site ID
    :param fipsToZones: (dict) dictionary mapping FIPS to zones
    :param fipsToPolys: (dict) dictionary mapping FIPS to poly shapes
    :param currZone: (string) ipm zone
    :param return_df: if True return data frame as output with columns cfs, capacs, sitenumbers.
                      If False returns tuple with three 1-d lists cfs, capacs, sitenumbers.
    :return: a data frame with three columns or a tuple with 3 lists:
             - CFs
             - capacs
             - sites IDs
    """

    if 'longitude' in metadata[0]:
        # wind has lat/long already given

        latCol, lonCol = metadata[0].index('latitude'), metadata[0].index('longitude')
        plantsInRegionOrZone = [row for row in metadata[1:] if locInZone(float(row[latCol]), float(row[lonCol]),
                                                                         currZone, fipsToZones, fipsToPolys)]
    else:
        # solar lat/lon is in filename

        nameCol = metadata[0].index('File')
        plantsInRegionOrZone = [row for row in metadata[1:] if locInZone(float(getCoordsFromFilename(row[nameCol])[0]),
                                                                         float(getCoordsFromFilename(row[nameCol])[1]),
                                                                         currZone, fipsToZones, fipsToPolys)]

    cfs = [float(row[cfCol]) for row in plantsInRegionOrZone]
    capacs = [float(row[capacCol]) for row in plantsInRegionOrZone]
    siteNumbers = [row[siteNumberOrFileCol] for row in plantsInRegionOrZone]

    if return_df:
        out = pd.DataFrame({'cfs': cfs, 'capacs': capacs, 'sitenumbers': siteNumbers})
    else:
        out = (cfs, capacs, siteNumbers)

    return out


def getWindSiteCfs_2(windDir, siteId, siteCapac, desiredTz, windGenDataYr, subHour=False, listout=True):
    """Read wind data for a single site

    optimized version of function getWindSiteCfs using pandas (75% faster)

    :param windDir: dir w/ wind data
    :param siteId: site ID to get gen data for
    :param siteCapac: wind site capac
    :param desiredTz: desired timezone
    :param windGenDataYr: year for wind gen data
    :param subHour: if True reads sub hourly data
    :return: 2 2d lists, both have first row = datetime. 1 2d list = hourly CFs, 2nd 2d list = subhourly CFs.
             Also row labels
    """
    column_types = {'siteId': 'float'}
    tzOffsetDict = {'UTCtoCST': -6, 'CSTtoCST': 0, 'ESTtoCST': -1, 'CSTtoEST': 1, 'UTCtoEST': -5,
                    'ESTtoEST': 0}
    timezoneOffset = tzOffsetDict['UTC' + 'to' + desiredTz]

    hourlyFile = 'powerhourly_' + siteId + '.csv'

    hourlyGen = pd.read_csv(os.path.join(windDir, 'hourlyPowerSERC', hourlyFile), dtype=column_types,
                            parse_dates=['Datetime'], infer_datetime_format=True)

    # correct time zone
    hourlyGen['Datetime'] = hourlyGen['Datetime'] + dt.timedelta(hours=timezoneOffset)

    # filter to current year
    hourlyGen = pd.DataFrame(hourlyGen[hourlyGen['Datetime'].dt.year == windGenDataYr])
    hourlyGen = hourlyGen.reset_index(drop=True)

    # convert to cfs
    hourlyGen[siteId] = hourlyGen[siteId]/siteCapac

    # rename columns
    hourlyGen = hourlyGen.rename(columns={'Datetime': 'datetime{}'.format(desiredTz), siteId: 'power(MW){}'.format(siteId)})

    # convert to lists
    if listout:
        hourlyGen = convert_pandas_to_list(hourlyGen, header=list(hourlyGen.columns))

    if subHour:
        subhourlyFile = 'powersubhourly_' + siteId + '.csv'
        subhourlyGen = pd.read_csv(os.path.join(windDir, 'subhourlyPowerSERC', subhourlyFile), dtype=column_types,
                                   parse_dates=['Datetime'], infer_datetime_format=True)

        # correct time zone
        subhourlyGen['Datetime'] = subhourlyGen['Datetime'] + dt.timedelta(hours=timezoneOffset)

        # filter to current year
        subhourlyGen = pd.DataFrame(subhourlyGen[subhourlyGen['Datetime'].dt.year == windGenDataYr])
        subhourlyGen = subhourlyGen.reset_index(drop=True)

        # convert to cfs
        subhourlyGen[siteId] = subhourlyGen[siteId] / siteCapac

        # rename columns
        subhourlyGen = subhourlyGen.rename(columns={'Datetime': 'datetime{}'.format(desiredTz),
                                                    siteId: 'power(MW){}'.format(siteId)})

        # convert to list
        if listout:
            subhourlyGen = convert_pandas_to_list(subhourlyGen, header=list(subhourlyGen.columns))
    else:
        subhourlyGen = None

    return hourlyGen, subhourlyGen


def getSolarSiteCfs_2(solarDir, siteFilename, datasetSiteCapac, siteTz, desiredTz, subHour=False, listout=True):
    """Read solar data for a single site

    Optimized version of getSolarSiteCfs. This function uses Pandas in order to get more efficient processing

    :param solarDir:
    :param siteFilename:
    :param datasetSiteCapac:
    :param siteTz:
    :param desiredTz:
    :param listout: True if final result should be converted to list (same format as original function)
    :param subHour: True if subhourly values should be returned
    :return:
    """

    column_types = {'Power(MW)': 'float'}

    tzOffsetDict = {'UTCtoCST': -6, 'CSTtoCST': 0, 'ESTtoCST': -1, 'CSTtoEST': 1, 'UTCtoEST': -5,
                    'ESTtoEST': 0}
    timezoneOffset = tzOffsetDict[siteTz + 'to' + desiredTz]

    if os.path.exists(os.path.join(solarDir, 'SERC', siteFilename)):
        subhourly = pd.read_csv(os.path.join(solarDir, 'SERC', siteFilename), dtype=column_types)

        # convert from string to datetime
        subhourly['LocalTime'] = pd.to_datetime(subhourly['LocalTime'], format='%m/%d/%y %H:%M')

        # correct time zone
        subhourly['LocalTime'] = subhourly['LocalTime'] + dt.timedelta(hours=timezoneOffset)

        # convert to hourly
        subhourly['Date_hour'] = subhourly['LocalTime'].dt.floor(freq='H')
        hourly = subhourly.groupby(['Date_hour'], as_index=False).aggregate(np.mean)
        del subhourly['Date_hour']

        # change column name
        hourly = hourly.rename(columns={'Date_hour': 'datetimeEST'})

        #convert to cfs
        hourly_cf = pd.DataFrame(hourly)
        hourly_cf['cf'] = hourly_cf['Power(MW)']/datasetSiteCapac
        del hourly_cf['Power(MW)']

        # convert timestamps to datetime
        hourly_cf['datetimeEST'] = hourly_cf['datetimeEST'].dt.to_pydatetime()

        # convert to lists
        if listout:
            hourly_cf = convert_pandas_to_list(hourly_cf, header=['datetimeEST', 'power(MW){}'.format(siteFilename)])

        if subHour:
            # change column name
            subhourly = subhourly.rename(columns={'LocalTime': 'datetimeEST'})

            #convert to cfs
            subhourly_cf = pd.DataFrame(subhourly)
            subhourly_cf['cf'] = subhourly_cf['Power(MW)']/datasetSiteCapac
            del subhourly_cf['Power(MW)']

            # convert timestamps to datetime
            subhourly_cf['datetimeEST'] = subhourly_cf['datetimeEST'].dt.to_pydatetime()

            # convert to lists
            if listout:
                subhourly_cf = convert_pandas_to_list(subhourly_cf, header=['datetimeEST',
                                                                            'power(MW){}'.format(siteFilename)])
        else:
            subhourly_cf = None

    else:
        hourly_cf, subhourly_cf = [], []

    return hourly_cf, subhourly_cf


def convert_pandas_to_list(df, header=None):
    """ converts a pandas data frame to a list

    To maintain compatibility with older functions, uses same dimension convention as the previous version. It does
    not follow the ordinary dimension convention of pd.tolist(). Each column in the data frame is a list (pd.tolist()
    transforms each row to a list).

    :param df: data frame
    :param header: list of strings with header. If None, no header is added
    """

    df_list = df.values.T.tolist()

    if header is not None:
        df_new = []
        for i in range(len(df_list)):
            df_new = df_new + [[header[i]] + df_list[i]]

    return df_new


def getFleetToCapacDict(idAndCapac):
    idCol = idAndCapac[0].index('Id')
    fleetCapacCol = idAndCapac[0].index('FleetCapacity')
    idToFleetCapac = dict()
    for row in idAndCapac[1:]:
        idToFleetCapac[row[idCol]] = row[fleetCapacCol]
    return idToFleetCapac


def timezoneOfSolarSite(solarFilename, currZone):
    # in KY or TN, which is half & half CST & EST
    kyLine = [(36.601261, -84.861318), (38.048865, -86.251172)]
    tnLine = [(36.601261, -85.076648), (34.997791, -85.605184)]
    gaAlBorder = [(34.9894, -85.644), (29.893, -84.854)]
    if currZone == 'S_C_TVA' or currZone == 'S_C_KY' or currZone == 'S_SOU':
        if currZone == 'S_C_TVA':
            line = tnLine
        elif currZone == 'S_C_KY':
            line = kyLine
        elif currZone == 'S_SOU':
            line = gaAlBorder

        (siteLat, siteLong) = getCoordsFromFilename(solarFilename)

        if siteEastOfLine(line, float(siteLat), float(siteLong)):
            tz = 'EST'
        else:
            tz = 'CST'

    elif currZone == 'S_VACA':
        tz = 'EST'

    return tz


def getCoordsFromFilename(solarFilename):
    latStart = solarFilename.index('_') + 1
    latEnd = solarFilename[latStart:].index('_')
    lat = solarFilename[latStart:(latStart + latEnd)]
    longStart = solarFilename.index('-')
    longEnd = solarFilename[longStart:].index('_')
    longitude = solarFilename[longStart:(longStart + longEnd)]
    return (lat, longitude)


# Long = x coord, lat = y coord
def siteEastOfLine(line, siteLat, siteLong):
    (deltaLat, deltaLong) = (line[0][0] - line[1][0], line[0][1] - line[1][1])
    lineSlope = deltaLat / deltaLong
    intercept = line[1][0] - lineSlope * line[1][1]  # b = y - mx
    longOnLineForSiteLat = (siteLat - intercept) / lineSlope  # x = (y-b)/m
    return siteLong > longOnLineForSiteLat  # long decreases (more negative) west across US


def trimNewRECFsToCEHours(zonalNewWindCFs, zonalNewSolarCFs, hoursForCE):
    """Trim hours of new renewable CFs to hours of CE simulation

    :param zonalNewWindCFs: hourly CFs for new wind & solar builds (1d lists)
    :param zonalNewSolarCFs: hourly CFs for new wind & solar builds (1d lists)
    :param hoursForCE: hours included in CE (1d list, 1-8760 basis)
    :return: zonal hourly CFs for new wind and solar builds only for CE hours as dict of zone:
            [CFs (1-8760 basis, 1d list)]
    """
    newWindCFsCEZonal, newSolarCFsCEZonal = dict(), dict()

    for gcm in hoursForCE:
        # (hr - 1) b/c hours in year start @ 1, not 0 like Python idx
        newWindCFsCEZonal[gcm] = {zone: {'Wind+{0:02d}'.format(int(i)): [zonalNewWindCFs[zone][i][hr - 1]
                                                                         for hr in hoursForCE[gcm]]
                                         for i in zonalNewWindCFs[zone]} for zone in zonalNewWindCFs}

        newSolarCFsCEZonal[gcm] = {zone: {'Solar PV+{0:02d}'.format(int(i)): [zonalNewSolarCFs[zone][i][hr - 1]
                                                                              for hr in hoursForCE[gcm]]
                                          for i in zonalNewSolarCFs[zone]} for zone in zonalNewSolarCFs}

    #        newWindCFsCEZonal[gcm] = {zone: [zonalNewWindCFs[zone][i][hr - 1] for hr in hoursForCE[gcm]]
    #                                  for zone in zonalNewWindCFs}
    #        newSolarCFsCEZonal[gcm] = {zone: [zonalNewSolarCFs[zone][i][hr - 1] for hr in hoursForCE[gcm]]
    #                                   for zone in zonalNewSolarCFs}

    return newWindCFsCEZonal, newSolarCFsCEZonal
