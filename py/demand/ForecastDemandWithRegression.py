# Michael Craig
# August 3, 2017
# Script forecasts demand using regressions devleoped by Francisco

import os, copy
import numpy as np
import pandas as pd
from collections import OrderedDict
from AuxFuncs import *
from demand.TempTransformation import *
from thermalderatings.AuxCurtailmentFuncs import read_netcdf
from thermalderatings.ModifyGeneratorCapacityWithWaterTData import find125GridMaurerLatLong


def forecastZonalDemandWithReg(yr, genparam, curtailparam):
    """ Forecast demand for given year using linear regression parameters

    This function uses the parameters from the previously fitted regressions to estimated future hourly load at each zone.

    :param yr: (integer) Current year of analysis
    :param genparam: object of class Generalparameters
    :param curtailparam: object of class Curtailmentparameters

    :return: - **zonalDemand:** nested dictionary with hourly load data for each zone in current year in each gcm. {gcm:zone:[hourly demand]}
             - **zonalTempDfs:** nested dictionary with data frames with meteo and load data for each zone in current year in each gcm. {gcm: zone: Df}
    """

    totalDemandDict, totalDemandDictDf = OrderedDict(), OrderedDict()

    if not genparam.incDemandCC:

        list_percent_local = copy.deepcopy(genparam.gcmranking)
        if np.min(list_percent_local) < 0 or np.max(list_percent_local) > 1:
            print('percentile values are not in the range [0, 1]!')
            print('Assuming they are ranks of the GCMs. Dividing by the number of GCMs')
            list_percent_local = [x / 20 for x in list_percent_local]

        # if running case without impacts in demand, gcmranking is list with PERCENTILES!
        df_demand_reference = getDemandReference(genparam, list_percent=list_percent_local, yearFixed=2015)

        for i, yr in enumerate(df_demand_reference['yearorig'].unique()):

            df_year = df_demand_reference[df_demand_reference['yearorig'] == yr]
            zonalDemand, zonalTempDfs = OrderedDict(), OrderedDict()

            for zone in genparam.ipmZones:
                zonalDemand[zone] = list(df_year.loc[df_year['zone'] == zone, 'load(MW)'].values)
                zonalTempDfs[zone] = df_year.loc[df_year['zone'] == zone]

            #key = 'ref{0:02d}'.format(i)
            key = curtailparam.listgcms[i]
            totalDemandDict[key], totalDemandDictDf[key] = zonalDemand, zonalTempDfs

    else:
        for (idx_gcm, gcm) in enumerate(curtailparam.listgcms):
            zonalDemand, zonalTempDfs = OrderedDict(), OrderedDict()

            for zone in genparam.ipmZones:

                (dataYr, tempCoefs, intCoefs, fixEffHr, fixEffYr, intercept,
                 holidays) = loadRegData(genparam, yr, zone, curtailparam, idx_gcm, netcdf=True)

                addTimeDummies(dataYr, str(yr), holidays)
                predictDemand(dataYr, tempCoefs, intCoefs, fixEffHr, fixEffYr, intercept, yr)
                #dataYr.to_csv(os.path.join(genparam.resultsDir, 'demandAndMetDf' + zone + str(yr) + '.csv'))

                zonalDemand[zone], zonalTempDfs[zone] = list(dataYr['load(MW)'].values), dataYr

            totalDemandDict[gcm], totalDemandDictDf[gcm] = zonalDemand, zonalTempDfs

    return totalDemandDict, totalDemandDictDf


def loadRegData(genparam, currYear, zone, curtailparam, idx_gcm, netcdf=False):
    """Load all necessary data into pandas DFs except intercept (just return value)

    :param genparam: object of class Generalparameters
    :param currYear: (int) current year of simulation
    :param zone: (string) ipm zone being simulated
    :param curtailparam: object of class Curtailmentparameters
    :param idx_gcm: (int) index of GCM being considered (within curtailparam.listgcms)
    :param netcdf: (boolean) true if format of meteo data is netcdf

    :return: - data: dataframe with weather data
             - tempCoefs: dataframe with temperature coefficients
             - intCoefs: dataframe with interactions coefficients (generally not used)
             - fixEffHr: dataframe with fixed effects for each hour of day
             - fixEffYr: data frame with fixed effects of each year
             - intercept: dataframe with intercept
             - holidays: dataframe with list of holidays
    """

    dataDir = os.path.join(genparam.dataRoot, 'DemandData')
    gcm = curtailparam.listgcms[idx_gcm]

    if netcdf:
        cellLat, cellLon = getCellForZone(dataDir, zone)
        cellLat, cellLon = find125GridMaurerLatLong(cellLat, cellLon)
        data = read_netcdf(cellLat, cellLon, currYear, curtailparam, gcm)

    else:
        data = pd.read_table(os.path.join(dataDir, 'meteo_memphis'), names=['precip', 'sh', 'rh', 'tC', 'p'])
        data = isolateYrData(data, currYear)

    # a few negative values of RH have appeared in the past. To avoid errors when computing demand set negative
    # values to 0.1%
    if np.sum((data['rh'] <= 0)) > 0:
        print('**WARNING:** There are negative values of RH in the data set. They were changed to 0.1%.')
        data['rh'] = data['rh'].apply(lambda x: 0.1 if x <= 0 else x)

    tempCoefs = pd.read_csv(os.path.join(dataDir, 'temp_coef_{}.csv'.format(zone)))
    fixEffHr = pd.read_csv(os.path.join(dataDir, 'fix_eff_hour_{}.csv'.format(zone)))
    fixEffYr = pd.read_csv(os.path.join(dataDir, 'fix_eff_year_{}.csv'.format(zone)))
    intercept = float(readCSVto2dList(os.path.join(dataDir, 'intercept_{}.csv'.format(zone)))[1][0])

    # check if interactions file exists
    if os.path.isfile(os.path.join(dataDir, 'interaction_coef_{}.csv'.format(zone))):
        intCoefs = pd.read_csv(os.path.join(dataDir, 'interaction_coef_{}.csv'.format(zone)))
    else:
        intCoefs = None

    holidays = pd.read_csv(os.path.join(dataDir, 'holiday_nerc.csv'))

    return data, tempCoefs, intCoefs, fixEffHr, fixEffYr, intercept, holidays


def getCellForZone(dataDir, zone):
    """This functions reads a file that maps each zone to a grid cell that has the METEO data that will be used

    :param dataDir: path to folder with demand data
    :param zone: string with symbol for subzone in SERC
    :return: tuple with cellLat, cellLon
    """

    map_zone_cell = pd.read_csv(os.path.join(dataDir, 'map_zone_cell.csv'))

    cellLat = map_zone_cell[map_zone_cell['zone'] == zone].iloc[0]['cellLat']
    cellLon = map_zone_cell[map_zone_cell['zone'] == zone].iloc[0]['cellLon']

    return cellLat, cellLon


def getStationForZone(dataDir, zone):
    """This functions reads a file that maps each zone to a representative city that has the METEO data

    :param dataDir: path to folder with demand data
    :param zone: string with symbol for subzone in SERC
    :return: name of city used for this zone
    """

    map_zone_cell = pd.read_csv(os.path.join(dataDir, 'map_zone_cell.csv'))

    city = map_zone_cell[map_zone_cell['zone'] == zone].iloc[0, 3]

    return city


def isolateYrData(data, yr):
    """Add datetime to each row of data, then isolate year and data of interest.

    :param data: dataframe with weather data
    :param yr: (int) current year
    :return: dataframe with data of current year
    """
    startDt, endDt = '1/1/1950 00:00:00', '12/31/2099 23:00:00'
    data['date'] = pd.date_range(start=startDt, end=endDt, freq='H')  # ,tz='US/Central')
    dataYr = data[data['date'].dt.year == yr][['date', 'rh', 'tC']]
    dataYr.reset_index(inplace=True)
    return dataYr


def addTimeDummies(dataYr, yr, holidays):
    """Add indicators for day of week, type of day (weekend versus weekday), and season to DF.

    :param dataYr: dataframe with data of current year
    :param yr: (int) current year
    :param holidays: dataframe with holidays
    """
    # Add whether weekday or weekend
    addWeekdayOrWeekend(dataYr, holidays)
    # Add which season
    seasons = {'Spring': ('3/15/' + yr, '6/15/' + yr), 'Fall': ('9/15/' + yr, '12/15/' + yr),
               'Summer': ('6/15/' + yr, '9/15/' + yr)}
    dataYr['season'] = 'Winter'  # use winter - don't need to worry about wrpa-around in yr

    for season in seasons:
        startDt, endDt = pd.to_datetime(seasons[season][0]), pd.to_datetime(seasons[season][1])
        dataYr.loc[(dataYr['date'] >= startDt) & (dataYr['date'] < endDt), 'season'] = season

    # Add hour of day
    dataYr['hour.of.day'] = dataYr['date'].dt.hour


def addWeekdayOrWeekend(dataYr, holidays):
    """Adds column that indicates whether weekday or weekend. Import list of holidays from Francisco that also are
    labeled as weekends.

    :param dataYr: dataframe with data of current year
    :param holidays: dataframe with holidays
    """
    dataYr['dayofweek'] = dataYr['date'].dt.dayofweek
    dataYr['type.day'] = 'weekday'
    dataYr.loc[dataYr['dayofweek'] >= 5, 'type.day'] = 'weekend'
    holidaysDt = pd.to_datetime(holidays['holiday']).dt.date
    dataYr.loc[dataYr['date'].dt.date.isin(holidaysDt), 'type.day'] = 'weekend'


def predictDemand(dataYr, tempCoefs, intCoefs, fixEffHr, fixEffYr, intercept, yr):
    """ Predict hourly demand using regression function

    Adds a column 'load(MW)' to data frame dataYr with simulated hourly load

    Predict value as:

    .. math::
       y=\\beta_0+\\beta_1*T_{bin}+\\beta_2*[T_{bin}*dp]+FE_{yr}+FE_{hr}

    All subfunctions return np array of 8760x1.

    :param dataYr: data frame with data for regression
    :param tempCoefs: temperature coefficients
    :param intCoefs: interaction coefficients
    :param fixEffHr: fixed effect for hour
    :param fixEffYr: fixed effect for year
    :param intercept: intercept of regression equation
    :param yr: current year
    """
    tempVals = getTempVals(dataYr, tempCoefs)

    if intCoefs is not None:
        interactionVals = getInteractionVals(dataYr, intCoefs)
    else:
        interactionVals = 0

    interceptVals = np.repeat(np.array(intercept), dataYr.shape[0])
    fixEffHrVals = getFixEffHrVals(dataYr, fixEffHr)
    fixEffYrVals = getFixEffYrVal(yr, fixEffYr, dataYr)

    dataYr['load(MW)'] = tempVals + interactionVals + interceptVals + fixEffHrVals + fixEffYrVals


def getTempVals(dataYr, tempCoefs):
    """Estimate temperature dependent load values by first putting T into bins (2d array), then multiplying by per-bin coefficients.

    :param dataYr: data frame with data for regression
    :param tempCoefs: dataframe with coefficients for temperature
    :return: numoy array of estimated load values
    """

    if 'tC' in dataYr.columns:
        temps = dataYr['tC']
    else:
        temps = dataYr['airT']

    # parse bins
    bins = [float(a.split(',')[1].replace(']', '').replace(')', '')) for a in np.array(tempCoefs['bin'])]
    # remove +infinity (last term)
    bins = bins[:-1]

    tempBins = createTempComponents(temps, bins)

    return np.dot(tempBins, np.array(tempCoefs['value']))


def getInteractionVals(dataYr, intCoefs):
    """Element-wise of binned T by dew pt, then dot w/ coefficients

    :param dataYr: data frame with data for regression
    :param intCoefs: dataframe with coefficients for interaction temp * humidity
    :return:
    """

    if 'tC' in dataYr.columns:
        temps, rhs = dataYr['tC'], np.array(dataYr['rh'])
    else:
        temps, rhs = dataYr['airT'], np.array(dataYr['rh'])

    tempBins = createTempComponents(temps)
    dewPt = convertRelHum2DewPoint(rhs, temps)
    return np.dot(np.multiply(dewPt, tempBins.T).T, np.array(intCoefs['value']))


def getFixEffHrVals(dataYr, fixEffHr):
    """Add FE hour vals (by season, time of day, and type of day) to DF, then return FEs

    :param dataYr: data frame with data for regression
    :param fixEffHr: dataframe with fixed effect for hour of day
    :return:
    """
    dataYrWithFE = dataYr.merge(fixEffHr, how='left', on=['season', 'type.day', 'hour.of.day'])
    return np.array(dataYrWithFE['value'])


def getFixEffYrVal(yr, fixEffYr, dataYr):
    """Return np array w/ FE for yr.

    :param yr: current year
    :param fixEffYr: estimated fixed effects
    :param dataYr: pandas data frame
    :return: np array with fixed effect for year yr.
    """
    if yr in fixEffYr['year'].values:
        val = fixEffYr.loc[fixEffYr['year'] == yr, 'value']
    else:
        # if year not included in training data set, use fixed effect for last year of training data
        val = fixEffYr.tail(1).iloc[0]['value']

    return np.repeat(np.array(val), dataYr.shape[0])


def getDemandReference(genparam, list_percent=(0.2, 0.5, 0.8), yearFixed=2015):

    df_dem_total = None

    for zone in genparam.ipmZones:

        start_time = time.time()
        print('Computing reference hourly demand for zone {} ...'.format(zone), flush=True)

        dataDir = os.path.join(genparam.dataRoot, 'DemandData')

        city = getStationForZone(dataDir, zone)
        c = pd.read_csv(os.path.join(dataDir, 'training', 'RH_airT_19790101_20161231.{}.csv'.format(city)))
        c = c.rename(columns={c.columns[0]: 'datetime'})
        c = c.rename(index=str, columns={'datetime': 'date', 'air_T': 'airT', 'rel_humid': 'rh'})
        c['date'] = pd.to_datetime(c['date'])
        c['yearorig'] = c['date'].dt.year
        c['year'] = yearFixed

        tempCoefs = pd.read_csv(os.path.join(dataDir, 'temp_coef_{}.csv'.format(zone)))
        fixEffHr = pd.read_csv(os.path.join(dataDir, 'fix_eff_hour_{}.csv'.format(zone)))
        fixEffYr = pd.read_csv(os.path.join(dataDir, 'fix_eff_year_{}.csv'.format(zone)))
        intercept = float(readCSVto2dList(os.path.join(dataDir, 'intercept_{}.csv'.format(zone)))[1][0])

        intCoefs = None

        holidays = pd.read_csv(os.path.join(dataDir, 'holiday_nerc.csv'))

        addTimeDummies(c, str(yearFixed), holidays)

        tempVals = getTempVals(c, tempCoefs)

        interactionVals = 0

        interceptVals = np.repeat(np.array(intercept), c.shape[0])
        fixEffHrVals = getFixEffHrVals(c, fixEffHr)
        fixEffYrVals = getFixEffYrVal(yearFixed, fixEffYr, c)

        c['load(MW)'] = tempVals + interactionVals + interceptVals + fixEffHrVals + fixEffYrVals
        c['zone'] = zone

        if df_dem_total is None:
            df_dem_total = c
        else:
            df_dem_total = pd.concat([df_dem_total, c])

        print('Done! Elapsed time: {}'.format(str_elapsedtime(start_time)))

    # downcast strings to categorical
    df_dem_total['zone'].astype('category')

    df_dem_total = df_dem_total[['date', 'zone', 'rh', 'airT', 'yearorig', 'year', 'load(MW)']]

    df_dem_total_2 = df_dem_total.groupby(['date']).agg({'load(MW)': np.sum}).reset_index()
    df_dem_total_2['year'] = df_dem_total_2['date'].dt.year

    df_dem_total_3 = df_dem_total_2.groupby(['year']).apply(lambda x: x[x['load(MW)'] == np.max(x['load(MW)'])]).reset_index(drop=True)

    # get quantiles of demand
    load_quantiles = df_dem_total_3['load(MW)'].quantile(list_percent, interpolation='nearest').values

    df_dem_total_4 = df_dem_total_3[df_dem_total_3['load(MW)'].isin(load_quantiles)]

    df_dem_total_5 = df_dem_total[df_dem_total['yearorig'].isin(df_dem_total_4['year'].values)]

    return df_dem_total_5
