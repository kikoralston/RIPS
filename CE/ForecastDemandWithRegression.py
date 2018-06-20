# Michael Craig
# August 3, 2017
# Script forecasts demand using regressions devleoped by Francisco

import os
import sys
import pandas as pd
from AuxFuncs import *
from TempTransformation import *


def forecastZonalDemandWithReg(yr, genparam):
    """ Forecast demand for given year using regression put together by Francisco

    This function runs the previously fitted regressions to estimated future hourly load at each zone. If the
    analysis area is 'test', then it uses user defined values of demand given in a file 'demand.csv' at the
    root of the data folder.

    :param yr: (integer) Current year of analysis
    :param genparam: object of class Generalparameters
    :return: zonalDemand: dictionary with hourly load data for each zone in current year
             zonalTempDfs: dictionary with data frames with meteo and load data for each zone in current year
    """
    zonalDemand, zonalTempDfs = dict(), dict()

    if genparam.analysisArea == 'test':
        # if analysisArea == 'test', assumes that there is a file 'demand.csv' at dataRoot, with hourly demand
        # for current year for each zone. the csv file must have one column for each zone
        # (first row has headers with zone names)

        fileNameWithDir = os.path.join(genparam.dataRoot, 'demand.csv')

        if os.path.isfile(fileNameWithDir):
            a = pd.read_csv(fileNameWithDir)
            # filter current year
            y = [int(date) == int(yr) for date in a['date']]
            a = a[y]
            a = a.reset_index(drop=True)
            del a['date']

            zonalDemand = a.to_dict(orient='list')
        else:
            print('Error! analysisArea is set to \'test\'. There must be a file \'demand.csv\' at the ' +
                  'root of the data folder')
            sys.exit()

    else:
        for zone in genparam.ipmZones:
            data, tempCoefs, intCoefs, fixEffHr, fixEffYr, intercept, holidays = loadRegData(genparam.dataRoot)
            dataYr = isolateYrData(data, yr)
            addTimeDummies(dataYr, str(yr), holidays)
            predictDemand(dataYr, tempCoefs, intCoefs, fixEffHr, fixEffYr, intercept, yr)
            dataYr.to_csv(os.path.join(genparam.resultsDir, 'demandAndMetDf' + zone + str(yr) + '.csv'))
            zonalDemand[zone], zonalTempDfs[zone] = list(dataYr['load(MW)'].values), dataYr

    return zonalDemand, zonalTempDfs


# Load all necessary data into pandas DFs exept inercept (just return value)
def loadRegData(dataRoot):

    dataDir = os.path.join(dataRoot, 'DemandData')

    data = pd.read_table(os.path.join(dataDir, 'meteo_memphis'), names=['precip', 'sh', 'rh', 'tC', 'p'])
    tempCoefs = pd.read_csv(os.path.join(dataDir, 'temp_coef.csv'))
    intCoefs = pd.read_csv(os.path.join(dataDir, 'interaction_coef.csv'))
    fixEffHr = pd.read_csv(os.path.join(dataDir, 'fix_eff_hour.csv'))
    fixEffYr = pd.read_csv(os.path.join(dataDir, 'fix_eff_year.csv'))
    intercept = float(readCSVto2dList(os.path.join(dataDir, 'intercept.csv'))[1][0])
    holidays = pd.read_csv(os.path.join(dataDir, 'holiday_nerc.csv'))
    return data, tempCoefs, intCoefs, fixEffHr, fixEffYr, intercept, holidays


# Add datetime to each row of data, then isolate year and data of interest.
def isolateYrData(data, yr):
    startDt, endDt = '1/1/1950 00:00:00', '12/31/2099 23:00:00'
    data['dt'] = pd.date_range(start=startDt, end=endDt, freq='H')  # ,tz='US/Central')
    dataYr = data[data['dt'].dt.year == yr][['dt', 'rh', 'tC']]
    dataYr.reset_index(inplace=True)
    return dataYr


# Add indicators for day of week, type of day (weekend versus weekday), and season to DF.
def addTimeDummies(dataYr, yr, holidays):
    # Add whether weekday or weekend
    addWeekdayOrWeekend(dataYr, holidays)
    # Add which season
    seasons = {'Spring': ('3/15/' + yr, '6/15/' + yr),
               'Fall': ('9/15/' + yr, '12/15/' + yr),
               'Summer': ('6/15/' + yr, '9/15/' + yr)}
    dataYr['season'] = 'Winter'  # use winter - don't need to worry about wrpa-around in yr
    for season in seasons:
        startDt, endDt = pd.to_datetime(seasons[season][0]), pd.to_datetime(seasons[season][1])
        dataYr.loc[(dataYr['dt'] >= startDt) & (dataYr['dt'] < endDt), 'season'] = season
    # Add hour of day
    dataYr['hour.of.day'] = dataYr['dt'].dt.hour


# Adds column that indicates whether weekday or weekend. Import list of holidays
# from Francisco that also are labeled as weekends.
def addWeekdayOrWeekend(dataYr, holidays):
    dataYr['dayofweek'] = dataYr['dt'].dt.dayofweek
    dataYr['type.day'] = 'weekday'
    dataYr.loc[dataYr['dayofweek'] >= 5, 'type.day'] = 'weekend'
    holidaysDt = pd.to_datetime(holidays['holiday']).dt.date
    dataYr.loc[dataYr['dt'].dt.date.isin(holidaysDt), 'type.day'] = 'weekend'


# Predict value as:
# y = beta*Tbin + alpha*Tbin*dewPt + FEyr + FEhr + intercept
# Note that Tbin*dewPt is element-wise, then multiplied via dot product into alpha.
# All subfucntions return np array of 8760x1.
def predictDemand(dataYr, tempCoefs, intCoefs, fixEffHr, fixEffYr, intercept, yr):
    tempVals = getTempVals(dataYr, tempCoefs)
    interactionVals = getInteractionVals(dataYr, intCoefs)
    interceptVals = np.repeat(np.array(intercept), dataYr.shape[0])
    fixEffHrVals = getFixEffHrVals(dataYr, fixEffHr)
    fixEffYrVals = getFixEffYrVal(yr, fixEffYr, dataYr)
    dataYr['load(MW)'] = tempVals + interactionVals + interceptVals + fixEffHrVals + fixEffYrVals


# Estimate temperature values by first putting T into bins (2d array), then multiplying by
# per-bin coefficients.
def getTempVals(dataYr, tempCoefs):
    temps = dataYr['tC']
    tempBins = createTempComponents(temps)
    return np.dot(tempBins, np.array(tempCoefs['value']))


# Element-wise of binned T by dew pt, then dot w/ coefficients
def getInteractionVals(dataYr, intCoefs):
    temps, rhs = dataYr['tC'], np.array(dataYr['rh'])
    tempBins = createTempComponents(temps)
    dewPt = convertRelHum2DewPoint(rhs, temps)
    return np.dot(np.multiply(dewPt, tempBins.T).T, np.array(intCoefs['value']))


# Add FE hour vals (by seaon, time of day, and type of day) to DF, then return FEs
def getFixEffHrVals(dataYr, fixEffHr):
    dataYrWithFE = dataYr.merge(fixEffHr, how='left', on=['season', 'type.day', 'hour.of.day'])
    return np.array(dataYrWithFE['value'])


# Return np array w/ FE for yr.
def getFixEffYrVal(yr, fixEffYr, dataYr):
    if yr in fixEffYr['year'].values:
        val = fixEffYr.loc[fixEffYr['year'] == yr, 'value']
    else:
        val = 0
    return np.repeat(np.array(val), dataYr.shape[0])

# forecastZonalDemandWithReg(2015,'pc',['S_C_TVA'],'C:\\Users\\mtcraig\\Desktop\\EPP Research\\PythonRIPSProject')
