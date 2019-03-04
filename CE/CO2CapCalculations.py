# Michael Craig
# October 4, 2016
# Project CPP CO2 emissions cap to given year (forward or backward) based on 2022
# & 2030 limits.

from AuxFuncs import *
import os
import json


def readInitialCondCO2(fname):
    """ Reads file with initial conditions of CO2 emission

    File should be a json with the following format

    {
       "startYr": "YYYY",
       "startEms": "XXXX"
    }

    or

    [
       {
          "startYr": "YYYY",
          "startEms": "XXXX"
       }
    ]

    :param fname: string with complete path to file
    :return: tuple with Start year and total emissions in start year (in short tons)
    """

    if os.path.exists(fname):
        with open(fname) as f:
            data = json.load(f)

        # check if it is a list
        if isinstance(data, list):
            data = data[0]

        startYr = int(data['startYr'])
        startEms = float(data['startEms'])
    else:

        startYr = 2015
        #startEms = 251929000
        startEms = 380000000

        print()
        print('File {0} not found!'.format(fname))
        print('Setting default values for CO2 emissions parameters!')
        print('Start Year: {0:4d}'.format(startYr))
        print('Start Emission: {0:,d} short tons CO2'.format(startEms))
        print()

    return startYr, startEms


def getCo2Cap(co2CapScenario, startEms=380000000):
    """
    Get CO2 cap for given scenario. Refs: For emissions limits, see: see Databases, CO2EmissionERCOT,
    UCBaseCase2015Output9April2017 folder, baseCaseCo2Emissions9April2017.xlsx.

    :param co2CapScenario: name of CO2 scenario
    :param startEms: total emissions (in short tons) in initial year
    :return:
    """
    if co2CapScenario == 'cpp':
        capYear, capEms = 2050, 0.5*startEms  # 50% redux
    elif co2CapScenario == 'deep':
        capYear, capEms = 2050, 0.8*startEms  # 80% redux
    elif co2CapScenario == 'none':
        capYear, capEms = 2050, float('inf')  # set to arbitrarily large value so effectively no cap

    return capYear, capEms  # short tons!


def interpolateCO2Cap(currYear, genparam):
    """Computes a CO2 cap for current year

    Data soruce for 2015 emissions: UC run w/out co2 price or storage. See xls file above.
    (data from EIA, "Texas electricity profile 2017", Table 7, electric power industry emissions estimates, 1990-2014).

    :param currYear: current year
    :param genparam: object of class Generalparameters
    :return: projected co2 emissions cap for current year (short tons)
    """
    startYr, startEms = readInitialCondCO2(fname=os.path.join(genparam.dataRoot, 'co2values.json'))

    endYr, endLimit = getCo2Cap(genparam.co2CapScenario, startEms)

    #if currYear == startYr: startEms *= 10  # if first year, don't want to enforce co2 cap, so just scale up

    (deltaYear, deltaEms) = (endYr - startYr, startEms - endLimit)
    emsReduxPerYear = deltaEms / deltaYear
    diffCurrYearFromStart = currYear - startYr

    co2Cap = (startEms) if diffCurrYearFromStart == 0 else (startEms - diffCurrYearFromStart * emsReduxPerYear)

    return co2Cap


################################################################################


##################### OLD FUNCTIONS #############################################
# Inputs: year to project cap to, 2022 and 2030 CPP CO2 emissions cap (short tons)
# Outputs: projected CO2 emissions cap (short tons)
def setCppCapBasedOnCurrYear(currYear, co2Cpp2022Limit, co2Cpp2030Limit):
    (deltaYear, deltaEms) = (2030 - 2022, co2Cpp2022Limit - co2Cpp2030Limit)
    emsReduxPerYear = deltaEms / deltaYear
    diffCurrYearFrom2022 = currYear - 2022
    return co2Cpp2022Limit - diffCurrYearFrom2022 * emsReduxPerYear


# Import 2022 & 2030 regional CPP mass limit for all states in region
# Output emissions are in short tons/yr, and include new source complement.
def calcRegionCPPLimits(states, dataRoot):

    cppLimitDir = os.path.join(dataRoot, 'Clean Power Plan')

    cppLimitName = 'StateCPPCO2EmsCaps5Oct2016.csv'

    stateCppLimits = readCSVto2dList(os.path.join(cppLimitDir, cppLimitName))
    stateCol = stateCppLimits[0].index('State')
    startYr, endYr = 2022, 2030
    (limit2022Col, limit2030Col) = (stateCppLimits[0].index(str(startYr)), stateCppLimits[0].index(str(endYr)))
    state2022Limits = [float(row[limit2022Col]) for row in stateCppLimits[1:] if row[stateCol] in states]
    state2030Limits = [float(row[limit2030Col]) for row in stateCppLimits[1:] if row[stateCol] in states]
    return (startYr, endYr, sum(state2022Limits), sum(state2030Limits))  # short tons/yr; includes new source complement


# Output emissions are in short tons/yr.
# Data source: see 'C:\Users\mtcraig\Desktop\EPP Research\Databases\DeepDecarbCO2LimitsTX'
def calcDeepDecarbLimits():
    startYr, endYr = 2014, 2050
    startEms, endEms = 280582138, 46980374
    return startYr, endYr, startEms, endEms
