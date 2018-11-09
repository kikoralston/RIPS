# Michael Craig
# October 4, 2016
# Functions related to electricity demand profile

import operator, os
from AuxFuncs import *


def scaleDemandForGrowthAndEE(baseDemand, annualDemandGrowth, demandYear, currYear):
    """ GET SCALED DEMAND BASED ON CURRENT YEAR

    This function is not used anymore

    :param baseDemand: initial demand values (1d list)
    :param annualDemandGrowth: annual demand growth (fraction)
    :param demandYear: year of initial demand values
    :param currYear: current CE year
    :return: demand values in current CE year (1d list)
    """
    demandScalar = (1 + annualDemandGrowth) ** (currYear - demandYear)
    return [val * demandScalar for val in baseDemand]


def getAggregateSolarAndWind(windCFs, ewdIdAndCapac, solarCFs, solarFilenameAndCapac):
    """ Aggregates existing solar and wind generation by zone

    These aggregated values will be used to compute a net demand value

    :param windCFs: wind CFs (2d list OR None if no wind/solar in zone)
    :param ewdIdAndCapac: list of solar/wind IDs and their capacities in fleet (2d list)
    :param solarCFs: solar CFs (2d list OR None if no wind/solar in zone)
    :param solarFilenameAndCapac:
    :return: hourly wind gen (1d list w/out headers)
             hourly solar gen (1d list w/out headers)
    """

    if windCFs is None:
        hourlyWindGen = [0 for val in range(len(hourlyDemand))]
    else:
        hourlyWindGen = getHourlyGenProfile(windCFs, ewdIdAndCapac)

    if solarCFs is None:
        hourlySolarGen = [0 for val in range(len(hourlyDemand))]
    else:
        hourlySolarGen = getHourlyGenProfile(solarCFs, solarFilenameAndCapac)

    return hourlyWindGen, hourlySolarGen


def getNetDemand(hourlyDemand, hourlyWindGen, hourlySolarGen):
    """ GET NET DEMAND AND REMOVE WIND & SOLAR FROM FLEET

    :param hourlyDemand: hourly demand values in zone (1d list w/out header)
    :param hourlyWindGen: hourly wind gen in zone (1d list w/out headers)
    :param hourlySolarGen: hourly solar gen in zone (1d list w/out headers)
    :return: net demand in zone (1d list w/out headers),
    """

    if len(hourlyWindGen) > 0 and len(hourlySolarGen) > 0:
        hourlyWindAndSolarGen = list(map(operator.add, hourlyWindGen, hourlySolarGen))
        return list(map(operator.sub, hourlyDemand, hourlyWindAndSolarGen))
    elif len(hourlyWindGen) > 0:
        return list(map(operator.sub, hourlyDemand, hourlyWindGen))
    elif len(hourlySolarGen) > 0:
        return list(map(operator.sub, hourlyDemand, hourlySolarGen))
    else:
        return hourlyDemand


def getHourlyGenProfile(cfs, idAndCapacs):
    """Computes aggregate Solar or Wind hourly generation profile

    :param cfs: CFs (2d list w/ header)
    :param idAndCapacs: unit ids and capacities (2d list w/ header)
    :return: hourly generation values (1d list w/out header)
    """
    (idCol, capacCol) = (idAndCapacs[0].index('Id'), idAndCapacs[0].index('FleetCapacity'))
    totalHourlyGen = []
    for idAndCapac in idAndCapacs[1:]:
        (unitID, capac) = (idAndCapac[idCol], idAndCapac[capacCol])
        cfRow = [row[1:] for row in cfs[1:] if row[0] == unitID]
        hourlyGen = [capac * cf for cf in cfRow[0]]
        if totalHourlyGen == []:
            totalHourlyGen = hourlyGen
        else:
            for hr in range(len(hourlyGen)): totalHourlyGen[hr] += hourlyGen[hr]
    return totalHourlyGen