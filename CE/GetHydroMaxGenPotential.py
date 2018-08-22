"""
Michael Craig
July 21, 2017
Import monthly hydropower generation from PNNL, and assign to each season considered in CE model.
"""

import os
from AuxFuncs import *
from GAMSAuxFuncs import createGenSymbol
from DemandFuncsCE import getHoursInMonths


def getHydroEPotential(fleet, demandZonal, repAndSpeHoursDict, currYear, dataRoot):
    """Computes maximum hydro potential for each hydro generator

    This function computes the maximum hydro potential for each hydro generator using PNNL data. Hydro potential is
    indexed by hydro unit and time block (special hours, summer, winter, ...). This max hydro potential will be used
    as upper bound in hourly hydro generation by the GAMS model. Scales potential

    July 2018: Currently we have only one hydro potential value per season. The objective is to decompose monthly
    potential to daily potential

    :param fleet: gen fleet
    :param demandZonal: zonal demand (dict)
    :param repAndSpeHoursDict: rep hrs per season or special block (dict of block:rep hrs, where seasons are labled
    by name and special hours are 'special')
    :param currYear: current year
    :param dataRoot: folder with input data
    :return: dict of season:genSymbol:generation potential
    """
    hydroPotentials = importHydroPotentialGen(currYear, dataRoot)
    hydroPotPerSeason = dict()
    plantCol, orisCol = fleet[0].index('PlantType'), fleet[0].index('ORIS Plant Code')
    zoneCol = fleet[0].index('Region Name')
    hydroUnits = [row for row in fleet[1:] if row[plantCol] == 'Hydro']
    for season in repAndSpeHoursDict:
        repHrs = repAndSpeHoursDict[season]
        months = getMonthsOfRepHrs(repHrs)
        monthlyDemandZonal = getMonthlyDemand(demandZonal, months)
        repHrsDemandZonal = getRepHrsDemand(demandZonal, repHrs)
        seasonDict, unitsNoData = dict(), list()

        for row in hydroUnits:
            oris, zone, genSymbol = row[orisCol], row[zoneCol], createGenSymbol(row, fleet[0])
            potential, hasData = getMonthsPotential(oris, hydroPotentials, months)

            if not hasData:
                unitsNoData.append(genSymbol)
            else:
                seasonDict[genSymbol] = potential * repHrsDemandZonal[zone] / monthlyDemandZonal[zone]

        if len(unitsNoData) > 0: assignPotentialsToMissingUnits(unitsNoData, fleet, seasonDict, repHrs)

        hydroPotPerSeason[season] = seasonDict

    return hydroPotPerSeason


def importHydroPotentialGen(currYear, dataRoot):
    """ Returns hydro potential gen for current year in 2d list w/ col 1 = plant IDs

    Read file with monthly hydro potential created by PNNL

    :param currYear: (integer)
    :param dataRoot: (string)
    :return: 2d list w/ col 1 = plant IDs
    """
    dataDir = os.path.join(dataRoot, 'HydroMonthlyDataPNNL')

    potentialsAllYears = readCSVto2dList(os.path.join(dataDir, 'monhydrogen_1550.csv'))
    # First row in data is months listed w/out years; assign years to each col
    startYr, endYr = [2015, 2050]
    yrIndex = (currYear - startYr) * 12 + 1  # +1 to account for first col = PlantID
    return [[row[0]] + row[yrIndex:yrIndex + 12] for row in potentialsAllYears]


def getMonthsOfRepHrs(repHrs):
    """ Function determines which months representative hours fall in, and outputs those months in a 1d list

    :param repHrs: 1d list of hours in a given season
    :return: 1d list of months
    """
    months = list()
    for month in range(1, 13):
        monthHours = getHoursInMonths([month])  # returns 1d list of hours in given month (1-8760 basis)
        if repHrs[0] in monthHours:
            months.append(month)
        elif repHrs[-1] in monthHours:
            months.append(month)
    return months


def getMonthlyDemand(demandZonal, months):
    """ Computes total accumulated demand (in MWh) in given months

    :param demandZonal: 1d list of hourly demand
    :param months: 1d list of months
    :return: dict of zone:total demand in those months
    """
    monthsHours = getHoursInMonths(months)
    demandInMonthsZonal = dict()
    for zone in demandZonal: demandInMonthsZonal[zone] = sumDemandInHours(demandZonal[zone], monthsHours)
    return demandInMonthsZonal


def getRepHrsDemand(demandZonal, repHrs):
    """Computes total accumulated demand (in MWh) in given subset of representative hours

    :param demandZonal:  1d list of demand zonal (dict of zone:hourly demand),
    :param repHrs: 1d list of rep hrs per season
    :return: dict of zone:total demand in those rep hrs
    """
    demandInRepHrsZonal = dict()
    for zone in demandZonal: demandInRepHrsZonal[zone] = sumDemandInHours(demandZonal[zone], repHrs)
    return demandInRepHrsZonal


def sumDemandInHours(demand, hours):
    """

    :param demand: 1d list of demand (0-8759 idx)
    :param hours: 1d list of hours (1-8760 idx).
    :return: sum of demand in those hours
    """
    return sum([demand[hr - 1] for hr in hours])


def getMonthsPotential(oris, hydroPotentials, months):
    """Gets total hydropower generation potential for month(s).

    :param oris: ORIS ID
    :param hydroPotentials: hydro potential generation (2d list) for curr year
    :param months: months of interest
    :return: total potential generation for month(s)
    """
    orisCol = hydroPotentials[0].index('PlantID')
    potentialOriss = [row[orisCol] for row in hydroPotentials]
    monthsPotentials, hasData = list(), True
    if oris in potentialOriss:
        unitRowIdx = potentialOriss.index(oris)
        allMonthsPotentials = hydroPotentials[unitRowIdx]  # leave ORIS label on so idx & month (1-12) align
        monthsPotentials = [float(allMonthsPotentials[month]) if float(allMonthsPotentials[month]) > 0 else 0 for month
                            in months]
    else:
        hasData = False
    return sum(monthsPotentials), hasData


def assignPotentialsToMissingUnits(unitsNoData, fleet, seasonDict, repHrs):
    """For ORIS units not in PNNL data, assign them potential for season based on
    average potential as CF of rest of fleet for that season. Modifies seasonDict to include other units

    :param unitsNoData: gen symbols (oris+unit ID) for units w/out data, ,
    :param fleet: 2d list with gen fleet
    :param seasonDict: dict of season:genSymbol:potential
    :param repHrs: 1d list of rep hrs per season
    """
    plantCol = fleet[0].index('PlantType')
    capacCol = fleet[0].index('Capacity (MW)')
    hydroUnits = [row for row in fleet[1:] if row[plantCol] == 'Hydro']
    genSymbols = [createGenSymbol(row, fleet[0]) for row in hydroUnits]
    capacs = [row[capacCol] for row in hydroUnits]
    fleetAvgCF = getFleetAverageCF(genSymbols, capacs, seasonDict, repHrs)
    for genSymbol in unitsNoData:
        seasonDict[genSymbol] = capacs[genSymbols.index(genSymbol)] * len(repHrs) * fleetAvgCF


def getFleetAverageCF(genSymbols, capacs, seasonDict, repHrs):
    """Gets average fleet CF for rep hours in given season.

    :param genSymbols: 1d list of all hydro gen symbols
    :param capacs: 1d list of all hydro capac,
    :param seasonDict: dict of season:gensymbol:potential,
    :param repHrs: rep hours (1d list)
    :return: average fleet CF in season
    """
    cfs = list()
    for (genSymbol, capac) in zip(genSymbols, capacs):
        if genSymbol in seasonDict:
            potential = seasonDict[genSymbol]
            cfs.append(potential / (capac * len(repHrs)))
    return sum(cfs) / len(cfs)


def getListOfDaysInMonth(m):
    """ Gets list with days of the year in given month m

    :param m: month (1-12 basis)
    :return: list of days of the year in given month (0-364 basis)
    """

    daysPerMonth = getDaysPerMonth()

    if m == 1:
        initialDay = 0
    else:
        initialDay = sum([daysPerMonth[month-1][1] for month in range(1, 13) if month < m])

    listDays = [initialDay + d for d in range(daysPerMonth[m-1][1])]

    return listDays


def compute_max_daily_hydro(fleet, currYear, dataRoot):
    """ Converts monthly hydro capacity for each plant to daily hydro capacity

    This conversion follows a simple linear rule with simulated water releases simulated by PNNL

    $$P^{MAX}_{day} = P_{month}^{PNNL}\frac{r_{day}^{PNNL}}{\sum_{day \in month} r_{day}^{PNNL}}$$

    :param fleet: gen fleet
    :param currYear: current year
    :param dataRoot: folder with input data
    :return: dictionary with daily hydro capacity {gensymbol: 1d list of daily capacities}
    """

    # reads monthly hydro potential for year currYear (in MWh)
    hydroPotentials = importHydroPotentialGen(currYear, dataRoot)

    # reads daily water releases for year currYear
    daily_releases = importHydroDailyReleases(currYear, dataRoot)

    plantCol, orisCol = fleet[0].index('PlantType'), fleet[0].index('ORIS Plant Code')
    zoneCol = fleet[0].index('Region Name')
    hydroUnits = [row for row in fleet[1:] if row[plantCol] == 'Hydro']

    dailyPotentialDict = dict()

    for row in hydroUnits:

        oris, zone, genSymbol = row[orisCol], row[zoneCol], createGenSymbol(row, fleet[0])
        daily_releases_unit, hasData = getReleasesForUnit(oris, daily_releases)

        dailyPotential = list()

        # loop through days of year (no leap years)
        for m in range(1, 13):
            potMonth, hasData = getMonthsPotential(oris, hydroPotentials, [m])

            daysInMonth = getListOfDaysInMonth(m)

            monthlyReleasesTot = sum([daily_releases_unit[d] for d in daysInMonth])

            dailyPotential = dailyPotential + [potMonth * (daily_releases_unit[d]/monthlyReleasesTot) for d in daysInMonth]

        dailyPotentialDict[genSymbol] = dailyPotential

    return dailyPotentialDict


def importHydroDailyReleases(currYear, dataRoot):
    """ Returns hydro daily releases for current year in 2d list w/ col 1 = plant IDs

    Read file with daily hydro releases created by PNNL (Assumes format is the same as the monthly file)

    :param currYear: (integer)
    :param dataRoot: (string)
    :return: 2d list w/ col 1 = plant IDs
    """
    dataDir = os.path.join(dataRoot, 'HydroMonthlyDataPNNL')

    releasesAllYears = readCSVto2dList(os.path.join(dataDir, 'monhydrogen_1550.csv'))

    # First row in data is months listed w/out years; assign years to each col
    startYr, endYr = [2015, 2050]

    yrIndex = (currYear - startYr) * 365 + 1  # +1 to account for first col = PlantID

    return [[row[0]] + row[yrIndex:yrIndex + 365] for row in releasesAllYears]


def getReleasesForUnit(oris, daily_releases):
    """Gets total hydropower generation potential for month(s).

    :param oris: ORIS ID
    :param daily_releases: hydro daily releases (2d list) for curr year
    :return: 1d list with daily release for hydro unit in curr year (index of first year is 0)
    """

    orisCol = daily_releases[0].index('PlantID')

    releasesOrisIDs = [row[orisCol] for row in daily_releases]

    daily_releases_unit, hasData = list(), True

    if oris in releasesOrisIDs:
        unitRowIdx = releasesOrisIDs.index(oris)
        daily_releases_unit = daily_releases[unitRowIdx]

        daily_releases_unit = daily_releases_unit[1:]  # remove label in first column (days are in 0-364 basis)

        daily_releases_unit = [float(v) for v in daily_releases_unit] # convert to numeric

    else:
        hasData = False

    return daily_releases_unit, hasData
