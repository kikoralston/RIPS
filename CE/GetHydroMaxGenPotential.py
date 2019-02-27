"""
Michael Craig
July 21, 2017
Import monthly hydropower generation from PNNL, and assign to each season considered in CE model.
"""

import os
from AuxFuncs import *
from GAMSAuxFuncs import createGenSymbol
from DemandFuncsCE import getHoursInMonths
import pandas as pd
import pickle as pk
import datetime as dt
import calendar
import tarfile
import shutil


def getHydroEPotential(fleet, demandZonal, repAndSpeHoursDict, currYear, genparam, curtailparam):
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
    :param genparam: object of class Generalparameters
    :param curtailparam: object of class Curtailmentparameters

    :return: nested dict of {gcm:{season:{genSymbol:generation potential}}}
    """

    hydroPotentials = {gcm: importHydroPotentialGen(currYear, genparam, gcm) for gcm in curtailparam.listgcms}

    hydroPotentialsTotal = dict()

    for gcm in repAndSpeHoursDict:
        hydroPotPerSeason = dict()
        plantCol, orisCol = fleet[0].index('PlantType'), fleet[0].index('ORIS Plant Code')
        zoneCol = fleet[0].index('Region Name')
        hydroUnits = [row for row in fleet[1:] if row[plantCol] == 'Hydro']
        for season in repAndSpeHoursDict[gcm]:
            repHrs = repAndSpeHoursDict[gcm][season]
            months = getMonthsOfRepHrs(repHrs)
            monthlyDemandZonal = getMonthlyDemand(demandZonal[gcm], months)
            repHrsDemandZonal = getRepHrsDemand(demandZonal[gcm], repHrs)
            seasonDict, unitsNoData = dict(), list()

            for row in hydroUnits:
                oris, zone, genSymbol = row[orisCol], row[zoneCol], createGenSymbol(row, fleet[0])
                potential, hasData = getMonthsPotential(oris, hydroPotentials[gcm], months)

                if not hasData:
                    unitsNoData.append(genSymbol)
                else:
                    seasonDict[genSymbol] = potential * repHrsDemandZonal[zone] / monthlyDemandZonal[zone]

            if len(unitsNoData) > 0: assignPotentialsToMissingUnits(unitsNoData, fleet, seasonDict, repHrs)

            hydroPotPerSeason[season] = seasonDict

        hydroPotentialsTotal[gcm] = hydroPotPerSeason

    return hydroPotentialsTotal


def importHydroPotentialGen(currYear, genparam, gcm, format='new'):
    """ Returns hydro potential gen for current year in 2d list w/ col 1 = plant IDs

    Read file with monthly hydro potential created by PNNL

    :param currYear: (integer)
    :param genparam: object of class Generalparameters
    :param gcm: (string) name of gcm model
    :param format: (string) 'new' refers to the new hydro potential file format from october 2018

    :return: 2d list w/ col 1 = plant IDs
    """
    dataDir = os.path.join(genparam.dataRoot, 'HydroMonthlyDataPNNL')

    if genparam.referenceCase:
        # get historical values

        a = pd.read_csv(os.path.join(dataDir, 'PlantLocationEIA2003_2016.txt'), sep='\t')
        a = a[list(a.columns[:17])]
        a = a.loc[~a['OrisId'].isna(), :]

        a = a[['OrisId', '1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12']]
        a = a.rename(columns={'OrisId': 'PlantID'})
        a['PlantID'] = a['PlantID'].astype(int)
        a = a.astype(str)

        out = [a.columns.tolist()] + a.values.tolist()

    else:
        name_gcm = gcm

        # reading flow file. substitute '_' to '.' and change 'RCP' to small caps
        name_gcm = name_gcm.replace('_', '.')
        name_gcm = name_gcm.replace('rcp', 'RCP')

        if format == 'new':
            with open(os.path.join(dataDir, 'monthlyhydropotential.pk'), 'rb') as f:
                allData = pk.load(f)

            potentialsAllYears = allData[name_gcm]

            potentialsCurrYear = potentialsAllYears[potentialsAllYears['months'].dt.year == currYear]
            potentialsCurrYear = potentialsCurrYear.drop(columns=['months'])
            potentialsCurrYear = potentialsCurrYear.reset_index(drop=True)
            potentialsCurrYear = potentialsCurrYear.T
            potentialsCurrYear = potentialsCurrYear.reset_index()

            out = potentialsCurrYear.astype(str)
            out = out.values.tolist()
            out = [['PlantID', '1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12']] + out

        else:
            potentialsAllYears = readCSVto2dList(os.path.join(dataDir, 'monhydrogen_1550.csv'))
            # First row in data is months listed w/out years; assign years to each col
            startYr, endYr = [2015, 2050]
            yrIndex = (currYear - startYr) * 12 + 1  # +1 to account for first col = PlantID
            out = [[row[0]] + row[yrIndex:yrIndex + 12] for row in potentialsAllYears]

    return out


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


def getDailyHydroPotentialsUC(fleetUC, hydroData, daysUC, ucYear):
    """ Compiles dict with daily hydro potentials

    Compiles dict with daily hydro potentials for individual hydro gens for UC run from list with monthly potentials

    :param fleetUC: 2d list with complete generator fleet
    :param hydroData: 2d list with monthly hydro potentials for simulation year
    :param daysUC: list with days of the year (1-365) that UC simulation will be performed
    :return: dict with {genID: hydro energy MWh}
    """

    headers = fleetUC[0]
    plantCol = headers.index('PlantType')
    orisCol = headers.index('ORIS Plant Code')
    capacityCol = headers.index('Capacity (MW)')
    hydroFleet = [headers] + [row for row in fleetUC if row[plantCol] == 'Hydro']

    hydroDailyPotential = dict()

    for rowFleet in hydroFleet[1:]:

        genSymbol = createGenSymbol(rowFleet, headers)

        # get row in hydroFleet of this generator
        rowHydroData = [row for row in hydroData[1:] if row[0] == rowFleet[orisCol]]

        if len(rowHydroData) > 0:
            rowHydroData = rowHydroData[0]

            capacity = float(rowFleet[capacityCol])

            hydroEnergyAux = 0
            for d in daysUC:
                date_aux = dt.date(ucYear, 1, 1) + dt.timedelta(days=d-1)
                imonth = date_aux.month
                ndays = calendar.monthrange(ucYear, imonth)[1]
                hydroEnergyAux = hydroEnergyAux + float(rowHydroData[imonth])/ndays

            hydroDailyPotential[genSymbol] = hydroEnergyAux
            #hydroDailyPotential[genSymbol] = hydroEnergyAux/(capacity*len(daysUC)*24)

        else:
            hydroDailyPotential[genSymbol] = 0

    return hydroDailyPotential


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


def importHydroDailyReleases(currYear, dataRoot, gcm):
    """ Returns hydro daily releases for current year in 2d list w/ col 1 = plant IDs

    Read file with daily hydro releases created by PNNL (Assumes format is the same as the monthly file)

    :param currYear: (integer)
    :param dataRoot: (string)
    :return: 2d list w/ col 1 = plant IDs
    """
    dataDir = os.path.join(dataRoot, 'HydroMonthlyDataPNNL')

    with open(os.path.join(dataDir, 'daily_releases_{}.pk'.format(gcm)), 'rb') as f:
        releasesAllYears = pk.load(f)


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


def processHydroPotential(pathin, pathout):
    """
    This routine processes the raw hydro potential data created by Nathalie Voisin (PNNL)

    October 2018

    This routine reads the txt files created by Nathalie Voisin (PNNL) for each hydro generator in SERC and processes
    it in order to create the monthly hydro potential values (in MWh) used by the CE model

    :return:
    """

    x = pd.read_table(os.path.join(pathin, 'PlantLocationEIA2003_2016.txt'),
                      usecols=['OrisId', 'Latitude', 'Longitude', 'Lt_final', 'LN_final',
                               '1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12'])

    # remove rows with ID NaN
    x = x[~(pd.isnull(x['OrisId']))]

    x['OrisId'] = x['OrisId'].astype(int)

    totaldata = {}

    for i, id in enumerate(x['OrisId']):
        labels = []
        values = []

        print('{0:4d}. Reading hydro gen {1}'.format(i, id))

        with open(os.path.join(pathin, '{0:d}.txt'.format(id))) as f:
            for i, line in enumerate(f):

                # remove break line character
                if line[-1] == '\n':
                    line = line[:-1]

                line = line.strip()

                if i % 2 == 0:
                    labels.append(line)
                else:
                    values.append(list(map(float, line.split(sep=' '))))

            labels = labels[1:]
            values = values[1:]

        out = pd.DataFrame()
        # expand to monthly values
        months = pd.date_range(start='2012-01-01', end='2100-12-01', freq='MS')

        out['months'] = months

        for i, l in enumerate(labels):
            gcm = l.split(sep='.regulated')[0]
            gcm = gcm.split(sep='serc.')[1]

            y = sum([list((x[x['OrisId'] == id].iloc[0, 5:17]).values*z) for z in values[i]], [])
            out[gcm] = y

        totaldata[id] = out

    # rearrange totaldata dictionary as {gcm: dataframe with monthly hydro potential}
    # data frame will have each column as one hydro generator and row as months

    print('Rearranging dictionary to have GCMs as keys...')

    totaldata2 = {}

    # get list of gcms from last data frame
    list_gcms = totaldata[id].columns
    list_gcms = list_gcms[1:]

    for g in list_gcms:
        out = pd.DataFrame()
        for id in x['OrisId']:
            out['months'] = months
            out[id] = totaldata[id][g]

        totaldata2[g] = out

    print('Saving final dictionary as a pickle file...')

    # save final dictionary to pickle file
    with open(os.path.join(pathout, 'monthlyhydropotential.pk'), 'wb') as f:
        pk.dump(totaldata2, f)

    print('Done!')

#    with open(os.path.join(pathout, 'monthlyhydropotential.pk'), 'rb') as f:
#        x = pk.load(f)


def processDailyReleases(pathin, pathout):
    """
    This routine processes the raw hydro release data created by Nathalie Voisin (PNNL) and Yifan Cheng (UW)

    February 2019

    :return:
    """

    # create tmp folder for processing data
    os.mkdir(os.path.join(pathout, 'tmp'))

    list_tarfiles = [f for f in os.listdir(pathin) if f[-4:] == '.tgz']

    for i, f in enumerate(list_tarfiles):

        print('{0:4d}. Processing daily releases file {1}'.format(i, f))

        tar = tarfile.open(os.path.join(pathin, f), "r:gz")

        # get list with names of txt files in tar
        listfiles = tar.getnames()

        # extract all files and close tar
        tar.extractall(path=os.path.join(pathout, 'tmp'))
        tar.close()

        df_total = None

        for f in listfiles:
            fname = os.path.join(pathout, 'tmp', f)
            df = pd.read_csv(fname)

            gcm_rcp = df['gcm_rcp'].iloc[0]
            del df['gcm_rcp']

            # convert string to datetime
            df['time'] = pd.to_datetime(df['time'], format='%Y-%m-%d')

            if df_total is None:
                df_total = pd.DataFrame(df)
            else:
                df_total = df_total.append(df, ignore_index=True)

        # save df_total to pickle file
        df_total.to_pickle(os.path.join(pathout, 'daily_releases_{}.pk'.format(gcm_rcp)))

        # delete extracted data
        shutil.rmtree(os.path.dirname(os.path.join(pathout, 'tmp', listfiles[0])))

    os.rmdir(os.path.join(pathout, 'tmp'))








