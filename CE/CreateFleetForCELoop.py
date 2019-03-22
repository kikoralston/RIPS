"""
Michael Craig
October 4, 2016
Create fleet for current CE loop by removing retired units from fleet.
Determines which units retire due to age, and also accounts for past retirements
for age and economic reasons.
"""

import os, copy
from AuxFuncs import *
from GAMSAuxFuncs import createGenSymbol


def createFleetForCurrentCELoop(genFleet, currYear, capacExpRetiredUnitsByAge, dataRoot, scenario):
    """CREATE FLEET FOR CURRENT CE LOOP

    :param genFleet: gen fleet (2d list)
    :param currYear: current CE year
    :param capacExpRetiredUnitsByAge: list of units retired each year (2d list)
    :param dataRoot:
    :param scenario:
    :return: gen fleet w/ retired units removed (2d list)
    """
    markAndSaveRetiredUnitsFromAge(genFleet, currYear, capacExpRetiredUnitsByAge, dataRoot, scenario)
    genFleetForCE = [genFleet[0]] + [row for row in genFleet[1:] if onlineAndNotRetired(row, genFleet[0], currYear)]
    return genFleetForCE


def markAndSaveRetiredUnitsFromAge(genFleet, currYear, capacExpRetiredUnitsByAge, dataRoot, scenario):
    """RETIRE UNITS BY AGE

    Marks units that retire in gen fleet list due to lifetime and modifies them in list.
    The function returns nothing, but modifies list genFleet.

    :param genFleet: gen fleet (2d list)
    :param currYear: current CE year
    :param capacExpRetiredUnitsByAge:
    :param dataRoot: list of units retired each year (2d list)
    :param scenario:
    """
    lifetimeByPlantTypeDict = importPlantTypeLifetimes(dataRoot, scenario)
    renewablePlantTypes = ['Geothermal', 'Hydro', 'Pumped Storage', 'Wind', 'Solar PV']
    onlineYearCol = genFleet[0].index('On Line Year')
    plantTypeCol = genFleet[0].index('PlantType')
    retiredByAgeCol = genFleet[0].index('YearRetiredByAge')
    retiredByCECol = genFleet[0].index('YearRetiredByCE')
    retiredUnitsByAge = []

    for row in genFleet[1:]:
        if row[retiredByAgeCol] == '' and row[retiredByCECol] == '':  # if not already retired by age or CE
            (onlineYear, plantType) = (row[onlineYearCol], row[plantTypeCol])
            lifetimePlantType = lifetimeByPlantTypeDict[plantType]
            if int(onlineYear) + lifetimePlantType < currYear:
                if plantType in renewablePlantTypes:
                    readdRenewablePlant(genFleet, row, currYear)  # readd to fleet before add retired year!

                row[retiredByAgeCol] = currYear
                retiredUnitsByAge.append(createGenSymbol(row, genFleet[0]))

    capacExpRetiredUnitsByAge.append(['UnitsRetiredByAge' + str(currYear)] + retiredUnitsByAge)


def importPlantTypeLifetimes(dataRoot, scenario):
    """Import lifetimes for each plant type

    This function imports lifetime data for each plant type

    :param dataRoot: string with complete path to root of data folder
    :param scenario: string with type of scenario being simulated
    :return: dictionary with lifetime for each plant type
    """
    lifetimeDir = os.path.join(dataRoot, 'NewPlantData')

    if scenario == 'coalret':
        lifetimeFilename = 'LifetimeValuesExistingPlants4Aug2016EarlyCoalRet.csv'
    else:
        lifetimeFilename = 'LifetimeValuesExistingPlants4Aug2016.csv'

    lifetimeData = readCSVto2dList(os.path.join(lifetimeDir, lifetimeFilename))
    lifetimeByPlantTypeDict = convert2dListToDictionaryWithIntVals(lifetimeData,
                                                                   'PlantType', 'Lifetime(yrs)')
    return lifetimeByPlantTypeDict


def convert2dListToDictionaryWithIntVals(list2d, keyHeader, valHeader):
    """Converts Zx2 2d list to dictionary

    Generic utility function that converts a Zx2 list to dictionary

    :param list2d: 2d list (2 cols)
    :param keyHeader: header of col w/ keys
    :param valHeader: header of col w/ vals
    :return: dictionary
    """
    dictResult = dict()
    (keyCol, valCol) = (list2d[0].index(keyHeader), list2d[0].index(valHeader))
    for row in list2d[1:]: dictResult[row[keyCol]] = int(row[valCol])
    return dictResult


def readdRenewablePlant(genFleet, row, currYear):
    """READ RENEWABLES THAT RETIRE DUE TO AGE

    If renewable retires due to age, automatically adds new unit to end of fleet
    that is same as old unit except for unit ID & online year.

    :param genFleet: gen fleet (2d list)
    :param row: row in genFleet of retired RE unit
    :param currYear: current CE year
    """
    (unitIdCol, onlineCol) = (genFleet[0].index('Unit ID'), genFleet[0].index('On Line Year'))
    newRow = copy.deepcopy(row)
    newRow[unitIdCol] += 'Replaced'
    newRow[onlineCol] = currYear
    genFleet.append(newRow)


def onlineAndNotRetired(genRow, headers, currYear):
    """REMOVE UNITS RETIRED FROM FLEET

    Checks if unit has already gone online and is not retired

    :param genRow: row of gen fleet (1d list)
    :param headers: headers of gen fleet (1d list)
    :param currYear: current year
    :return: True if online & not retired, False otherwise
    """
    (ipmRetiredCol, ceRetiredCol) = (headers.index('Retirement Year'), headers.index('YearRetiredByCE'))
    retiredByAgeCol = headers.index('YearRetiredByAge')
    onlineCol = headers.index('On Line Year')
    if int(genRow[onlineCol]) > currYear:
        return False  # some units don't come online until 2020
    elif int(genRow[ipmRetiredCol]) < currYear:
        return False  # units flagged as retiring by IPM
    elif genRow[retiredByAgeCol] != '' and int(genRow[retiredByAgeCol]) <= currYear:
        return False  # units retired due to age
    elif genRow[ceRetiredCol] != '' and int(genRow[ceRetiredCol]) <= currYear:
        return False  # units retired by CE
    else:
        return True
