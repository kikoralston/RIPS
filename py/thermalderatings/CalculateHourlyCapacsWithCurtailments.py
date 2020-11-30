# Michael Craig, 17 May 2017
# Calculate hourly capacity including curtailments for existing generators.
# Outputs dict of gen:hourly capacities and a 2d list w/ generator ID followed by curtailments.

from GAMSUtil.GAMSAuxFuncs import *
import copy


def calculateHourlyCapacsWithCurtailments(genFleet, hrlyCurtailmentsAllGensInTgtYr, currYear):
    """Compute hourly available capacity in MW after deratings for existing fleet

    This function uses the capacity deratings in % of capacity computed by function :fun:`determineHrlyCurtailmentsForExistingGens`
    and computes time series of available capacity (in MW) for all generators in the fleet.

    :param genFleet: 2d list with Generation Fleet data matrix
    :param hrlyCurtailmentsAllGensInTgtYr: dict mapping each gen to 2d list of datetime for year of run to hourly net
                                           capacity curtailments (fraction of total capacity)
    :param currYear: integer representing year being simulated
    :return: (dict) a dictionary {gcm: {gen symbol: [hourly available capacity in MW]}}
    """
    genHourlyCapacs = dict()

    fuelTypeCol, capacCol = genFleet[0].index('Modeled Fuels'), genFleet[0].index('Capacity (MW)')

    for gcm in hrlyCurtailmentsAllGensInTgtYr:
        auxDict = dict()

        if len(hrlyCurtailmentsAllGensInTgtYr[gcm]) > 0:
            lenCurtailments = len(hrlyCurtailmentsAllGensInTgtYr[gcm][list(hrlyCurtailmentsAllGensInTgtYr[gcm].keys())[0]])
        else:
            lenCurtailments = 8760

        for row in genFleet[1:]:

            genSymbol = createGenSymbol(row, genFleet[0])

            if row[fuelTypeCol] == 'Wind' or row[fuelTypeCol] == 'Solar':
                hourlyCapacs = [float(row[capacCol])] * lenCurtailments
            elif genSymbol in hrlyCurtailmentsAllGensInTgtYr[gcm]:
                hrlyCurtailments = copy.deepcopy(hrlyCurtailmentsAllGensInTgtYr[gcm][genSymbol])
                hourlyCapacs = subtractCurtailmentsFromCapac(hrlyCurtailments, float(row[capacCol]), genSymbol)
            else:  # no curtailments
                hourlyCapacs = [float(row[capacCol])] * lenCurtailments  # has header

            auxDict[genSymbol] = hourlyCapacs

        genHourlyCapacs[gcm] = auxDict

    return genHourlyCapacs


def subtractCurtailmentsFromCapac(hrlyCurtailments, capac, genSymbol):
    """Computes net available capacity (in MW) of individual generator

    :param hrlyCurtailments: list with hourly curtailments for generator in current year (fraction of total capacity)
    :param capac: capacity with generator (MW)
    :param genSymbol: string with generator symbol (ORIS+UNIT)
    :return: list of Hourly curtailments in MW
    """
    hourlyCapacs = [float(val) * capac for val in hrlyCurtailments]

    if min(hourlyCapacs) < 0:
        print('Hourly capacity < 0 for ', genSymbol, ' at idx ', hourlyCapacs.index(min(hourlyCapacs)))
        hourlyCapacs = [max(0, val) for val in hrlyCurtailments]

    return hourlyCapacs


def calculateHourlyTechCapacsWithCurtailments(newTechsCE, hrlyCurtailmentsAllTechsInTgtYr, currYear, ptCurtailedAll):
    """Compute hourly available capacity of new thermal generators

    This function uses the capacity deratings in % of capacity computed by function :fun:`determineHrlyCurtailmentsForNewTechs`
    and computes time series of available capacity (in MW) for all candidate technologies.

    :param newTechsCE: (2d list) list with data for new technologies
    :param hrlyCurtailmentsAllTechsInTgtYr: (dict) a dictionary {gcm: {(techSymbol, cell): [hourly deratings as fraction of capacity]}}
    :param currYear: (int) current year
    :param ptCurtailedAll: (list) list with names of types of plants that are curtailed
    :return: (dict) a dictionary {gcm: {(techSymbol, cell): [hourly available capacity in MW]}}
    """
    thermalTechHourlyCapacs = dict()
    (fuelTypeCol, capacCol) = (newTechsCE[0].index('FuelType'), newTechsCE[0].index('Capacity(MW)'))
    plantTypeCol = newTechsCE[0].index('TechnologyType')
    techSymbols = [createTechSymbol(row, newTechsCE[0], ptCurtailedAll) for row in newTechsCE]

    for gcm in hrlyCurtailmentsAllTechsInTgtYr:
        auxDict = dict()
        for (techSymbol, cell) in hrlyCurtailmentsAllTechsInTgtYr[gcm]:
            if getTechFromTechSymbol(techSymbol) in ptCurtailedAll:
                newTechsRow = newTechsCE[techSymbols.index(techSymbol)]
                hrlyCurtailments = hrlyCurtailmentsAllTechsInTgtYr[gcm][(techSymbol, cell)]
                hourlyCapacs = subtractCurtailmentsFromCapac(hrlyCurtailments,
                                                             float(newTechsRow[capacCol]), techSymbol + cell)
                auxDict[(techSymbol, cell)] = hourlyCapacs

        thermalTechHourlyCapacs[gcm] = auxDict

    return thermalTechHourlyCapacs
