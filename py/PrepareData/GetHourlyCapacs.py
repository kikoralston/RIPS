#Michael Craig, 17 May 2017
#Calculates hourly capacity curtailments of non-RE generators (existing
#and new techs) for CE model,
#then returns a dictionary of gen:hourly capac that accounts for curtailments
#for only those hours included in CE model.

from GAMSUtil.GAMSAuxFuncs import *
from thermalderatings.CalculateHourlyCapacsWithCurtailments import (calculateHourlyTechCapacsWithCurtailments,
                                                                    calculateHourlyCapacsWithCurtailments)


def getHourlyNonRECapacsForCE(genFleetForCE, hrlyCurtailmentsAllGensInTgtYr, hoursForCE, currYear):
    """Get hourly capacities for existing gens for CE run

    This function prepares a list with hourly capacities (in MW) of existing non renewable generators

    :param genFleetForCE: (2d list) generator fleet
    :param hrlyCurtailmentsAllGensInTgtYr: (dict) dictionary with hourly capacities (in % of nominal capacity) of existing generators
    :param hoursForCE: (list) list of hours in the year that will be included in CE run
    :param currYear: (int) current year of simulation
    :return: (dictionary) nested dictionary {gcm: {gen id: [hourly capacs]}}. lists with hourly capacities includes only hours of CE simulation
    """
    #Get hourly capacities
    hourlyCapacsAllGens = calculateHourlyCapacsWithCurtailments(genFleetForCE,
                                                                hrlyCurtailmentsAllGensInTgtYr, currYear)
    #Narrow down capacities to hours for CE
    hourlyCapacsCE = dict()

    for gcm in hourlyCapacsAllGens:
        auxDict = dict()
        for row in genFleetForCE[1:]:
            genSymbol = createGenSymbol(row, genFleetForCE[0])
            hourlyCapacs = [hourlyCapacsAllGens[gcm][genSymbol][hr-1] for hr in hoursForCE[gcm]]
            auxDict[genSymbol] = hourlyCapacs
        hourlyCapacsCE[gcm] = auxDict

    return hourlyCapacsCE


def getHourlyCurtailedTechCapacsForCE(newTechsCE,hrlyCurtailmentsAllTechsInTgtYr, hoursForCE,currYear,
                                      plantTypesCurtailed):
    """ Trim hours of techs capacity factors to hours of CE run

    :param newTechsCE: (2d list) list with parameters of candidate techs
    :param hrlyCurtailmentsAllTechsInTgtYr: (dict) dictionary with hourly capacities (in % of nominal capacity) of candidate techs in each location
    :param hoursForCE: (list) list of hours in the year that will be included in CE run
    :param currYear: (int) current year of simulation
    :param plantTypesCurtailed: (list) names of types of plants that suffer climate-related deratings
    :return: dict of {(tech+cooltype,cell): [hourly capacs]}. lists with hourly capacities includes only hours of CE simulation
    """
    #Get hourly capacities
    thermalTechHourlyCapacs = calculateHourlyTechCapacsWithCurtailments(newTechsCE, hrlyCurtailmentsAllTechsInTgtYr,
                                                                        currYear, plantTypesCurtailed)
    #Narrow down capacities to hours for CE
    thermalTechHourlyCapacsCE = dict()

    for gcm in thermalTechHourlyCapacs:
        auxDict = dict()
        for key in thermalTechHourlyCapacs[gcm]:
            hourlyCapacs = [thermalTechHourlyCapacs[gcm][key][hr-1] for hr in hoursForCE[gcm]]
            auxDict[key] = hourlyCapacs

        thermalTechHourlyCapacsCE[gcm] = auxDict

    return thermalTechHourlyCapacsCE


def getHourlyCapacitiesForDays(fleetUC, hourlyCapacsAllGens, hoursForUC):
    """Get hourly capacity values for UCED simulation day

    This function isolates hourly capacities for generators just to the hours of the year included in the current
    daily UCED simulation.

    :param fleetUC: (2-d list) generator fleet
    :param hourlyCapacsAllGens:
    :param hoursForUC: (1-d list) list with hours in year (1-8760) of current UCED run
    :return: (dict) dictionary with hourly capacity values (in MW) for each generator, i.e. {gen: [capacs]}
    """
    (hourlyCapacsUC, hourlyCapacsUCList) = (dict(), [])

    # In order to treat cases where the hours of the year extend beyond 8760
    # repeat hours of last day of the year
    hoursForUCAux = [((8760 - 24 + (hr % 24 if hr % 24 > 0 else 24)) if hr > 8760 else hr) for hr in hoursForUC]

    for row in fleetUC[1:]:
        genSymbol = createGenSymbol(row, fleetUC[0])

        hourlyCapacsDay = [hourlyCapacsAllGens[genSymbol][hr-1] for hr in hoursForUCAux]

        hourlyCapacsUC[genSymbol] = hourlyCapacsDay
        hourlyCapacsUCList.append([genSymbol] + hourlyCapacsDay)

    # write2dListToCSV(hourlyCapacsUCList,'hourlyCapacsUCDay' + str(hoursForUC[0]//24+1) + '.csv')
    return hourlyCapacsUC


