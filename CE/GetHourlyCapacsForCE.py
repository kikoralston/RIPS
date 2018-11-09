#Michael Craig, 17 May 2017
#Calculates hourly capacity curtailments of non-RE generators (existing
#and new techs) for CE model,
#then returns a dictionary of gen:hourly capac that accounts for curtailments
#for only those hours included in CE model.

from AuxFuncs import *
from GAMSAuxFuncs import *
from CalculateHourlyCapacsWithCurtailments import (calculateHourlyTechCapacsWithCurtailments,
                calculateHourlyCapacsWithCurtailments)
import os


def getHourlyNonRECapacsForCE(genFleetForCE, hrlyCurtailmentsAllGensInTgtYr, hoursForCE,currYear):
    """ GET HOURLY CAPACITIES FOR EXISTING GENS FOR CE

    :param genFleetForCE:
    :param hrlyCurtailmentsAllGensInTgtYr:
    :param hoursForCE:
    :param currYear:
    :return:
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


def getHourlyCurtailedTechCapacsForCE(newTechsCE,hrlyCurtailmentsAllTechsInTgtYr,
                            hoursForCE,currYear,plantTypesCurtailed):
    """ TRIM HOURS OF TECH CFS TO HOURS OF CE

    :param newTechsCE:
    :param hrlyCurtailmentsAllTechsInTgtYr:
    :param hoursForCE:
    :param currYear:
    :param plantTypesCurtailed:
    :return: dict of {(tech+cooltype,cell): [hourly capacs]}
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

