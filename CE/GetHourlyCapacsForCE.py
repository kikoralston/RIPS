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

########### GET HOURLY CAPACITIES FOR EXISTING GENS FOR CE #####################
def getHourlyNonRECapacsForCE(genFleetForCE,hrlyCurtailmentsAllGensInTgtYr,
                                hoursForCE,currYear):
    #Get hourly capacities
    (hourlyCapacsAllGens,hourlyCapacsAllGensList) = calculateHourlyCapacsWithCurtailments(genFleetForCE,
                                                            hrlyCurtailmentsAllGensInTgtYr,currYear)
    #Narrow down capacities to hours for CE
    (hourlyCapacsCE,hourlyCapacsCEList) = (dict(),[])
    for row in genFleetForCE[1:]:
        genSymbol = createGenSymbol(row,genFleetForCE[0])
        hourlyCapacs = [hourlyCapacsAllGens[genSymbol][hr-1] for hr in hoursForCE]
        hourlyCapacsCE[genSymbol] = hourlyCapacs
        hourlyCapacsCEList.append([genSymbol] + hourlyCapacs)
    return hourlyCapacsCE,hourlyCapacsCEList

########### TRIM HOURS OF TECH CFS TO HOURS OF CE #####################
#Returns dict of (tech+cooltype,cell):[hourly capacs]
def getHourlyCurtailedTechCapacsForCE(newTechsCE,hrlyCurtailmentsAllTechsInTgtYr,
                            hoursForCE,currYear,plantTypesCurtailed):
    #Get hourly capacities
    thermalTechHourlyCapacs,thermalTechHourlyCapacsList = calculateHourlyTechCapacsWithCurtailments(newTechsCE,
                    hrlyCurtailmentsAllTechsInTgtYr,currYear,plantTypesCurtailed)
    #Narrow down capacities to hours for CE
    (thermalTechHourlyCapacsCE,hourlyThermalTechCapacsCEList) = (dict(),[])
    for key in thermalTechHourlyCapacs:
        hourlyCapacs = [thermalTechHourlyCapacs[key][hr-1] for hr in hoursForCE]
        thermalTechHourlyCapacsCE[key] = hourlyCapacs
        hourlyThermalTechCapacsCEList.append([key] + hourlyCapacs)
    return thermalTechHourlyCapacsCE,hourlyThermalTechCapacsCEList,thermalTechHourlyCapacsList

