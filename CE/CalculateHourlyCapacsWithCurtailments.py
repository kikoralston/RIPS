#Michael Craig, 17 May 2017
#Calculate hourly capacity including curtailments for existing generators.
#Outputs dict of gen:hourly capacities and a 2d list w/ generator ID followed by curtailments.

from AuxFuncs import *
from GAMSAuxFuncs import *
import os

########### GET EXISTING GEN HOURLY CAPACITIES W/ THERMAL CONSTRAINTS ##########
#Return dictionary of gen symbol to 1d list of hourly capacity for year
def calculateHourlyCapacsWithCurtailments(genFleet,hrlyCurtailmentsAllGensInTgtYr,currYear):
    (genHourlyCapacs,genHourlyCapacsList) = (dict(),[])
    (fuelTypeCol,capacCol) = (genFleet[0].index('Modeled Fuels'),genFleet[0].index('Capacity (MW)'))
    lenCurtailments = len(hrlyCurtailmentsAllGensInTgtYr[list(hrlyCurtailmentsAllGensInTgtYr.keys())[0]])
    for row in genFleet[1:]:
        genSymbol = createGenSymbol(row,genFleet[0])
        if row[fuelTypeCol] == 'Wind' or row[fuelTypeCol] == 'Solar': 
            hourlyCapacs = [float(row[capacCol])] * lenCurtailments
        elif genSymbol in hrlyCurtailmentsAllGensInTgtYr:
            hrlyCurtailments = hrlyCurtailmentsAllGensInTgtYr[genSymbol] 
            hourlyCapacs = subtractCurtailmentsFromCapac(hrlyCurtailments,float(row[capacCol]),genSymbol)
        else: #no curtailments
            hourlyCapacs = [float(row[capacCol])] * lenCurtailments #has header 
        genHourlyCapacs[genSymbol] = hourlyCapacs
        genHourlyCapacsList.append([genSymbol] + hourlyCapacs)
    return (genHourlyCapacs,genHourlyCapacsList)

#Hourly curtailments are fraction fo total capacity. 
def subtractCurtailmentsFromCapac(hrlyCurtailments,capac,genSymbol):
    hourlyCapacs = [float(val)*capac for val in hrlyCurtailments]
    if min(hourlyCapacs)<0: 
        print('Hourly capacity < 0 for ',genSymbol,' at idx ',hourlyCapacs.index(min(hourlyCapacs)))
        hourlyCapacs = [max(0,val) for val in hrlyCurtailments]
    return hourlyCapacs

########### GET NEW TECH HOURLY CAPACITIES W/ THERMAL CONSTRAINTS ##############
def calculateHourlyTechCapacsWithCurtailments(newTechsCE,hrlyCurtailmentsAllTechsInTgtYr,
            currYear,ptCurtailedAll):
    (thermalTechHourlyCapacs,thermalTechHourlyCapacsList) = (dict(),[])
    (fuelTypeCol,capacCol) = (newTechsCE[0].index('FuelType'),newTechsCE[0].index('Capacity(MW)'))
    plantTypeCol = newTechsCE[0].index('TechnologyType')
    techSymbols = [createTechSymbol(row,newTechsCE[0],ptCurtailedAll) for row in newTechsCE]
    for (techSymbol,cell) in hrlyCurtailmentsAllTechsInTgtYr:
        if getTechFromTechSymbol(techSymbol) in ptCurtailedAll:
            newTechsRow = newTechsCE[techSymbols.index(techSymbol)]
            hrlyCurtailments = hrlyCurtailmentsAllTechsInTgtYr[(techSymbol,cell)] 
            hourlyCapacs = subtractCurtailmentsFromCapac(hrlyCurtailments,
                                        float(newTechsRow[capacCol]),techSymbol+cell)
            thermalTechHourlyCapacs[(techSymbol,cell)] = hourlyCapacs
            thermalTechHourlyCapacsList.append([(techSymbol,cell)] + hourlyCapacs)
    return thermalTechHourlyCapacs,thermalTechHourlyCapacsList


