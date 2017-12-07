#Michael Craig
#July 21, 2017
#Import monthly hydropower generation from PNNL, and assign to each season
#considered in CE model.

import copy, os
from AuxFuncs import *
from GAMSAuxFuncs import createGenSymbol
from DemandFuncsCE import getHoursInMonths

#Inputs: gen fleet, zonal demand (dict), rep hrs per season or special block (dict of block:rep hrs,
#where seasons are labled by name and special hours are 'special'), curr year, and run location.
#Outputs: dict of season:genSymbol:generation potential
def getHydroEPotential(fleet,demandZonal,repAndSpeHoursDict,currYear,runLoc):
    hydroPotentials = importHydroPotentialGen(currYear,runLoc)
    hydroPotPerSeason = dict()
    plantCol,orisCol = fleet[0].index('PlantType'),fleet[0].index('ORIS Plant Code')
    zoneCol = fleet[0].index('Region Name')
    hydroUnits = [row for row in fleet[1:] if row[plantCol] == 'Hydro']
    for season in repAndSpeHoursDict:
        repHrs = repAndSpeHoursDict[season]
        months = getMonthsOfRepHrs(repHrs)
        monthlyDemandZonal = getMonthlyDemand(demandZonal,months)
        repHrsDemandZonal = getRepHrsDemand(demandZonal,repHrs)
        seasonDict,unitsNoData = dict(),list()
        for row in hydroUnits:
            oris,zone,genSymbol = row[orisCol],row[zoneCol],createGenSymbol(row,fleet[0])
            potential,hasData = getMonthsPotential(oris,hydroPotentials,months)
            if hasData == False: unitsNoData.append(genSymbol)
            else: seasonDict[genSymbol] = potential*repHrsDemandZonal[zone]/monthlyDemandZonal[zone]
        if len(unitsNoData)>0: assignPotentialsToMissingUnits(unitsNoData,fleet,seasonDict,repHrs)
        hydroPotPerSeason[season] = seasonDict
    return hydroPotPerSeason

#Returns hydro potential gen for current year in 2d list w/ col 1 = plant IDs
def importHydroPotentialGen(currYear,runLoc):
    if runLoc == 'pc': dataDir = 'C:\\Users\\mtcraig\\Desktop\\EPP Research\\Databases\\HydroMonthlyDataPNNL'
    else: dataDir = os.path.join('Data','HydroMonthlyDataPNNL')
    potentialsAllYears = readCSVto2dList(os.path.join(dataDir,'monhydrogen_1550.csv'))
    #First row in data is months listed w/out years; assign years to each col
    startYr,endYr = [2015,2050]
    yrIndex = (currYear-startYr)*12+1 #+1 to account for first col = PlantID
    return [[row[0]] + row[yrIndex:yrIndex+12] for row in potentialsAllYears]

#Inputs: 1d list of hours in a given season. 
#Function determines which months those hours fall in, and outputs those months in a 1d list.
def getMonthsOfRepHrs(repHrs):
    months = list()
    for month in range(1,13):
        monthHours = getHoursInMonths([month]) #returns 1d list of hours in given month
        if repHrs[0] in monthHours: months.append(month)
        elif repHrs[-1] in monthHours: months.append(month)
    return months

#Inputs: 1d list of hourly demand, 1d list of months
#Outputs: total demand in those months
def getMonthlyDemand(demandZonal,months):
    monthsHours = getHoursInMonths(months)
    demandInMonthsZonal = dict()
    for zone in demandZonal: demandInMonthsZonal[zone] = sumDemandInHours(demandZonal[zone],monthsHours)
    return demandInMonthsZonal

#Inputs: 1d list of demand zonal (dict of zone:hourly demand), 1d list of rep hrs per season
#Outputs: dict of zone:total demand in those rep hrs
def getRepHrsDemand(demandZonal,repHrs):
    demandInRepHrsZonal = dict()
    for zone in demandZonal: demandInRepHrsZonal[zone] = sumDemandInHours(demandZonal[zone],repHrs)
    return demandInRepHrsZonal

#Inputs: 1d lists of demand (0-8759 idx) and hours (1-8760 idx).
#Outputs: sum of demand in those hours
def sumDemandInHours(demand,hours):
    return sum([demand[hr-1] for hr in hours])

#Gets total hydropower generation potential for month(s).
#Inputs: ORIS ID, hydro potential generation (2d list) for curr year, months of interest
#Outputs: total potential generation for month(s)
def getMonthsPotential(oris,hydroPotentials,months):
    orisCol = hydroPotentials[0].index('PlantID')
    potentialOriss = [row[orisCol] for row in hydroPotentials]
    monthsPotentials,hasData = list(),True
    if oris in potentialOriss: 
        unitRowIdx = potentialOriss.index(oris)
        allMonthsPotentials = hydroPotentials[unitRowIdx] #leave ORIS label on so idx & month (1-12) align
        monthsPotentials = [float(allMonthsPotentials[month]) if float(allMonthsPotentials[month])>0 else 0 for month in months]
    else: hasData = False
    return sum(monthsPotentials),hasData

#For ORIS units not in PNNL data, assign them potential for season based on
#average potential as CF of rest of fleet for that season.
#Inputs: gen symbols (oris+unit ID) for units w/out data, gen fleet, dict of season:genSymbol:potential
#Outputs: modify seasonDict to include other units
def assignPotentialsToMissingUnits(unitsNoData,fleet,seasonDict,repHrs):
    plantCol = fleet[0].index('PlantType')
    capacCol = fleet[0].index('Capacity (MW)')
    hydroUnits = [row for row in fleet[1:] if row[plantCol] == 'Hydro']
    genSymbols = [createGenSymbol(row,fleet[0]) for row in hydroUnits]
    capacs = [row[capacCol] for row in hydroUnits]
    fleetAvgCF = getFleetAverageCF(genSymbols,capacs,seasonDict,repHrs)
    for genSymbol in unitsNoData: 
        seasonDict[genSymbol] = capacs[genSymbols.index(genSymbol)]*len(repHrs)*fleetAvgCF

#Gets average fleet CF for rep hours in given season.
#Inputs: 1d list of all hydro gen symbols and capac, dict of season:gensymbol:potential, rep hours (1d list)
#Outputs: average fleet CF in season
def getFleetAverageCF(genSymbols,capacs,seasonDict,repHrs):
    cfs = list()
    for (genSymbol,capac) in zip(genSymbols,capacs):
        if genSymbol in seasonDict:
            potential = seasonDict[genSymbol]
            cfs.append(potential/(capac*len(repHrs)))
    return sum(cfs)/len(cfs)


