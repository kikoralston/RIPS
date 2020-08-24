#Michael Craig
#October 12, 2016
#Remove hydro (normal + pumped storage) units from fleet and subtract their monthly average generation
#from demand profile.

import copy, operator, os
from AuxFuncs import *
from SetupGeneratorFleet import *

###### UNIVERSAL PARAMETERS ####################################################
def setKeyParameters():
    #GENERAL PARAMETERS
    testModel = False #use dummy test system; not currently working
    onlyTennessee = False
    states = ['North Carolina','South Carolina','Georgia','Mississippi','Alabama',
                'Kentucky','Tennessee'] 
    statesAbbrev = ['NC','SC','GA','MS','AL','KY','TN']
    if onlyTennessee==True: (states,statesAbbrev) = (['Tennessee'],['TN'])
    powerSystems =  ['S_SOU','S_VACA','S_C_KY','S_C_TVA','S_D_WOTA',
                     'S_D_N_AR','S_D_AMSO','S_D_REST'] 
    compressFleet = True
    fuelPricesTimeSeries = importFuelPrices()
    (co2CppSerc2022Limit,co2CppSerc2030Limit) = calcRegionCPPLimit(states) #short tons/yr; includes new source complement
    #THERMAL CURTAILMENT PARAMETERS
    processRBMData = False
    #RENEWABLE CAPACITY FACTOR PARAMETERS
    tzAnalysis = 'EST'
    projectName = 'rips'
    #CAPACITY EXPANSION PARAMETERS
    runCapacExp = True
    capacExpFilename = 'CEChrono15Aug2016.gms' #CEChrono23June2016,CEChrono23June2016EffCapac,CEChrono23June2016EffCapacWithWeights
    (startYear,endYear,yearStepCE) = (2015,2025,5)
    retirementCFCutoff = .2 #retire units w/ CF lower than given value
    daysPerSeason = 3
    selectCurtailDays = True
    planningReserveMargin = 0.15 #fraction of peak demand
    discountRate = 0.07 #fraction    
    allowCoalWithoutCCS = False
    onlyNSPSUnits = True
    cellsEligibleForNewPlants = 'WithGens' #'All' (all cells) or 'WithGen' (only cells already with gen inside)
    #UNIT COMMITMENT PARAMETERS
    calculateCO2Price = True
    #CONVERSION PARAMETERS
    scaleMWtoGW = 1000
    scaleDollarsToThousands = 1000
    scaleLbToShortTon = 2000
    #TEMPORARY PARAMETERS
    annualDemandGrowth = 0.001 #fraction per year
    return (testModel,annualDemandGrowth,startYear,yearStepCE,endYear,daysPerSeason,selectCurtailDays,
            states,statesAbbrev,powerSystems,planningReserveMargin,discountRate,fuelPricesTimeSeries,
            scaleMWtoGW,scaleDollarsToThousands,co2CppSerc2022Limit,co2CppSerc2030Limit,
            compressFleet,runCapacExp,allowCoalWithoutCCS,capacExpFilename,onlyNSPSUnits,
            cellsEligibleForNewPlants,calculateCO2Price,retirementCFCutoff,processRBMData,scaleLbToShortTon,
            tzAnalysis,projectName)

########### GET FUEL PRICES ####################################################
def importFuelPrices():
    fuelPriceDir = 'C:\\Users\\mtcraig\\Desktop\\EPP Research\\Databases\\FuelPricesCapacityExpansion'
    fuelFileName = 'FuelPriceTimeSeries2Aug2016.csv'
    return readCSVto2dList(os.path.join(fuelPriceDir,fuelFileName))

########### GET CPP MASS LIMITS ################################################
def calcRegionCPPLimit(states):
    stateCppLimits = importStateCppLimits()
    stateCol = stateCppLimits[0].index('State')
    (limit2022Col,limit2030Col) = (stateCppLimits[0].index('2022'),stateCppLimits[0].index('2030'))
    state2022Limits = [float(row[limit2022Col]) for row in stateCppLimits[1:] if row[stateCol] in states]
    state2030Limits = [float(row[limit2030Col]) for row in stateCppLimits[1:] if row[stateCol] in states]
    return (sum(state2022Limits),sum(state2030Limits))

#State CPP limits in short tons
def importStateCppLimits():
    cppLimitDir = 'C:\\Users\\mtcraig\\Desktop\\EPP Research\\Power Plant Rules\\Clean Power Plan\\FinalRule'
    cppLimitName = 'StateCPPCO2EmsCaps5Oct2016.csv'
    return readCSVto2dList(os.path.join(cppLimitDir,cppLimitName))
################################################################################

########### MASTER FUNCTION ####################################################
def setKeyParametersUnique():
    eia923dir = 'C:\\Users\\mtcraig\\Desktop\\EPP Research\\Databases\\EIAForm923\\2011to2015Compiled'
    eia923years = list(range(2011,2015))
    return (eia923years,eia923dir)

#Inputs: gen fleet, demand profile (1d list of hourly demand values)
#Outputs: gen fleet w/out hydro units, demand profile w/ average monthly 
#hydro generation subtracted (1d list of hourly values)
def getHistoricHydroGeneration():
    #Get fleet
    (testModel,annualDemandGrowth,startYear,yearStepCE,endYear,daysPerSeason,selectCurtailDays,
            states,statesAbbrev,powerSystems,planningReserveMargin,discountRate,fuelPricesTimeSeries,
            scaleMWtoGW,scaleDollarsToThousands,co2CppSerc2022Limit,co2CppSerc2030Limit,compressFleet,
            runCapacExp,allowCoalWithoutCCS,capacExpFilename,onlyNSPSUnits,cellsEligibleForNewPlants,
            calculateCO2Price,retirementCFCutoff,processRBMData,scaleLbToShortTon,
            tzAnalysis,projectName) = setKeyParameters()
    genFleet = getInitialFleetAndDemand(testModel,states,powerSystems,endYear,
                                                  startYear,fuelPricesTimeSeries,compressFleet,
                                                  statesAbbrev)
    #Get average hydro monthly gen
    (eia923years,eia923dir) = setKeyParametersUnique()
    (hydroFleetRows,genFleetNoHydro) = getHydroRows(genFleet) #hydroFleetRows has fleet header & full fleet rows w/ hydro plants
    (orisIDtoCapac,orisIdToLatLon) = getHydroOrisIdsAndCapacs(hydroFleetRows)
    hydroAvgMonthlyGen = getTotalHydroAvgMonthlyGen(orisIDtoCapac,eia923years,eia923dir)
    hydroAvgMonthlyGenList = convertDictTo2dList(hydroAvgMonthlyGen,orisIDtoCapac,orisIdToLatLon)
    write2dListToCSV(hydroAvgMonthlyGenList,'hydroAvgAndPotMonthlyGen.csv')
################################################################################

########### GET GEN FLEET #####################################################
def getInitialFleetAndDemand(testModel,states,powerSystems,endYear,startYear,
                            fuelPricesTimeSeries,compressFleet,statesAbbrev):
    genFleet = setupGeneratorFleet(testModel,states,powerSystems,endYear,startYear,
                                fuelPricesTimeSeries,compressFleet)
    print('Built initial generator fleet')
    # demandProfile = importHourlyDemand(testModel,states,statesAbbrev) #1d list, MWh
    # write2dListToCSV([demandProfile],'hourlyDemandInitial.csv')
    # print('Imported initial demand profile')    
    return genFleet
################################################################################

########### GET HYDRO ROWS #####################################################
#Returns hydro units in separate fleet, and fleet without hydro units
def getHydroRows(fleet):
    genFleetNoHydro = copy.deepcopy(fleet)
    heads = copy.copy(genFleetNoHydro[0])
    fuelTypeCol = heads.index('Modeled Fuels')
    hydroFleetRows = [heads]
    idxs = []
    for idx in range(len(genFleetNoHydro)-1,0,-1):
        if genFleetNoHydro[idx][fuelTypeCol] == 'Hydro': 
            hydroRow = genFleetNoHydro.pop(idx)
            hydroFleetRows.append(hydroRow)
    return (hydroFleetRows,genFleetNoHydro)
################################################################################

########### GET HYDRO ORIS ID AND CAPACITIES FROM FLEET ########################
#Return dictionary mapping ORIS ID to capacity
def getHydroOrisIdsAndCapacs(hydroFleetRows):
    (orisCol,capacCol) = (hydroFleetRows[0].index('ORIS Plant Code'),
                        hydroFleetRows[0].index('Capacity (MW)'))
    (latCol,lonCol) = (hydroFleetRows[0].index('Latitude'),
                        hydroFleetRows[0].index('Longitude'))
    orisIDtoCapac = dict()
    orisIdToLatLon = dict()
    for row in hydroFleetRows[1:]:
        if row[orisCol] in orisIDtoCapac: orisIDtoCapac[row[orisCol]] += float(row[capacCol])
        else: orisIDtoCapac[row[orisCol]] = float(row[capacCol])
        (lat,lon) = (row[latCol],row[lonCol])
        orisIdToLatLon[row[orisCol]] = (lat,lon)
    return (orisIDtoCapac,orisIdToLatLon)
################################################################################

########### GET HISTORIC EIA923 GENERATION FOR EACH HYDRO UNIT #################
#Get total average monthly hydro generation by getting monthly generation per year
#for each unit, then adding average values.
#Inputs: dict mapping ORIS id to capac of all hydro units, years of EIA 923
#data to use, & dir of that data.
#Outputs: 1d list (len=12) of average total hydro generation per month
def getTotalHydroAvgMonthlyGen(orisIDtoCapac,eia923years,eia923dir):
    (orisIDtoMonthlyGen,orisIDtoMonthlyGenCount) = (dict(),dict())
    for orisId in orisIDtoCapac: 
        orisIDtoMonthlyGen[orisId] = []
        orisIDtoMonthlyGenCount[orisId] = []
    for year in eia923years:
        (orisIDtoMonthlyGen,orisIDtoMonthlyGenCount) = getMonthlyGenInYear(orisIDtoMonthlyGen,
                                                        orisIDtoMonthlyGenCount,year,eia923dir)
    return getAverageMonthlyGen(orisIDtoMonthlyGen,orisIDtoMonthlyGenCount)
     
#For each hydro unit in fleet, get monthly generation in a given year.
#Inputs: dict mapping oris ID to list of monthly gen, dict mapping ORIS ID to 
#list of count of gen values per month, year of analysis, dir w/ EIA 923 data
#Outputs: dict mapping oris ID to list of monthly gen, dict mapping ORIS ID to 
#list of count of gen values per month
def getMonthlyGenInYear(orisIDtoMonthlyGen,orisIDtoMonthlyGenCount,year,eia923dir):
    numMonths = 12
    (idCol,idLabel) = (0,'Plant Id')
    genFile = 'gen' + str(year) + '.csv'
    genData = readCSVto2dList(os.path.join(eia923dir,genFile))
    firstColVals = [row[idCol] for row in genData]
    headsRow = firstColVals.index(idLabel) #this has detailed header rows; 1 row up are overarching headers, hence -1 in next line
    netGenFirstCol = genData[headsRow-1].index('Electricity Net Generation (MWh)')
    if 'Reported Fuel Type Code' in genData[headsRow]: fuelCol = genData[headsRow].index('Reported Fuel Type Code') 
    else: fuelCol = genData[headsRow].index('Reported\nFuel Type Code') 
    for row in genData[headsRow+1:]:
        (orisId,fuel) = (row[idCol],row[fuelCol])
        if orisId in orisIDtoMonthlyGen and fuel == 'WAT':
                monthlyGen = [toNum(row[idx]) for idx in range(netGenFirstCol,netGenFirstCol+numMonths)]
                if orisIDtoMonthlyGen[orisId]==[]:
                    orisIDtoMonthlyGen[orisId] = monthlyGen
                    orisIDtoMonthlyGenCount[orisId] = [1]*len(monthlyGen)
                else:
                    orisIDtoMonthlyGen[orisId] = list(map(operator.add,orisIDtoMonthlyGen[orisId],monthlyGen))
                    orisIDtoMonthlyGenCount[orisId] = list(map(operator.add,orisIDtoMonthlyGenCount[orisId],[1]*len(monthlyGen)))
    return (orisIDtoMonthlyGen,orisIDtoMonthlyGenCount)

#For each unit, calculate average monthly gen, then return dict 
#mapping ORIS ID to average monthly gen.
#Inputs: dict mapping oris ID to list of monthly gen, dict mapping ORIS ID to 
#list of count of gen values per month
def getAverageMonthlyGen(orisIDtoMonthlyGen,orisIDtoMonthlyGenCount):
    combinedAverageGen = []
    orisIDtoAvgMonthlyGen = dict()
    for orisId in orisIDtoMonthlyGen:
        (monthlyGen,count) = (orisIDtoMonthlyGen[orisId],orisIDtoMonthlyGenCount[orisId])
        averageMonthlyGen = [monthlyGen[idx]/count[idx] for idx in range(len(monthlyGen))]
        orisIDtoAvgMonthlyGen[orisId] = averageMonthlyGen
    return orisIDtoAvgMonthlyGen
################################################################################

#Convert dict to 2d list
def convertDictTo2dList(hydroAvgMonthlyGen,orisIDtoCapac,orisIdToLatLon):
    months = ['Jan','Feb','Mar','Apr','May','Jun','Jul',
                            'Aug','Sep','Oct','Nov','Dec']
    daysPerMonth = [31,28,31,30,31,30,31,31,30,31,30,31]
    avgHeads = [month + 'AvgGen(MWh)' for month in months]
    potentialHeads = [month + 'PotentialGen(MWh)' for month in months]
    heads = ['OrisId','Latitude','Longitude'] + avgHeads + potentialHeads 
    hydroAvgAndPotentialGenList = [copy.copy(heads)]
    for orisId in hydroAvgMonthlyGen:
        averageMonthlyGen = hydroAvgMonthlyGen[orisId]
        (lat,lon) = orisIdToLatLon[orisId]
        capacity = orisIDtoCapac[orisId]
        potentialGen = [float(capacity) * 24 * val for val in daysPerMonth]
        #Some average monthly gen values > potential gen - cap them
        for idx in range(len(averageMonthlyGen)):
            if averageMonthlyGen[idx] > potentialGen[idx]: averageMonthlyGen[idx] = potentialGen[idx]
        hydroAvgAndPotentialGenList.append([orisId] + [lat] + [lon] + averageMonthlyGen + potentialGen)
    return hydroAvgAndPotentialGenList

########### HELPER FUNCTION ####################################################
#Converts a string w/ commas in it to a float
def toNum(s):
    if s=='.': return 0
    numSegments = s.split(',')
    result = ""
    for segment in numSegments:
        result += segment
    return float(result)
################################################################################

getHistoricHydroGeneration()