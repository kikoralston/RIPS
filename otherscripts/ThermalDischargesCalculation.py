
from datetime import date, timedelta
import os, copy
from SetupGeneratorFleet import *

################# SET PARAMETERS ###############################################
def setParameters():
    #Parameters for this analysis
    years = [yr for yr in range(1995,2016)]
    # states = ['al']
    states = ['al','nc','sc','ky','tn','ms','ga']
    baseDir = 'C:\\Users\\mtcraig\\Desktop\\CEMSData'
    dirSave = os.path.join(baseDir,'OnceThroughThermalDischarges')
    dirReadBase = os.path.join(baseDir,'CSVFiles')
    #Parameters for fleet
    #GENERAL PARAMETERS
    testModel = False 
    fleetStates = ['North Carolina','South Carolina','Georgia','Mississippi','Alabama',
                'Kentucky','Tennessee'] 
    powerSystems =  ['S_SOU','S_VACA','S_C_KY','S_C_TVA','S_D_WOTA',
                     'S_D_N_AR','S_D_AMSO','S_D_REST'] 
    (startYear,endYear) = (2015,2015)
    compressFleet = True
    fuelPricesTimeSeries = importFuelPrices()
    return (testModel,fleetStates,powerSystems,endYear,startYear,fuelPricesTimeSeries,
        compressFleet,years,states,dirSave,dirReadBase)

def importFuelPrices():
    fuelPriceDir = 'C:\\Users\\mtcraig\\Desktop\\EPP Research\\Databases\\FuelPricesCapacityExpansion'
    fuelFileName = 'FuelPriceTimeSeries2Aug2016.csv'
    return readCSVto2dList(os.path.join(fuelPriceDir,fuelFileName))
################################################################################

################# MASTER FUNCTION ##############################################
def masterFunction():
    (testModel,fleetStates,powerSystems,endYear,startYear,fuelPricesTimeSeries,
        compressFleet,years,states,dirSave,dirReadBase) = setParameters()
    (fleetOnceThrough,onceThroughGenIds,onceThroughOrisIds,genIdToHr) = getOnceThroughPlants(testModel,fleetStates,
                                     powerSystems,endYear,startYear,fuelPricesTimeSeries,compressFleet)
    for year in years:
        (dischargeDataList,orisIdToRow,dateToCol) = setupDischargeDataList(year,fleetOnceThrough,onceThroughOrisIds)
        dischargeDataList = saveDischargeData(states,year,dischargeDataList,fleetOnceThrough,onceThroughGenIds,
                            orisIdToRow,dateToCol,genIdToHr,dirReadBase)
        if not os.path.exists(dirSave): os.makedirs(dirSave)
        write2dListToCSV(dischargeDataList,os.path.join(dirSave,'thermalDischargesCoalPlants' + str(year) + '.csv'))
################################################################################

################# GET ONCE THROUGH PLANTS ######################################
def getOnceThroughPlants(testModel,states,powerSystems,endYear,startYear,
                                fuelPricesTimeSeries,compressFleet):
    genFleetFull = setupGeneratorFleet(testModel,states,powerSystems,endYear,startYear,
                                fuelPricesTimeSeries,compressFleet)
    (fleetOnceThrough,onceThroughGenIds,onceThroughOrisIds,genIdToHr) = isolateOnceThroughPlants(genFleetFull)
    (coalPlantsOnceThrough,onceThroughGenIds,onceThroughOrisIds,genIdToHr) = isolateOnceThroughPlants(genFleetFull,True)
    write2dListToCSV(fleetOnceThrough,'onceThroughAllPlants.csv')
    write2dListToCSV(coalPlantsOnceThrough,'onceThroughCoalPlants.csv')
    return (coalPlantsOnceThrough,onceThroughGenIds,onceThroughOrisIds,genIdToHr)

def isolateOnceThroughPlants(genFleetFull,onlyCoalPlants=False):
    coolingTypeCol = genFleetFull[0].index('Cooling Tech')
    orisCol = genFleetFull[0].index('ORIS Plant Code')
    genIdCol = genFleetFull[0].index('Unit ID')
    hrCol = genFleetFull[0].index('Heat Rate (Btu/kWh)')
    plantTypeCol = genFleetFull[0].index('PlantType')
    onceThroughGenIds = set()
    onceThroughOrisIds = set()
    genIdToHr = dict()
    if onlyCoalPlants == False:
        onceThroughRows = [row for row in genFleetFull if 'once through' in row[coolingTypeCol]]
    else:
        onceThroughRows = [row for row in genFleetFull if ('once through' in row[coolingTypeCol]
                                                        and row[plantTypeCol] == 'Coal Steam')]
    for row in onceThroughRows:
        genId = createGenId(row[orisCol],row[genIdCol])
        onceThroughGenIds.add(genId)
        onceThroughOrisIds.add(str(row[orisCol]))
        genIdToHr[genId] = float(row[hrCol])
    fleetOnceThrough = [copy.copy(genFleetFull[0])] + onceThroughRows
    return (fleetOnceThrough,onceThroughGenIds,onceThroughOrisIds,genIdToHr)

def createGenId(orisId,unitId):
    return str(orisId) + '+' + str(unitId)
################################################################################

################# SET UP EMPTY LIST FOR RESULTS ################################
def setupDischargeDataList(year,fleetOnceThrough,onceThroughOrisIds):
    genHeads = ['OrisId','Lat','Long']
    (datesList,dateToCol) = datesInYear(year,genHeads)
    heads = genHeads + datesList
    dischargeDataList = [heads]
    orisCol = fleetOnceThrough[0].index('ORIS Plant Code')
    latCol = fleetOnceThrough[0].index('Latitude')
    longCol = fleetOnceThrough[0].index('Longitude')
    orisIds = [str(row[orisCol]) for row in fleetOnceThrough]
    (ctr,orisIdToRow) = (1,dict()) #start @ 1 for headers
    for orisId in onceThroughOrisIds:
        fleetRowNum = orisIds.index(orisId)
        (lat,lon) = (fleetOnceThrough[fleetRowNum][latCol],fleetOnceThrough[fleetRowNum][longCol])
        genData = [orisId,lat,lon]
        orisIdToRow[orisId] = ctr
        ctr += 1
        dischargeDataList.append(genData + ['NoData'] * len(datesList))
    return (dischargeDataList,orisIdToRow,dateToCol)

def datesInYear(year,genHeads):  
    startDate = date(year,1,1)
    endDate = date(year,12,31)
    deltaDays = endDate - startDate
    (dates,dateToCol) = ([],dict())
    ctr = 0
    for dayInc in range(deltaDays.days + 1):
        currDay = startDate + timedelta(days=dayInc)
        dayStr = createDayStr(currDay)
        dateToCol[dayStr] = ctr + len(genHeads)
        dates.append(dayStr)
        ctr += 1
    return (dates,dateToCol)

#Takes in date datetime object
def createDayStr(currDay):
    currMonth = currDay.month
    if currMonth < 10: monthStr = '0' + str(currMonth)
    else: monthStr = str(currMonth)
    tempDay = currDay.day
    if tempDay < 10: dayStr = '0' + str(tempDay)
    else: dayStr = str(tempDay)
    return monthStr + '-' + dayStr + '-' + str(currDay.year)
################################################################################

################# GET ONCE THROUGH PLANTS ######################################
def saveDischargeData(states,year,dischargeDataList,fleetOnceThrough,onceThroughGenIds,orisIdToRow,dateToCol,
                        genIdToHr,dirReadBase):
    # print(onceThroughGenIds)
    print('Processing year ' + str(year))
    months = getMonthStrs()
    for state in states:
        matchedUnitIds = set()
        for month in months:
            newFile = True
            currFile = str(year) + state + month + '.csv'
            with open(os.path.join(dirReadBase,str(year),currFile),'rt') as f:
                fReader = csv.reader(f, delimiter=',')
                for row in fReader:
                    if newFile == True:
                        orisCol = row.index('ORISPL_CODE')
                        genIdCol = row.index('UNITID')
                        for idx in range(len(row)): #not all heat input headesr are the same
                            if 'HEAT_INPUT' in row[idx]: heatInputCol = idx
                        dateCol = row.index('OP_DATE')
                        newFile = False
                    else:
                        (orisId,unitId) = (row[orisCol],createGenId(row[orisCol],row[genIdCol]))
                        if unitId in onceThroughGenIds:

                            matchedUnitIds.add(unitId)
                            currDate = row[dateCol]
                            (discRow,discCol) = (orisIdToRow[orisId],dateToCol[currDate])
                            # print(discRow,discCol)
                            if row[heatInputCol] != '': 
                                currDischarge = estimateThermalDischarge(float(row[heatInputCol]),genIdToHr[unitId])
                            else: 
                                currDischarge = 0
                            # print(currDischarge)
                            if dischargeDataList[discRow][discCol] == 'NoData':
                                dischargeDataList[discRow][discCol] = currDischarge
                            else:
                                dischargeDataList[discRow][discCol] += currDischarge
        print('Processed ' + state)
        # print('Matched units:',matchedUnitIds)
    return dischargeDataList 

def estimateThermalDischarge(heatInput,heatRate):
    effNumerator = 3412 #Btu/hr
    plantEff = effNumerator/heatRate
    wasteHeatAsThermalDischarge = 0.7
    return wasteHeatAsThermalDischarge * (1 - plantEff) * heatInput

#Create list of month suffixes on CEMS data files w/ leading zeros
def getMonthStrs():
    monthStrs = []
    for num in range(1,13):
        if num<10: monthStrs.append('0' + str(num))
        else: monthStrs.append(str(num))
    return monthStrs    
################################################################################

masterFunction()