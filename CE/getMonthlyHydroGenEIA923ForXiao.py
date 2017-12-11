#Michael Craig
#Feb 14, 2017
#Get monthly generation data for hydropower plants in SERC region from 2011-2015
#for Xiao.

import copy, operator, os
from AuxFuncs import *

def setKeyParameters():
    eia923dir = 'C:\\Users\\mtcraig\\Desktop\\EPP Research\\Databases\\EIAForm923\\2011to2015Compiled'
    eia923years = list(range(2011,2015))
    genFleetName = 'genFleetPreCompression.csv'
    return (eia923years,eia923dir,genFleetName)

########### MASTER FUNCTION ####################################################
#Inputs: gen fleet, demand profile (1d list of hourly demand values)
#Outputs: gen fleet w/out hydro units, demand profile w/ average monthly 
#hydro generation subtracted (1d list of hourly values)
def getHydroMonthlyGenValues():
    (eia923years,eia923dir,genFleetName) = setKeyParameters()
    genFleet = readCSVto2dList(genFleetName)
    (hydroFleetRows) = getHydroRows(genFleet) #hydroFleetRows has fleet header & full fleet rows w/ hydro plants
    saveMonthlyGen(hydroFleetRows,eia923years,eia923dir)
    # demandMinusHydroGen = removeHydroGenFromDemand(hydroFleetRows,demandProfile,eia923years,eia923dir)
    # return (genFleetNoHydro,demandMinusHydroGen)
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
            hydroFleetRows.append(genFleetNoHydro[idx])
    return (hydroFleetRows)
################################################################################

########### SAVE MONTHLY GEN FOR ALL EIA YEARS #################################
def saveMonthlyGen(hydroFleetRows,eia923years,eia923dir):
    hydroOrisIDs = getHydroOrisIdsAndCapacs(hydroFleetRows)
    print('**',hydroOrisIDs,len(hydroOrisIDs))
    hydroMonthlyGen = initializeMonthlyGenList(hydroOrisIDs,eia923years)
    print('**',hydroMonthlyGen)
    for year in eia923years:
        getMonthlyGenInYear(hydroMonthlyGen,year,eia923dir)
    write2dListToCSV(hydroMonthlyGen,'hydroMonthlyGenForXiao.csv')

#Return list of hydro ORIS IDs
def getHydroOrisIdsAndCapacs(hydroFleetRows):
    (orisCol,capacCol) = (hydroFleetRows[0].index('ORIS Plant Code'),hydroFleetRows[0].index('Capacity (MW)'))
    hydroOrisIDs = list()
    for row in hydroFleetRows[1:]: 
        if row[orisCol] not in hydroOrisIDs:
            hydroOrisIDs.append(row[orisCol])
    return hydroOrisIDs

def initializeMonthlyGenList(hydroOrisIDs,eia923years):
    hydroMonthlyGen = [['YearMonth/HydroORISIDs']]
    for oris in hydroOrisIDs: hydroMonthlyGen.append([oris])
    print(hydroMonthlyGen)
    months = ['jan','feb','mar','apr','may','jun','jul','aug','sep','oct','nov','dec']
    for year in eia923years:
        for month in months: 
            hydroMonthlyGen[0].extend([month + str(year)]) 
            # hydroMonthlyGen[0].extend([month + str(year)] + [0 for i in range(len(hydroOrisIDs))]) 
    return hydroMonthlyGen
     
def getMonthlyGenInYear(hydroMonthlyGen,year,eia923dir):
    numMonths = 12
    matchFound = [1] + [0 for i in range(1,len(hydroMonthlyGen))]
    hydroMonthlyGenRowLabels = [row[0] for row in hydroMonthlyGen]
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
        if orisId in hydroMonthlyGenRowLabels and fuel == 'WAT':
            monthlyGenRow = hydroMonthlyGenRowLabels.index(orisId)
            hydroMonthlyGen[monthlyGenRow].extend([toNum(row[idx]) for idx in 
                                        range(netGenFirstCol,netGenFirstCol+numMonths)])
            matchFound[monthlyGenRow] = 1
    for idx in range(len(matchFound)):
        if matchFound[idx] == 0: hydroMonthlyGen[idx].extend([0 for val in range(numMonths)])
    return (hydroMonthlyGen)
################################################################################

########### HELPER FUNCTION ####################################################
#Converts a string w/ commas in it to a float
def toNum(s):
    if s == '.': return 0
    else: 
        numSegments = s.split(',')
        result = ""
        for segment in numSegments:
            result += segment
        return float(result)

#Read CSV to 2d list
#Input: full file name including dir (str)
#Output: 2d list
def readCSVto2dList(fileNameWithDir):
    with open(fileNameWithDir,'r') as f:
        f = csv.reader(f)
        f = list(f)
    return f

#Write 2d list to CSV
#Input: 2d list, full file name including dir & .csv (str)
def write2dListToCSV(list2d,fileNameWithDir):
    fullFileName=fileNameWithDir
    with open(fullFileName, 'w', newline='') as csvfile:
        w = csv.writer(csvfile)
        w.writerows(list2d)
################################################################################

getHydroMonthlyGenValues()