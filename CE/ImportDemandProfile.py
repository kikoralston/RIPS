#Michael Craig, 3 August 2016
#Import hourly demand profile for given set of states

#Takes regional hourly demand profile from IPM, and scales up/down total annual
#demand based on 2014 state-level total electricity demand data from EIA Form 861,
#which is scaled to 2015 data based on annual US electricity demand
#per Monthly Energy Review. Output is in MWh

import os,csv,copy
from AuxFuncs import *

#Outputs dictionary of zone:hourly zonal demand in MWh
def importZonalHourlyDemand(states,stateAbbrevs,ipmZones):
    #Get IPM hourly demand
    ipmHourlyProfiles = importIPMData()
    firstDataRow = findFirstDataRow(ipmHourlyProfiles,'Region')
    ipmHourlyProfiles = copy.deepcopy(ipmHourlyProfiles[firstDataRow:])
    zonalDemand = dict()
    for zone in ipmZones:
        zonalDemand[zone] = getIPMZonalHourlyDemand(ipmHourlyProfiles,zone)
    return zonalDemand

def getIPMZonalHourlyDemand(ipmHourlyProfiles,zone):
    (zoneCol,hour1Col,hour24Col) = (ipmHourlyProfiles[0].index('Region'),
                                  ipmHourlyProfiles[0].index('Hour 1'),
                                  ipmHourlyProfiles[0].index('Hour 24'))
    regionDemandRows = [row for row in ipmHourlyProfiles if row[zoneCol]==zone]
    ipmRegionHourlyDemand = []
    for row in regionDemandRows:
        ipmRegionHourlyDemand += row[hour1Col:hour24Col+1] 
    return [float(removeComma(val)) for val in ipmRegionHourlyDemand]

def importIPMData():
    currDir = 'C:\\Users\\mtcraig\\Desktop\\EPP Research\\Databases\\EPA IPM Clean Power Plan Files'
    currFilename = 'IPM513BaseCaseLoadDurationCurves_27Oct0214.csv'
    return readCSVto2dList(os.path.join(currDir,currFilename))

def removeComma(s):
    if ',' not in s: return s
    else:
        idx = s.find(',')
        return removeComma(s[:idx] + s[idx+1:])

def findFirstDataRow(list2d,keywordFirstDataRow):
    (notFoundFirstRow,ctr) = (True,0)
    while notFoundFirstRow == True:
        if keywordFirstDataRow in list2d[ctr]: notFoundFirstRow=False
        else: ctr+=1
    return ctr














################################################################################
####### MASTER FUNCTION ########################################################
################################################################################
#Outputs 1d list w/ hourly demand in MWh
def importHourlyDemand(testModel,states,stateAbbrevs):
    if testModel == False:
        statesToIPMRegions = {'North Carolina':'S_VACA','South Carolina':'S_VACA',
                            'Georgia':'S_SOU','Alabama':'S_SOU','Mississippi':'S_SOU',
                            'Tennessee':'S_C_TVA','Kentucky':'S_C_KY'}
        #Get IPM hourly demand
        ipmRegions = mapStatesToIPMRegions(states,statesToIPMRegions)
        ipmCombinedRegionHourlyDemand = getIPMHourlyDemand(ipmRegions) #MWh
        #Get national 2015 demand
        usDemandYears = (2014,2015)
        kWhToMWh = 1/1000
        (us2014Demand,us2015Demand) = getUSDemand(usDemandYears,kWhToMWh) #MWh
        #Get state-level 2014 demand
        stateDemandYear = 2014
        statesCombined2014Demand = getStatesDemand(states,stateAbbrevs,stateDemandYear) #MWh
        #Scale state-level 2014 demand to 2015
        scaledStatesCombined2014Demand = scaleStatesCombined2014Demand(statesCombined2014Demand,
                                                                    us2014Demand,us2015Demand)
        #Scale IPM hourly demand to state-level 2015 demand
        hourlyDemandProfile = scaleIPMHourlyDemandTo2015RegionTotal(ipmCombinedRegionHourlyDemand,
                                                                        scaledStatesCombined2014Demand)
    else: 
        hourlyDemandProfile = importTestDemandProfile()
    return hourlyDemandProfile

################################################################################
####### GET IPM HOURLY DEMAND ##################################################
################################################################################   
def mapStatesToIPMRegions(states,statesToIPMRegions):
    ipmRegions = set()
    for state in states: ipmRegions.add(statesToIPMRegions[state])
    return ipmRegions

def getIPMHourlyDemand(ipmRegions):
    ipmHourlyProfiles = importIPMData()
    firstDataRow = findFirstDataRow(ipmHourlyProfiles,'Region')
    ipmHourlyProfilesTrim = copy.deepcopy(ipmHourlyProfiles[firstDataRow:])
    ipmCombinedRegionHourlyDemand = None
    for region in ipmRegions:
        ipmRegionHourlyDemand = getIPMRegionHourlyDemand(region,ipmHourlyProfilesTrim)
        if ipmCombinedRegionHourlyDemand==None: 
            ipmCombinedRegionHourlyDemand = [float(removeComma(val)) for val in ipmRegionHourlyDemand]
        else: ipmCombinedRegionHourlyDemand = [ipmCombinedRegionHourlyDemand[idx] + 
                                               float(removeComma(ipmRegionHourlyDemand[idx]))
                                               for idx in range(len(ipmCombinedRegionHourlyDemand))]
    return ipmCombinedRegionHourlyDemand

def getIPMRegionHourlyDemand(region,ipmHourlyProfilesTrim):
    (regionCol,hour1Col,hour24Col) = (ipmHourlyProfilesTrim[0].index('Region'),
                                      ipmHourlyProfilesTrim[0].index('Hour 1'),
                                      ipmHourlyProfilesTrim[0].index('Hour 24'))
    regionDemandRows = [row for row in ipmHourlyProfilesTrim if row[regionCol]==region]
    ipmRegionHourlyDemand = []
    for row in regionDemandRows:
        ipmRegionHourlyDemand += row[hour1Col:hour24Col+1] 
    return ipmRegionHourlyDemand

def removeComma(s):
    if ',' not in s: return s
    else:
        idx = s.find(',')
        return removeComma(s[:idx] + s[idx+1:])

def findFirstDataRow(list2d,keywordFirstDataRow):
    (notFoundFirstRow,ctr) = (True,0)
    while notFoundFirstRow == True:
        if keywordFirstDataRow in list2d[ctr]: notFoundFirstRow=False
        else: ctr+=1
    return ctr

################################################################################
####### GET NATIONAL 2014 AND 2015 DEMAND ######################################
################################################################################   
def getUSDemand(usDemandYears,kWhToMWh):
    usAnnualDemand = importUSDemand()
    firstDataRow = findFirstDataRow(usAnnualDemand,'Annual Total')
    usAnnualDemandTrim = copy.deepcopy(usAnnualDemand[firstDataRow:])
    #Data is in million kWh; convert to MWh
    us2014Demand = getUSDemandForYear(usAnnualDemandTrim,usDemandYears[0])*kWhToMWh*1E6
    us2015Demand = getUSDemandForYear(usAnnualDemandTrim,usDemandYears[1])*kWhToMWh*1E6
    return (us2014Demand,us2015Demand)

def getUSDemandForYear(usAnnualDemandTrim,demandYear):
    yearCol = usAnnualDemandTrim[0].index('Annual Total')
    demandCol = usAnnualDemandTrim[0].index('Electricity End Use, Total')
    rowLabels = [row[yearCol] for row in usAnnualDemandTrim]
    return float(usAnnualDemandTrim[rowLabels.index(str(demandYear))][demandCol])

################################################################################
####### GET TOTAL 2014 STATE DEMANDS ###########################################
################################################################################   
def getStatesDemand(states,stateAbbrevs,stateDemandYear):
    statesDemandAllYears = importStatesDemand()
    firstDataRow = findFirstDataRow(statesDemandAllYears,'Year')
    statesDemandAllYearsTrim = copy.deepcopy(statesDemandAllYears[firstDataRow:])
    (yearCol,stateCol,demandCol) = (statesDemandAllYearsTrim[0].index('Year'),
                statesDemandAllYearsTrim[0].index('State'),statesDemandAllYearsTrim[0].index('Total'))
    states2014Demand = statesDemandAllYearsTrim[0] + [row for row in statesDemandAllYearsTrim if 
                                                        row[yearCol] == '2014']
    stateRowLabels = [row[stateCol] for row in states2014Demand]
    statesCombined2014Demand = 0
    for stateAbbrev in stateAbbrevs:
        stateRow = stateRowLabels.index(stateAbbrev)
        statesCombined2014Demand += float(removeComma(states2014Demand[stateRow][demandCol]))
    return statesCombined2014Demand

################################################################################
####### SCALE DEMAND TO 2015 ###################################################
################################################################################   
def scaleStatesCombined2014Demand(statesCombined2014Demand,us2014Demand,us2015Demand):
    usChange2014to2015 = (us2015Demand-us2014Demand)/us2014Demand
    return statesCombined2014Demand*(1+usChange2014to2015)

def scaleIPMHourlyDemandTo2015RegionTotal(ipmCombinedRegionHourlyDemand,scaledStatesCombined2014Demand):
    totalIPMDemand = sum(ipmCombinedRegionHourlyDemand)
    diffDemandIPMtoStates = (scaledStatesCombined2014Demand - totalIPMDemand)/totalIPMDemand
    return [val*(1+diffDemandIPMtoStates) for val in ipmCombinedRegionHourlyDemand]

################################################################################
####### IMPORT FUNCTIONS #######################################################
################################################################################
def importIPMData():
    currDir = 'C:\\Users\\mtcraig\\Desktop\\EPP Research\\Databases\\EPA IPM Clean Power Plan Files'
    currFilename = 'IPM513BaseCaseLoadDurationCurves_27Oct0214.csv'
    return readCSVto2dList(os.path.join(currDir,currFilename))

def importUSDemand():
    currDir = 'C:\\Users\\mtcraig\\Desktop\\EPP Research\\Databases\\EIA861andMER'
    currFilename = 'Table_7.6_Electricity_End_Use_annual.csv'
    return readCSVto2dList(os.path.join(currDir,currFilename))

def importStatesDemand():
    currDir = 'C:\\Users\\mtcraig\\Desktop\\EPP Research\\Databases\\EIA861andMER'
    currFilename = 'sales_annual_2016.csv'
    return readCSVto2dList(os.path.join(currDir,currFilename))

def importTestDemandProfile():
    currDir = 'C:\\Users\\mtcraig\\Desktop\\EPP Research\\Databases\\CETestFleet'
    currFilename = 'testDemandProfile.csv'
    return readCSVto2dList(os.path.join(currDir,currFilename))

################################################################################
####### TEST FUNCTIONS #########################################################
################################################################################

def testCall():
    states = ['North Carolina','South Carolina','Georgia','Mississippi','Alabama',
                    'Kentucky','Tennessee'] 
    stateAbbrevs = ['NC','SC','GA','MS','AL','KY','TN']
    scaledRegionalHourlyDemand = importHourlyDemand(states,stateAbbrevs)   
    #Spot checked values - see sales_annual_2016_comparedToIPM 

# testCall()