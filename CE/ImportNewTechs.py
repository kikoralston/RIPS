# Michael Craig
# October 4, 2016
# Function imports data for new technologies eligible for construction in capacity expansion model

import os
from AuxFuncs import *


def getNewTechs(currYear, genparam, reserveparam):
    """Import new tech data
    Import new tech data by loading forecasts of tech data + 2d list w/ missing data, and fill in that 2d list w/
    forecast data for current year. Inputs used from the parameter objects are the flags indicating which
    techs to import or not import

    :param currYear: current year of analysis
    :param genparam: object of class General Parameters
    :param reserveparam: object of class Reserve Parameters
    :return: 2d list w/ headers
    """
    # Set directory data is in

    newPlantDataDir = os.path.join(genparam.dataRoot, 'NewPlantData', 'ATB')

    # Set which new tech file, based on whether in special scenario
    techFrmwrkFile, techDataFile = 'NewTechFrameworkWithCoolingTechs.csv', 'newTechDataATB.csv'
    # if scenario == 'solar': 
    #     techFrmwrkFile,techDataFile = 'NewTechFrameworkLowSolarCost.csv','newTechDataATBLowSolarCost.csv'
    # elif scenario == 'nuclear':
    #     techFrmwrkFile,techDataFile = 'NewTechFramework.csv','newTechDataATBLowNukeCost.csv'
    # else: 
    #     techFrmwrkFile,techDataFile = 'NewTechFramework.csv','newTechDataATB.csv'
    newTechsCEFilename = os.path.join(newPlantDataDir, techFrmwrkFile)
    newTechsCE = readCSVto2dList(newTechsCEFilename)

    # Filter out certain units
    if not genparam.allowCoalWithoutCCS: newTechsCE = removeCoalWithoutCCSFromNewTechs(newTechsCE)
    if genpara.onlyNSPSUnits: newTechsCE = removeUnitsNotCompliantWithNSPS(newTechsCE)
    if not genparam.permitOncethru: newTechsCE = removeOnceThroughUnits(newTechsCE)

    addRegUpOfferCostAndElig(newTechsCE, genparam.regUpCostCoeffs)
    inputValuesForCurrentYear(newTechsCE, newPlantDataDir, currYear, techDataFile)

    if genparam.incITC: modRECapCostForITC(newTechsCE, currYear)

    modCostsForCoolingTechs(newPlantDataDir, newTechsCE)

    return newTechsCE


# Remove coal w/out CCS from new tech data
# Input: new techs (2d list w/ headers)
# Outputs: new techs w/out coal w/out CCS (2d list)
def removeCoalWithoutCCSFromNewTechs(newTechsCE):
    techTypeCol = newTechsCE[0].index('TechnologyType')
    newTechsCENoCoalWithoutCCS = [row for row in newTechsCE if row[techTypeCol] != 'Coal Steam']
    return newTechsCENoCoalWithoutCCS


# Remove new techs not compliant w/ NSPS
# Input: new techs (2d list w/ headers)
# Outputs: new techs w/out units not compliant w/ NSPS (2d list w/ headers)
def removeUnitsNotCompliantWithNSPS(newTechsCE):
    nspsCompliantCol = newTechsCE[0].index('NSPSCompliant')
    newTechsNSPSCompliant = [newTechsCE[0]] + [row for row in newTechsCE[1:] if row[nspsCompliantCol] == 'Yes']
    return newTechsNSPSCompliant


# Remove units w/ once through cooling
def removeOnceThroughUnits(newTechsCE):
    coolTechCol = newTechsCE[0].index('Cooling Tech')
    newTechsCE = [newTechsCE[0]] + [row for row in newTechsCE[1:] if row[coolTechCol] != 'once through']


# Add reg offer cost and eligiblity to new techs
def addRegUpOfferCostAndElig(newTechsCE, regUpCostCoeffs):
    newTechsCE[0].extend(['RegOfferCost($/MW)', 'RegOfferElig'])
    plantTypeCol = newTechsCE[0].index('TechnologyType')
    for row in newTechsCE[1:]:
        if row[plantTypeCol] in regUpCostCoeffs:
            regCost, regOffer = regUpCostCoeffs[row[plantTypeCol]], 1
        else:
            regCost, regOffer = 0, 0
        row.extend([regCost, regOffer])


# Add values for HR, capex, FOM, and VOM based on current year
def inputValuesForCurrentYear(newTechsCE, newPlantDataDir, currYear, techDataFile):
    newPlantData = readCSVto2dList(os.path.join(newPlantDataDir, techDataFile))
    (forecastTechCol, forecastParamCol) = (newPlantData[0].index('TechnologyType'),
                                           newPlantData[0].index('Parameter'))
    newTechsTechCol = newTechsCE[0].index('TechnologyType')
    if str(currYear) in newPlantData[0]:
        yearCol = newPlantData[0].index(str(currYear))
    else:
        yearCol = len(newPlantData[0]) - 1
    for row in newPlantData[1:]:
        (currTech, currParam, currVal) = (row[forecastTechCol], row[forecastParamCol], row[yearCol])
        # newTechsTechs = [row[newTechsTechCol] for row in newTechsCE]
        if currTech in [row[newTechsTechCol] for row in newTechsCE]:
            rowIdxs = [idx for idx in range(len(newTechsCE)) if newTechsCE[idx][newTechsTechCol] == currTech]
            newTechsCECol = newTechsCE[0].index(currParam)
            for idx in rowIdxs: newTechsCE[idx][newTechsCECol] = float(currVal)


# Account for ITC in RE cap costs
# http://programs.dsireusa.org/system/program/detail/658
def modRECapCostForITC(newTechsCE, currYear):
    windItc, windItcYear = .21, 2020  # wind ITC expires at 2020; .21 is average of 2016-2019 ITCs
    solarItcInit, solarItcIndef, solarItcYear = .3, .1, 2020  # solar ITC doesn't expire, but goes from .3 to .1
    if currYear <= windItcYear: modRECost(newTechsCE, windItc, 'Wind')
    if currYear <= solarItcYear:
        modRECost(newTechsCE, solarItcInit, 'Solar PV')
    else:
        modRECost(newTechsCE, solarItcIndef, 'Solar PV')


def modRECost(newTechsCE, itc, plantType):
    ptCol = newTechsCE[0].index('TechnologyType')
    capexCol = newTechsCE[0].index('CAPEX(2012$/MW)')
    ptRow = [row[0] for row in newTechsCE].index(plantType)
    newTechsCE[ptRow][capexCol] *= (1 - itc)


# Modify costs for thermal plants based on cooling tech. Data from IECM.
# Inputs: new techs (2d list)
def modCostsForCoolingTechs(newPlantDataDir, newTechsCE):
    coolingCosts = readCSVto2dList(os.path.join(newPlantDataDir, 'CoolingTechCostData_IECM_2017.csv'))
    coolCoolCol, techCoolCol = coolingCosts[0].index('Cooling Tech'), coolingCosts[0].index('TechnologyType')
    capCoolCol = coolingCosts[0].index('Total capital requirement ($/kWnet)')
    fomCoolCol = coolingCosts[0].index('Fixed O&M ($/yr)')
    capacCoolCol = coolingCosts[0].index('Capacity(MW)')
    hrPenCoolCol = coolingCosts[0].index('Net heat rate penalty')
    coolTechCol, techTechCol = newTechsCE[0].index('Cooling Tech'), newTechsCE[0].index('TechnologyType')
    capexTechCol = newTechsCE[0].index('CAPEX(2012$/MW)')
    hrTechCol, fomTechCol = newTechsCE[0].index('HR(Btu/kWh)'), newTechsCE[0].index('FOM(2012$/MW/yr)')
    for row in newTechsCE[1:]:
        if row[coolTechCol] != 'NA' and row[coolTechCol] != 'once through':
            tech, cool = row[techTechCol], row[coolTechCol]
            costRow = [row for row in coolingCosts if row[coolCoolCol] == cool and row[techCoolCol] == tech][0]
            row[capexTechCol] += float(costRow[capCoolCol]) * 1000  # convert cost from $/kW to $/MW
            row[fomTechCol] += float(costRow[fomCoolCol]) / float(costRow[capacCoolCol])  # go from $/yr to $/MW/yr
            row[hrTechCol] *= (1 + float(costRow[hrPenCoolCol]))
