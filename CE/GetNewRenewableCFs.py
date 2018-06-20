# Michael Craig
# October 4, 2016
# Determines hourly CFs for potential new wind & solar builds in a similar way
# as hourly CFs for existing wind & solar are determined, then trims
# those CFs to hours included in the CE model.

import operator, copy
from GetRenewableCFs import getRenewableCFs
from AuxFuncs import *


########### GET NEW WIND AND SOLAR CFS #########################################
# Determines hourly CFs for potential new wind & solar plants by getting average
# CF of input assumed incremental capacity of wind & solar.
# Inputs: gen fleet (2d list), net demand (1d list), list of state & state abbrevs,
# curr year of CE run, name of model (UC or CE)
# Outputs: hourly CFs for new wind and solar builds (1d list)
def getNewWindAndSolarCFs(genFleet, currYear, modelName, tzAnalysis, dataRoot, resultsDir,
                          windGenDataYr, currZone, fipsToZones, fipsToPolys):
    windCapacInCurrFleet = getPlantTypeCapacInFleet(genFleet, 'Wind')
    solarCapacInCurrFleet = getPlantTypeCapacInFleet(genFleet, 'Solar PV')

    (additionalWind, additionalSolar) = (5, 5)  # 3000, 3000

    totalWindCapac = windCapacInCurrFleet + additionalWind
    totalSolarCapac = solarCapacInCurrFleet + additionalSolar

    genFleetWithNewRE = addNewREToFleet(genFleet, totalWindCapac, totalSolarCapac, currZone)

    (windCFs, windCfsDtHr, windCfsDtSubhr, windIdAndCapac, solarCFs, solarCfsDtHr, solarCfsDtSubhr,
     solarFilenameAndCapac) = getRenewableCFs(genFleetWithNewRE, windCapacInCurrFleet,
                                              solarCapacInCurrFleet, tzAnalysis, dataRoot, windGenDataYr,
                                              currZone, fipsToZones, fipsToPolys)

    avgWindCF = getCapacWtdCF(rotate(windCfsDtHr), windIdAndCapac, tzAnalysis)
    avgWindCFSubhr = getCapacWtdCF(rotate(windCfsDtSubhr), windIdAndCapac, tzAnalysis)
    avgSolarCF = getCapacWtdCF(rotate(solarCfsDtHr), solarFilenameAndCapac, tzAnalysis)
    avgSolarCFSubhr = getCapacWtdCF(rotate(solarCfsDtSubhr), solarFilenameAndCapac, tzAnalysis)

    return (avgWindCF, avgWindCFSubhr, avgSolarCF, avgSolarCFSubhr, windCfsDtHr, windCfsDtSubhr,
            windIdAndCapac, solarCfsDtHr, solarCfsDtSubhr, solarFilenameAndCapac, additionalWind,
            additionalSolar)


# Inputs: gen fleet (2d list), plant type (str)
# Outputs: sum of capacity of all plants in fleet of that plant type
def getPlantTypeCapacInFleet(fleet, plantType):
    plantTypeCol = fleet[0].index('PlantType')
    capacCol = fleet[0].index('Capacity (MW)')
    capacsOfPlantType = [float(row[capacCol]) for row in fleet[1:] if row[plantTypeCol] == plantType]
    return sum(capacsOfPlantType)


# Temporarily add input fixed capacity of new wind and solar (jointly) to fleet.
# Adds all wind (new + existing) to single arbitrary state, since get
# wind and solar additions at best places across region.
# Inputs: gen fleet (2d list), capacity of wind and solar to add, list of states
# Output: gen fleet w/ wind & solar units added
def addNewREToFleet(genFleet, totalWindCapac, totalSolarCapac, currZone):
    genFleetWithNewRE = copy.deepcopy(genFleet)
    headers = genFleetWithNewRE[0]
    (zoneCol, plantTypeCol, fuelTypeCol, capacCol) = (headers.index('Region Name'),
                                                      headers.index('PlantType'), headers.index('Modeled Fuels'),
                                                      headers.index('Capacity (MW)'))
    plantTypes = ['Wind', 'Solar PV']
    # First remove existing wind and solar rows
    genFleetWithNewRE = [row for row in genFleetWithNewRE if row[plantTypeCol] not in plantTypes]
    # Now add all wind and solar (new + existing) to single state
    plantTypeToFuelType = {'Wind': 'Wind', 'Solar PV': 'Solar'}
    plantTypeToCapac = {'Wind': totalWindCapac, 'Solar PV': totalSolarCapac}
    emptyRow = ['' for x in range(len(headers))]
    for plantType in plantTypes:
        emptyRow[plantTypeCol] = plantType
        emptyRow[fuelTypeCol] = plantTypeToFuelType[plantType]
        emptyRow[capacCol] = plantTypeToCapac[plantType]
        rowToAdd = copy.deepcopy(emptyRow)
        rowToAdd[zoneCol] = currZone  # random state
        genFleetWithNewRE.append(rowToAdd)
    return genFleetWithNewRE


# Gets capacity-weighted CFs for set of renewable units. Used here to get
# capacity-weighted CF for input assumed incremental capacity of wind or solar.
# Inputs: hourly CFs for each RE unit (2d list), unit IDs and capacities (2d list)
# Outputs: hourly capacity-weighted CFs (1d list)
def getCapacWtdCF(allCfs, idAndCapacs, tzAnalysis):
    cfIdCol = allCfs[0].index('datetime' + tzAnalysis)
    (idAndCapacsIdCol, idAndCapacsCapacCol) = (idAndCapacs[0].index('Id'), idAndCapacs[0].index('FleetCapacity'))
    allCapacs = [row[idAndCapacsCapacCol] for row in idAndCapacs[1:]]
    allIds = [row[idAndCapacsIdCol] for row in idAndCapacs[1:]]
    totalCapac = sum(allCapacs)
    cfIds = [row[cfIdCol] for row in allCfs]
    capacWtCfs = [0] * len(allCfs[0][1:])
    for idx in range(1, len(cfIds)):
        cfId = cfIds[idx]
        idAndCapacRow = allIds.index(cfId)
        fleetCapac = allCapacs[idAndCapacRow]
        capacWt = fleetCapac / totalCapac
        cfs = allCfs[idx][1:]
        cfsWtd = [cf * capacWt for cf in cfs]
        capacWtCfs = list(map(operator.add, capacWtCfs, cfsWtd))
    return capacWtCfs


########### TRIM HOURS OF NEW RENEWABLE CFS TO HOURS OF CE #####################
# Inputs: hourly CFs for new wind & solar builds (1d lists), hours included in CE
# (1d list, 1-8760 basis)
# Outputs: zonal hourly CFs for new wind and solar builds only for CE hours
# as dict of zone: [CFs (1-8760 basis, 1d list)]
def trimNewRECFsToCEHours(zonalNewWindCFs, zonalNewSolarCFs, hoursForCE):
    newWindCFsCEZonal, newSolarCFsCEZonal = dict(), dict()
    for zone in zonalNewWindCFs:
        newWindCFsCEZonal[zone] = [zonalNewWindCFs[zone][hr - 1] for hr in
                                   hoursForCE]  # -1 b/c hours in year start @ 1, not 0 like Python idx
        newSolarCFsCEZonal[zone] = [zonalNewSolarCFs[zone][hr - 1] for hr in hoursForCE]  # -1 b/c
    return (newWindCFsCEZonal, newSolarCFsCEZonal)
