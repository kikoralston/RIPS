# Michael Craig
# October 4, 2016
# Determines hourly CFs for potential new wind & solar builds in a similar way
# as hourly CFs for existing wind & solar are determined, then trims
# those CFs to hours included in the CE model.

import operator, copy
from GetRenewableCFs import getRenewableCFs
from AuxFuncs import *


def getNewWindAndSolarCFs(genFleet, currYear, modelName, tzAnalysis, dataRoot, resultsDir, windGenDataYr, currZone,
                          fipsToZones, fipsToPolys, ncores_py=1):
    """ GET NEW WIND AND SOLAR CFS

    Determines hourly CFs for potential new wind & solar plants by getting average CF of input assumed incremental
    capacity of wind & solar.

    :param genFleet: gen fleet (2d list), net demand (1d list), list of state & state abbrevs,
    :param currYear: curr year of CE run,
    :param modelName: name of model (UC or CE)
    :param tzAnalysis:
    :param dataRoot:
    :param resultsDir:
    :param windGenDataYr:
    :param currZone: current zone
    :param fipsToZones:
    :param fipsToPolys:
    :return: hourly CFs for new wind and solar builds (1d list)
    """
    windCapacInCurrFleet = getPlantTypeCapacInFleet(genFleet, 'Wind')
    solarCapacInCurrFleet = getPlantTypeCapacInFleet(genFleet, 'Solar PV')

    (additionalWind, additionalSolar) = (5, 5)  # 3000, 3000

    totalWindCapac = windCapacInCurrFleet + additionalWind
    totalSolarCapac = solarCapacInCurrFleet + additionalSolar

    genFleetWithNewRE = addNewREToFleet(genFleet, totalWindCapac, totalSolarCapac, currZone)

    (windCFs, windCfsDtHr, windCfsDtSubhr, windIdAndCapac, solarCFs, solarCfsDtHr, solarCfsDtSubhr,
     solarFilenameAndCapac) = getRenewableCFs(genFleetWithNewRE, windCapacInCurrFleet, solarCapacInCurrFleet,
                                              tzAnalysis, dataRoot, windGenDataYr, currZone, fipsToZones,
                                              fipsToPolys, ncores_py=ncores_py)

    avgWindCF = getCapacWtdCF(rotate(windCfsDtHr), windIdAndCapac, tzAnalysis)
    avgWindCFSubhr = getCapacWtdCF(rotate(windCfsDtSubhr), windIdAndCapac, tzAnalysis)
    avgSolarCF = getCapacWtdCF(rotate(solarCfsDtHr), solarFilenameAndCapac, tzAnalysis)
    avgSolarCFSubhr = getCapacWtdCF(rotate(solarCfsDtSubhr), solarFilenameAndCapac, tzAnalysis)

    return (avgWindCF, avgWindCFSubhr, avgSolarCF, avgSolarCFSubhr, windCfsDtHr, windCfsDtSubhr,
            windIdAndCapac, solarCfsDtHr, solarCfsDtSubhr, solarFilenameAndCapac, additionalWind,
            additionalSolar)


def getPlantTypeCapacInFleet(fleet, plantType):
    """

    :param fleet: gen fleet (2d list),
    :param plantType: plant type (str)
    :return: sum of capacity of all plants in fleet of that plant type
    """
    plantTypeCol = fleet[0].index('PlantType')
    capacCol = fleet[0].index('Capacity (MW)')
    capacsOfPlantType = [float(row[capacCol]) for row in fleet[1:] if row[plantTypeCol] == plantType]
    return sum(capacsOfPlantType)


def addNewREToFleet(genFleet, totalWindCapac, totalSolarCapac, currZone):
    """ Temporarily add input fixed capacity of new wind and solar (jointly) to fleet.

    Temporarily add input fixed capacity of new wind and solar (jointly) to fleet. Adds all wind (new + existing)
    to single arbitrary state, since get wind and solar additions at best places across region.

    :param genFleet: gen fleet (2d list),
    :param totalWindCapac: capacity of wind to add
    :param totalSolarCapac: capacity of solar to add
    :param currZone: current zone
    :return: gen fleet w/ wind & solar units added
    """
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


def getCapacWtdCF(allCfs, idAndCapacs, tzAnalysis):
    """ Gets capacity-weighted CFs for set of renewable units.

    Gets capacity-weighted CFs for set of renewable units. Used here to get capacity-weighted CF for input
    assumed incremental capacity of wind or solar.

    :param allCfs: hourly CFs for each RE unit (2d list),
    :param idAndCapacs: unit IDs and capacities (2d list)
    :param tzAnalysis: timezone of analysis
    :return: hourly capacity-weighted CFs (1d list)
    """
    capacWtCfs = None

    if allCfs is not None:
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


def trimNewRECFsToCEHours(zonalNewWindCFs, zonalNewSolarCFs, hoursForCE):
    """TRIM HOURS OF NEW RENEWABLE CFS TO HOURS OF CE

    :param zonalNewWindCFs: hourly CFs for new wind & solar builds (1d lists)
    :param zonalNewSolarCFs: hourly CFs for new wind & solar builds (1d lists)
    :param hoursForCE: hours included in CE (1d list, 1-8760 basis)
    :return: zonal hourly CFs for new wind and solar builds only for CE hours as dict of zone:
            [CFs (1-8760 basis, 1d list)]
    """
    newWindCFsCEZonal, newSolarCFsCEZonal = dict(), dict()

    for gcm in hoursForCE:
            # (hr - 1) b/c hours in year start @ 1, not 0 like Python idx
        newWindCFsCEZonal[gcm] = {zone: [zonalNewWindCFs[zone][hr - 1] for hr in hoursForCE[gcm]]
                                  for zone in zonalNewWindCFs}
        newSolarCFsCEZonal[gcm] = {zone: [zonalNewSolarCFs[zone][hr - 1] for hr in hoursForCE[gcm]]
                                   for zone in zonalNewWindCFs}

    return newWindCFsCEZonal, newSolarCFsCEZonal
