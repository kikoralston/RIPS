# Michael Craig
# October 4, 2016
# Combine all wind and solar units in fleet together (separately for
# each plant type), then remove all but combined unit from fleet.

from SetupFleet.SetupGeneratorFleet import *


# Inputs: gen fleet (2d list)
def combineWindAndSolarToSinglePlant(fleetUC, ipmZones, dataRoot):
    for zone in ipmZones:
        combineWindOrSolarPlants(fleetUC, zone, dataRoot, 'Wind', 'Wind')
        combineWindOrSolarPlants(fleetUC, zone, dataRoot, 'Solar', 'Solar PV')


# Adds new combined unit, then removes other units
# Inputs: gen fleet (2d list), fuel type to combine, plant type to combine
def combineWindOrSolarPlants(fleetUC, zone, dataRoot, fuelType, plantType):
    plantCol = fleetUC[0].index('PlantType')
    zoneCol = fleetUC[0].index('Region Name')
    rowIdxs = [idx for idx in range(len(fleetUC)) if (fleetUC[idx][plantCol] == plantType and
                                                      fleetUC[idx][zoneCol] == zone)]
    if len(rowIdxs) > 0:
        newRow = [''] * len(fleetUC[0])
        addParametersToNewWindOrSolarRow(fleetUC, newRow, rowIdxs, fuelType, plantType, zone, dataRoot)
        fleetUC.append(newRow)
        for idx in reversed(rowIdxs): fleetUC.pop(idx)


def addParametersToNewWindOrSolarRow(fleetUC, newRow, rowIdxs, fuelType, plantType, zone, dataRoot):
    """ Adds parameters to new wind or solar row

    :param fleetUC: gen fleet
    :param newRow: new gen row (fill values in)
    :param rowIdxs: row indices of units that are being combined into new gen row
    :param fuelType: fuel & plant type of units being combined.
    :param plantType: fuel & plant type of units being combined.
    :param zone:
    :param dataRoot:
    """
    addStateZoneOrisFuelOnlineYearAndPlantType(fleetUC, newRow, fuelType, plantType, zone)
    addRegEligAndCost(fleetUC, newRow, rowIdxs[0])
    (capacCol, hrCol) = (fleetUC[0].index('Capacity (MW)'), fleetUC[0].index('Heat Rate (Btu/kWh)'))
    fuelPriceCol = fleetUC[0].index('FuelPrice($/MMBtu)')
    totalCapac = sum([float(fleetUC[idx][capacCol]) for idx in rowIdxs])
    (newRow[capacCol], newRow[hrCol], newRow[fuelPriceCol]) = (totalCapac, 0, 0)
    (noxCol, so2Col, co2Col) = (fleetUC[0].index('NOxEmRate(lb/MMBtu)'),
                                fleetUC[0].index('SO2EmRate(lb/MMBtu)'),
                                fleetUC[0].index('CO2EmRate(lb/MMBtu)'))
    (newRow[noxCol], newRow[so2Col], newRow[co2Col]) = (0, 0, 0)

    # Fill in rand adder col w/ average of other rows
    randAdderCol = fleetUC[0].index('RandOpCostAdder($/MWh)')
    randAdders = [float(fleetUC[idx][randAdderCol]) for idx in rowIdxs]
    capacs = [float(fleetUC[idx][capacCol]) for idx in rowIdxs]
    capacFracs = [val / sum(capacs) for val in capacs]
    avgRandAdder = sum([randAdders[idx] * capacFracs[idx] for idx in range(len(randAdders))])
    newRow[randAdderCol] = avgRandAdder

    # Add UC, VOM & FOM parameters
    tempFleet = [fleetUC[0], newRow]
    ucHeaders = ['MinDownTime(hrs)', 'RampRate(MW/hr)', 'MinLoad(MW)', 'StartCost($)']
    addUCValues(tempFleet, ucHeaders, importPhorumData(dataRoot))
    addVomAndFomValues(tempFleet, importVomAndFomData(dataRoot))

    # fill in CO2 emission value (ton/GWh)
    co2TonGwhCol = fleetUC[0].index('CO2EmRate(ton/GWh)')
    newRow[co2TonGwhCol] = 0

    # Fill in random added error of emission col w/ average of other rows
    randAdderEmCol = fleetUC[0].index('CO2EmRandErr(ton/GWh)')
    randAdderEm = [float(fleetUC[idx][randAdderEmCol]) for idx in rowIdxs]
    capacs = [float(fleetUC[idx][capacCol]) for idx in rowIdxs]
    capacFracs = [val / sum(capacs) for val in capacs]
    avgRandAdder = sum([randAdderEm[idx] * capacFracs[idx] for idx in range(len(randAdderEm))])
    newRow[randAdderEmCol] = avgRandAdder


# Copy down reg eligibility & cost from first wind and solar row
def addRegEligAndCost(fleetUC, newRow, firstOtherRow):
    regElig, regCost = fleetUC[0].index('RegOfferElig'), fleetUC[0].index('RegOfferCost($/MW)')
    newRow[regElig] = fleetUC[firstOtherRow][regElig]
    newRow[regCost] = fleetUC[firstOtherRow][regCost]
