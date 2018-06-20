#Michael Craig
#October 4, 2016
#Functions for adding sets to GAMS database. Used for CE & UC models.

from GAMSAuxFuncs import *


def addGeneratorSets(db, genFleet):
    """ADD GENERATOR SETS

    Add gen sets & subsets

    :param db:
    :param genFleet:
    :return:
    """
    genSymbols = isolateGenSymbols(genFleet, '')
    (genSetName, genSetDescription, genSetDimension) = ('egu', 'existing generators', 1)
    genSet = addSet(db, genSymbols, genSetName, genSetDescription, genSetDimension)

    windGenSymbols = isolateGenSymbols(genFleet, 'Wind')
    (windGenSetName, windGenSetDescription, windGenSetDimension) = ('windegu', 'existing wind generators', 1)
    windGenSet = addSet(db, windGenSymbols, windGenSetName, windGenSetDescription, windGenSetDimension)

    solarGenSymbols = isolateGenSymbols(genFleet, 'Solar PV')
    (solarGenSetName, solarGenSetDescription, solarGenSetDimension) = ('solaregu', 'existing solar generators', 1)
    solarGenSet = addSet(db, solarGenSymbols, solarGenSetName, solarGenSetDescription, solarGenSetDimension)

    hydroGenSymbols = isolateGenSymbols(genFleet, 'Hydro')
    (hydroName, hydroDesc, hydroDim) = ('hydroegu', 'existing hydro generators', 1)
    hydroGenSet = addSet(db, hydroGenSymbols, hydroName, hydroDesc, hydroDim)

    pumpHydroGenSymbols = isolateGenSymbols(genFleet, 'Pumped Storage')
    (phname, phdesc, phdim) = ('pumphydroegu', 'existing hydro generators', 1)
    pumpHydroGenSet = addSet(db, pumpHydroGenSymbols, phname, phdesc, phdim)

    return genSet, genSymbols, hydroGenSet, hydroGenSymbols, pumpHydroGenSet, pumpHydroGenSymbols


def isolateGenSymbols(genFleet,genPlantType):
    """Gen symbols: g1, g2, etc., where # = row in gen fleet

    :param genFleet:
    :param genPlantType:
    :return:
    """
    if genPlantType == '': #all generators
        genSymbols = [createGenSymbol(row, genFleet[0]) for row in genFleet[1:]]
    else:
        plantTypeCol = genFleet[0].index('PlantType')
        genSymbols = [createGenSymbol(row, genFleet[0]) for row in genFleet[1:]
                      if row[plantTypeCol] == genPlantType]
    return genSymbols


def addHourSet(db, hours):
    """ADD HOUR SETS

    Add all hours

    :param db:
    :param hours:
    :return:
    """
    hourSymbols = [createHourSymbol(hour) for hour in hours]
    (hourSetName, hourSetDescrip, hourSetDim) = ('h', 'hour', 1)
    hourSet = addSet(db, hourSymbols, hourSetName, hourSetDescrip, hourSetDim)
    return (hourSet, hourSymbols)


def addHourSeasonSubsets(db, repHrsBySeason):
    """Define season subsets of hours

    Inputs: GAMS db, dict of (season:rep hrs)

    :param db:
    :param repHrsBySeason:
    """
    for season in repHrsBySeason:
        seasonHours = repHrsBySeason[season]
        seasonHourSymbols = [createHourSymbol(hour) for hour in seasonHours]
        hourSubsetName = createHourSubsetName(season)
        addSeasonHourSubset(db, hourSubsetName, season, seasonHourSymbols)


def createHourSubsetName(subsetPrefix):
    return subsetPrefix + 'h'


def addSeasonHourSubset(db, hourSubsetName, season, seasonHourSymbols):
    (seasonSetName, seasonSetDescrip, seasonSetDim) = (hourSubsetName, 'hours in ' + season, 1)
    hourSeasonSet = addSet(db, seasonHourSymbols, seasonSetName, seasonSetDescrip, seasonSetDim)


def addHourSpecialSubset(db, specialHrs):
    """Define special subsets of hours

    Inputs: GAMS db, 1d list of special hours

    :param db:
    :param specialHrs:
    """
    specialHourSymbols = [createHourSymbol(hour) for hour in specialHrs]
    hourSubsetName = createHourSubsetName('special')
    (specialSetName,specialSetDescrip,specialSetDim) = (hourSubsetName,'hours in special days',1)
    hourSpecialSet = addSet(db,specialHourSymbols,specialSetName,specialSetDescrip,specialSetDim)


def addPeakHourSubset(db, peakDemandHourZonal):
    """Define peak demand hour subset

    :param db:
    :param peakDemandHourZonal:
    :return:
    """
    peakHrSymbols = [createHourSymbol(peakDemandHourZonal[zone]) for zone in peakDemandHourZonal]

    # convert to set and back to list to keep only unique values
    peakHrSymbols = list(set(peakHrSymbols))

    hourSubsetName = createHourSubsetName('peak')
    (specialSetName, specialSetDescrip, specialSetDim) = (hourSubsetName, 'peak hour', 1)
    peakHourSet = addSet(db, peakHrSymbols, specialSetName, specialSetDescrip, specialSetDim)
    return peakHourSet, peakHrSymbols


def addZoneSets(db, ipmZoneNums):
    """ADD ZONE SETS

    :param db:
    :param ipmZoneNums:
    :return:
    """
    zoneSymbols = [createZoneSymbol(num) for num in ipmZoneNums]
    (zoneName, zoneDesc, zoneDim) = ('z', 'zones', 1)
    zoneSet = addSet(db, zoneSymbols, zoneName, zoneDesc, zoneDim)
    return zoneSet, zoneSymbols

def addLineSets(db, lines):
    """ADD LINE SETS

    :param db:
    :param lines:
    :return:
    """
    (lName, lDesc, lDim) = ('l', 'lines', 1)
    lineSet = addSet(db, lines, lName, lDesc, lDim)
    return lineSet


def addCellSet(db,cellsToZones):
    """ADD CELLS FOR NEW TECHS SET

    :param db:
    :param cellsToZones:
    :return:
    """
    (cellname,celldesc,celldim) = ('c', 'cells',1)
    cellset = addSet(db,[cell for cell in cellsToZones],cellname,celldesc,celldim)
    return cellset


def addNewTechsSets(db, newTechsCE, plantTypesCurtailed):
    """ADD NEW TECH SETS

    Inputs: GAMS db, new techs (2d list)

    :param db:
    :param newTechsCE:
    :param plantTypesCurtailed:
    :return:
    """
    (techSymbols, techCurtailedSymbols,
     renewTechSymbols, techNotCurtailedSymbols) = isolateTechSymbols(newTechsCE, plantTypesCurtailed)

    #Add all techs
    (techSetName, techSetDescrip, techSetDim) = ('tech', 'techs for expansion', 1)
    techSet = addSet(db, techSymbols, techSetName, techSetDescrip, techSetDim)

    #Add subset of techs that can be curtailed
    (techCurtName, techCurtDesc, techCurtDim) = ('techcurtailed', 'curtailed techs for expansion', 1)
    techCurtailedSet = addSet(db, techCurtailedSymbols, techCurtName, techCurtDesc, techCurtDim)

    #Add subset of techs that are wind or solar
    (renewTechSetName, renewTechSetDescrip, renewTechSetDim) = ('techrenew', 'renewable techs for expansion', 1)
    renewTechSet = addSet(db, renewTechSymbols, renewTechSetName, renewTechSetDescrip, renewTechSetDim)

    #Add all other techs
    (techNotCurtName, techNotCurtDes, techNotCurtDim) = ('technotcurtailed', 'techs not curtailed', 1)
    techNotCurtailedSet = addSet(db, techNotCurtailedSymbols, techNotCurtName, techNotCurtDes, techNotCurtDim)

    return (techSet, techSymbols, techCurtailedSet, techCurtailedSymbols, renewTechSet, renewTechSymbols,
            techNotCurtailedSet, techNotCurtailedSymbols)


def isolateTechSymbols(newTechsCE, ptCurtailed):
    """Takes in new techs (2d list), and returns tech types as: all types, renew (wind or solar) types, or types
    that can be curtailed

    :param newTechsCE:
    :param ptCurtailed:
    :return:
    """
    techSymbols = [createTechSymbol(row,newTechsCE[0], ptCurtailed) for row in newTechsCE[1:]]

    # Get RE as marked as RE
    renewCol = newTechsCE[0].index('ThermalOrRenewable')
    techRenewSymbols = [createTechSymbol(row, newTechsCE[0], ptCurtailed)
                        for row in newTechsCE[1:] if row[renewCol] == 'renewable']

    # Get curtailed as eligible to be curtailed
    techCol = newTechsCE[0].index('TechnologyType')
    techCurtailedSymbols = [createTechSymbol(row, newTechsCE[0], ptCurtailed) for row in newTechsCE[1:]
                            if row[techCol] in ptCurtailed]

    # Get techs not in renew or techcurtailed
    techNotCurtailedSymbols = [symb for symb in techSymbols if (symb not in (techRenewSymbols+techCurtailedSymbols))]

    return techSymbols, techCurtailedSymbols, techRenewSymbols, techNotCurtailedSymbols


def addSet(db, setSymbols, setName, setDescription, setDim):
    """GENERIC FUNCTION TO ADD SET TO DATABASE

    Adds set to GAMS db

    :param db:
    :param setSymbols:
    :param setName:
    :param setDescription:
    :param setDim:
    :return:
    """
    addedSet = db.add_set(setName, setDim, setDescription)
    for symbol in setSymbols:
        addedSet.add_record(symbol)

    return addedSet