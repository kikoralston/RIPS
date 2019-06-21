#Michael Craig
#October 4, 2016
#Functions for adding sets to GAMS database. Used for CE & UC models.

from GAMSAuxFuncs import *


def addGeneratorSets(db, genFleet):
    """ADD GENERATOR SETS

    Add gen sets & subsets

    :param db: gams database object
    :param genFleet: 2d list ith generator fleet data
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


def isolateGenSymbols(genFleet, genPlantType):
    """Gen symbols: g1, g2, etc., where # = row in gen fleet

    :param genFleet: 2d list with generator fleet data
    :param genPlantType: string with PlantType
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

    :param db: GAMS database object
    :param hours: 1d list with hours to be included on set
    :return:
    """
    hourSymbols = [createHourSymbol(hour) for hour in hours]
    (hourSetName, hourSetDescrip, hourSetDim) = ('h', 'hour', 1)
    hourSet = addSet(db, hourSymbols, hourSetName, hourSetDescrip, hourSetDim)

    return hourSet, hourSymbols


def addHourSeasonSubsets(db, repHrsBySeason):
    """Define season subsets of hours

    Inputs: GAMS db, dict of (season:rep hrs)

    :param db: GAMS database object
    :param repHrsBySeason: dict of (season:rep hrs)
    """
    seasonList = repHrsBySeason[list(repHrsBySeason.keys())[0]].keys()

    for season in seasonList:
        hourSubsetName = createHourSubsetName(season)
        (seasonSetName, seasonSetDescrip, seasonSetDim) = (hourSubsetName, 'hours in ' + season, 2)
        seasonHourSymbols = []
        for gcm in repHrsBySeason:
            seasonHours = repHrsBySeason[gcm][season]
            seasonHourSymbols = seasonHourSymbols + [[gcm, createHourSymbol(hour)] for hour in seasonHours]

        # creates a 2d set for this season
        hourSeasonSet = addSet(db, seasonHourSymbols, seasonSetName, seasonSetDescrip, seasonSetDim)


def createHourSubsetName(subsetPrefix):
    return subsetPrefix + 'h'


def addHourSpecialSubset(db, specialHrs):
    """Define special subsets of hours

    Inputs: GAMS db, 1d list of special hours

    :param db: GAMS database object
    :param specialHrs: 1d list of special hours
    """
    specialHourSymbols = [[gcm, createHourSymbol(hour)] for gcm in specialHrs for hour in specialHrs[gcm]]
    hourSubsetName = createHourSubsetName('special')
    (specialSetName, specialSetDescrip, specialSetDim) = (hourSubsetName, 'hours in special days', 2)

    hourSpecialSet = addSet(db, specialHourSymbols, specialSetName, specialSetDescrip, specialSetDim)


def addPeakHourSubset(db, peakDemandHourZonal, genparam):
    """Define peak demand hour subset

    :param db: GAMS database object
    :param peakDemandHourZonal:
    :return:
    """
    # create new dict which converts zone names (e.g. S_C_TVA) to zone symbols (e.g, z1) and gets hour symbol
    auxDict = dict()
    for gcm in peakDemandHourZonal:
        auxDict2 = dict()
        for z in peakDemandHourZonal[gcm]:
            zoneSymbol = createZoneSymbol(genparam.ipmZoneNums[genparam.ipmZones.index(z)])
            auxDict2[zoneSymbol] = createHourSymbol(peakDemandHourZonal[gcm][z])

        auxDict[gcm] = auxDict2

    peakHrSymbols = [[gcm, zone, auxDict[gcm][zone]] for gcm in auxDict for zone in auxDict[gcm]]

    # convert to set and back to list to keep only unique values
    # peakHrSymbols = list(set(peakHrSymbols))

    hourSubsetName = createHourSubsetName('peak')
    (peakHrSetName, peakHrSetDescrip, peakHrSetDim) = (hourSubsetName, 'mapping of gcm of peak hour in each zone', 3)
    peakHourSet = addSet(db, peakHrSymbols, peakHrSetName, peakHrSetDescrip, peakHrSetDim)

    return peakHourSet, peakHrSymbols


def addZoneSets(db, ipmZoneNums):
    """ADD ZONE SETS

    :param db: GAMS database object
    :param ipmZoneNums:
    :return:
    """
    zoneSymbols = [createZoneSymbol(num) for num in ipmZoneNums]
    (zoneName, zoneDesc, zoneDim) = ('z', 'zones', 1)
    zoneSet = addSet(db, zoneSymbols, zoneName, zoneDesc, zoneDim)
    return zoneSet, zoneSymbols


def addLineSets(db, lines):
    """ADD LINE SETS

    :param db: GAMS database object
    :param lines:
    :return:
    """
    (lName, lDesc, lDim) = ('l', 'lines', 1)
    lineSet = addSet(db, lines, lName, lDesc, lDim)
    return lineSet


def addCellSet(db,cellsToZones):
    """ADD CELLS FOR NEW TECHS SET

    :param db: GAMS database object
    :param cellsToZones:
    :return:
    """
    (cellname,celldesc,celldim) = ('c', 'cells',1)
    cellset = addSet(db,[cell for cell in cellsToZones],cellname,celldesc,celldim)
    return cellset


def addNewTechsSets(db, newTechsCE, plantTypesCurtailed, blocksWind=None, blocksSolar=None):
    """ADD NEW TECH SETS

    Inputs: GAMS db, new techs (2d list)

    :param db: GAMS database object
    :param newTechsCE:
    :param plantTypesCurtailed:
    :return:
    """
    (techSymbols, techCurtailedSymbols,
     renewTechSymbols, techNotCurtailedSymbols) = isolateTechSymbols(newTechsCE, plantTypesCurtailed, blocksWind,
                                                                     blocksSolar)

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


def isolateTechSymbols(newTechsCE, ptCurtailed, blocksWind=None, blocksSolar=None):
    """Takes in new techs (2d list), and returns tech types as: all types, renew (wind or solar) types, or types
    that can be curtailed

    :param newTechsCE: 2d list with new techs data
    :param ptCurtailed: 1-d list with plant types that are curtailed
    :param blocksWind: 1-d list with indexes of blocks of wind
    :param blocksSolar: 1-d list with indexes of blocks of solar
    :return:
    """
    techSymbols = [createTechSymbol(row, newTechsCE[0], ptCurtailed) for row in newTechsCE[1:]]

    # create different wind and solar tech types according to aggregated block

    if blocksWind is not None:
        techre = techSymbols.pop(techSymbols.index('Wind'))
        techSymbols = techSymbols + blocksWind

    if blocksSolar is not None:
        techre = techSymbols.pop(techSymbols.index('Solar PV'))
        techSymbols = techSymbols + blocksSolar

    # Get RE as marked as RE
    #renewCol = newTechsCE[0].index('ThermalOrRenewable')
    #techRenewSymbols = [createTechSymbol(row, newTechsCE[0], ptCurtailed)
    #                    for row in newTechsCE[1:] if row[renewCol] == 'renewable']
    techRenewSymbols = [t for t in techSymbols if 'Solar PV' in t or 'Wind' in t]

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

    :param db: GAMS database object
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