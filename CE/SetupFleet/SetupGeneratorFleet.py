# Michael Craig, 6 June 2016

import os, copy, operator, random
from AuxFuncs import *
from thermalderatings.CurtailmentRegressions import getCoolType


def setupGeneratorFleet(currYear, genparam, reserveparam):
    """Imports fleet data

    This is the main function to compile the existing power plant fleet that will be used in the simulation

    Uses pre defined parameters of this analysis (states, year, & power system for analysis) to import data
    of generator fleet. Imports fleet, isolates fleet to given state and power system, removes retired units, and
    adds emissions rates and cooling information.

    Constructs base generator fleet w/ cooling information, emissions rates, VOM, fuel price, and unit commitment
    parameters. Then the fleet is compressed.

    DATA SOURCES: fleet - NEEDS; cooling information - EIA860; emissions rates - eGRID; VOM & UC parameters - PHORUM.

    :param currYear: (int) year of fleet
    :param genparam: object of class :mod:`Generalparameters`
    :param reserveparam: object of class :mod:`Reserveparameters`
    :return: (2d list) generator fleet data
    """
    # Import entire current fleet

    baseGenFleet = importNEEDSFleet(genparam.dataRoot)

    # Slim down fleet based on region & year of analysis
    stateColName = "State Name"
    isolateGensInStates(baseGenFleet, genparam.states, stateColName)
    isolateGensInPowerSystem(baseGenFleet, genparam.ipmZones)
    # removeRetiredUnits(baseGenFleet,retirementYearScreen)
    # Modify / add fleet parameters
    addEmissionsRates(baseGenFleet, genparam.states, genparam.dataRoot)
    addCoolingTechnologyAndSource(baseGenFleet, genparam.states, genparam.dataRoot)
    addLatLong(baseGenFleet, genparam.states, genparam.dataRoot)
    # Compress fleet to get rid of tiny units
    if genparam.compressFleet:
        # write2dListToCSV(baseGenFleet,'genFleetPreCompression.csv')
        baseGenFleet = performFleetCompression(baseGenFleet, genparam.ipmZones, genparam.ptCurtailedAll)
    # Add VOM & FOM values
    vomAndFomData = importVomAndFomData(genparam.dataRoot)
    addVOMandFOM(baseGenFleet, vomAndFomData)
    # Add PHORUM-based UC parameters
    phorumData = importPhorumData(genparam.dataRoot)
    addUnitCommitmentParameters(baseGenFleet, phorumData)
    # Add fuel prices
    addFuelPrices(baseGenFleet, currYear, genparam.fuelPricesTimeSeries)
    # Add random value to rows that will be included in op cost
    addRandomOpCostAdder(baseGenFleet, genparam.ocAdderMin, genparam.ocAdderMax)
    # Add reg offer costs and reg offer eligibility
    addRegResOfferAndElig(baseGenFleet, reserveparam.regUpCostCoeffs)

    # add standard coolDesignT using same source as new plants
    newPlantDataDir = os.path.join(genparam.dataRoot, 'NewPlantData', 'ATB')
    coolingCosts = readCSVto2dList(os.path.join(newPlantDataDir, 'CoolingTechCostData_IECM_2017.csv'))

    coolCoolCol, techCoolCol = coolingCosts[0].index('Cooling Tech'), coolingCosts[0].index('TechnologyType')
    coolDesignTcol = coolingCosts[0].index('coolingDesignT')

    coolTechCol, plantTypeCol = baseGenFleet[0].index('Cooling Tech'), baseGenFleet[0].index('PlantType')

    # add column with standard cooling design temperature for each tech
    baseGenFleet[0] = baseGenFleet[0] + ['coolingDesignT']

    for row in baseGenFleet[1:]:
        tech, cool = row[plantTypeCol], row[coolTechCol]

        # change description of cooling type to the one compatible with new techs and regressions
        cool = getCoolType(cool)

        # finds row in data base of new techs that matches tech and cool
        coolingRow = [rr for rr in coolingCosts if rr[coolCoolCol] == cool and rr[techCoolCol] == tech]

        if len(coolingRow) > 0:
            coolingRow = coolingRow[0]
            row.append(coolingRow[coolDesignTcol])
        else:
            row.append('NA')

    # add columns with CO2 emissions in ton/GWh and random error over original value
    baseGenFleet[0] = baseGenFleet[0] + ['CO2EmRate(ton/GWh)'] + ['CO2EmRandErr(ton/GWh)']

    co2LbCol = baseGenFleet[0].index('CO2EmRate(lb/MMBtu)')
    hrCol = baseGenFleet[0].index('Heat Rate (Btu/kWh)')

    for row in baseGenFleet[1:]:
        co2ratevalue, hrvalue = float(row[co2LbCol]), float(row[hrCol])

        co2rateTonvalue = (co2ratevalue / 2000) * hrvalue
        co2emRandomError = random.uniform(0, 0.05)

        row.append(co2rateTonvalue)
        row.append(co2emRandomError)

    return baseGenFleet


def addCoolingTechnologyAndSource(baseGenFleet, statesForAnalysis, dataRoot):
    """ADD COOLING TECHNOLOGY AND SOURCE TO FLEET

    Imports cooling map and data from EIA 860, then adds them generator fleet

    :param baseGenFleet: (2d list) generator fleet (see :py:func:`.importNEEDSFleet`)
    :param statesForAnalysis: states for analysis (1d list)
    :param dataRoot: (string) path to root of data folder
    """
    [assocData, equipData] = get860Data(statesForAnalysis, dataRoot)
    unitsToCoolingIDMap = mapUnitsToCoolingID(assocData, baseGenFleet, equipData)
    addCoolingInfoToFleet(baseGenFleet, equipData, unitsToCoolingIDMap)


def addCoolingInfoToFleet(baseGenFleet, equipData, unitsToCoolingIDMap):
    """Adds cooling technology and source to fleet

    :param baseGenFleet: (2d list) generator fleet (see :func:`.importNEEDSFleet`)
    :param equipData: (2d list) cooling equipment data (EIA 860)
    :param unitsToCoolingIDMap: (dictionary) map of generator ID to cooling ID
    """
    coolingHeaders = ['Cooling Tech', 'Cooling Source']
    baseGenFleet[0].extend(coolingHeaders)
    fleetHeadersMap = mapHeadersToCols(baseGenFleet)
    baseORISCol = fleetHeadersMap['ORIS Plant Code']
    baseUnitCol = fleetHeadersMap['Unit ID']
    for idx in range(1, len(baseGenFleet)):
        genRow = baseGenFleet[idx]
        (orisID, genID) = (genRow[baseORISCol], genRow[baseUnitCol])
        coolingID = unitsToCoolingIDMap[orisID + '+' + genID]
        if coolingID != 'NoMatch':
            [coolingTech, coolingSource] = getCoolingTechAndSource(equipData, orisID, coolingID)
        else:
            [coolingTech, coolingSource] = (coolingID, coolingID)
        baseGenFleet[idx].extend([coolingTech, coolingSource])


def getCoolingTechAndSource(equipData, orisID, coolingID):
    """Gets cooling tech for given ORIS and cooling ID. Works for generators in NEEDS that have a corresponding
    unit in EIA860.

    :param equipData: (2d list) cooling tech data
    :param orisID: (str) oris ID
    :param coolingID: (str) cooling ID
    :return: tuple with cooling technology and source for input generator
    """
    equipHeadersMap = mapHeadersToCols(equipData)
    equipORISCol = equipHeadersMap['Plant Code']
    equipCoolingIDCol = equipHeadersMap['Cooling ID']
    equipCoolingTechCol = equipHeadersMap['Cooling Type 1']
    equipCoolingSourceCol = equipHeadersMap['Cooling Water Source']
    (equipORISIDs, equipCoolingIDs) = (colTo1dList(equipData, equipORISCol),
                                       colTo1dList(equipData, equipCoolingIDCol))
    equipRow = search2Lists(equipORISIDs, equipCoolingIDs, orisID, coolingID)
    coolingTech = equipData[equipRow][equipCoolingTechCol]
    coolingTechMap = getCoolingTechMap()
    coolingTech = coolingTechMap[coolingTech]
    if coolingTech == 'dry cooling':
        coolingSource = 'None'
    else:
        coolingSource = equipData[equipRow][equipCoolingSourceCol]
    retireCol = equipHeadersMap['Cooling Status']
    return [coolingTech, coolingSource]


def getCoolingTechMap():
    """Maps 2-letter codes to comprehensible cooling techs

    :return: dictionary mapping cooling tech abbrev to full name (e.g. {'DC': 'dry cooling'})
    """
    coolingTechMap = {'DC': 'dry cooling', 'OC': 'once through with pond',
                      'ON': 'once through no pond', 'RC': 'recirculating with pond or canal',
                      'RF': 'recirculating with tower', 'RI': 'recirculating with tower',
                      'RN': 'recirculating with tower', 'HT': 'helper tower',
                      'OT': 'other', 'HRC': 'hybrid pond or canal with dry cooling',
                      'HRF': 'hybrid tower with dry cooling',
                      'HRI': 'hybrid tower with dry cooling'}
    return coolingTechMap


def mapUnitsToCoolingID(assocData, baseGenFleet, equipData):
    """Maps NEEDS generator IDs to EIA860 cooling IDs

    :param assocData: (2d list) cooling association data from EIA860
    :param baseGenFleet: (2d list) base generator fleet from NEEDS (see :func:`.importNEEDSFleet`)
    :param equipData: (2d list) map of NEEDS generators to EIA860 cooling ID or 'NoMatch'
    :return: dictionary mapping NEEDS generator IDs to EIA860 cooling IDs
    """
    assocHeadersMap = mapHeadersToCols(assocData)
    fleetHeadersMap = mapHeadersToCols(baseGenFleet)
    assocORISCol = assocHeadersMap['Plant Code']
    assocBoilerCol = assocHeadersMap['Boiler ID']
    assocCoolingIDCol = assocHeadersMap['Cooling ID']
    baseORISCol = fleetHeadersMap['ORIS Plant Code']
    baseUnitCol = fleetHeadersMap['Unit ID']
    (assocORISIDs, assocBlrIDs) = (colTo1dList(assocData, assocORISCol),
                                   colTo1dList(assocData, assocBoilerCol))
    mapBaseGensToCoolingID = dict()
    for idx in range(1, len(baseGenFleet)):
        (baseORIS, baseUnit) = (baseGenFleet[idx][baseORISCol],
                                baseGenFleet[idx][baseUnitCol])
        assocRow = search2Lists(assocORISIDs, assocBlrIDs, baseORIS, baseUnit)
        if assocRow == False:  # no matching oris-gen row
            if baseORIS not in assocORISIDs:  # no matching oris row
                assocCoolingID = 'NoMatch'
            else:  # matching oris row
                assocRow = assocORISIDs.index(baseORIS)
                assocRowLinkedToRetiredUnit = getRetirementStatusOfAssocRow(assocRow,
                                                                            assocData, equipData)
                if assocRowLinkedToRetiredUnit == False:
                    assocCoolingID = assocData[assocRow][assocCoolingIDCol]
                else:
                    foundMatchingRow = False
                    while assocRowLinkedToRetiredUnit == True:
                        restOfAssocOris = assocORISIDs[assocRow + 1:]
                        if baseORIS in restOfAssocOris:
                            assocRow += restOfAssocOris.index(baseORIS) + 1
                            assocRowLinkedToRetiredUnit = getRetirementStatusOfAssocRow(assocRow,
                                                                                        assocData, equipData)
                            if assocRowLinkedToRetiredUnit == False:
                                foundMatchingRow = True
                        else:
                            assocRowLinkedToRetiredUnit = False
                    if foundMatchingRow:
                        assocCoolingID = assocData[assocRow][assocCoolingIDCol]
                    else:
                        assocCoolingID = 'NoMatch'
        else:  # found matching oris-gen row
            # assocCoolingID = assocData[assocRow][assocCoolingIDCol]
            assocRowLinkedToRetiredUnit = getRetirementStatusOfAssocRow(assocRow,
                                                                        assocData, equipData)
            if assocRowLinkedToRetiredUnit == False:  # not retired, so done
                assocCoolingID = assocData[assocRow][assocCoolingIDCol]
            else:  # retired, so keep looking
                (restOfAssocOris, restOfAssocBlrs) = (assocORISIDs[assocRow + 1:], assocBlrIDs[assocRow + 1:])
                assocRow += search2Lists(restOfAssocOris, restOfAssocBlrs, baseORIS, baseUnit) + 1
                if assocRow == False:  # no more oris-gen row matches
                    assocCoolingID = 'NoMatch'
                else:  # found another oris-gen row match
                    assocRowLinkedToRetiredUnit = getRetirementStatusOfAssocRow(assocRow,
                                                                                assocData, equipData)
                    if assocRowLinkedToRetiredUnit:  # found antoher retired cooling type - so quit
                        assocCoolingID = 'NoMatch'
                    else:  # found non-retired match
                        assocCoolingID = assocData[assocRow][assocCoolingIDCol]
        mapBaseGensToCoolingID[baseORIS + '+' + baseUnit] = assocCoolingID
    return mapBaseGensToCoolingID


def getRetirementStatusOfAssocRow(assocRow, assocData, equipData):
    """Check if cooling tech associated with cooling ID for ORIS ID match is retired. Returns true if unit is retired.

    :param assocRow:
    :param assocData:
    :param equipData:
    :return: (boolean) True if unit is retired. False otherwise
    """
    assocHeadersMap = mapHeadersToCols(assocData)
    assocOrisIDCol = assocHeadersMap['Plant Code']
    assocCoolingIDCol = assocHeadersMap['Cooling ID']
    (orisID, coolingID) = (assocData[assocRow][assocOrisIDCol],
                           assocData[assocRow][assocCoolingIDCol])
    equipHeadersMap = mapHeadersToCols(equipData)
    equipORISCol = equipHeadersMap['Plant Code']
    equipCoolingIDCol = equipHeadersMap['Cooling ID']
    equipRetiredCol = equipHeadersMap['Cooling Status']
    (equipORISIDs, equipCoolingIDs) = (colTo1dList(equipData, equipORISCol),
                                       colTo1dList(equipData, equipCoolingIDCol))
    equipRow = search2Lists(equipORISIDs, equipCoolingIDs, orisID, coolingID)
    retiredStatus = equipData[equipRow][equipRetiredCol]

    return retiredStatus == 'RE'


def get860Data(statesForAnalysis, dataRoot):
    """Imports EIA860 cooling equipment and association data, and isolates equipment data to units in states of analysis

    :param statesForAnalysis: (1d list) states for analysis
    :param dataRoot: (string) path to directory with input data
    :return: (2d list) EIA860 cooling IDs and technologies
    """
    [assocData, equipData] = import860data(dataRoot)
    # First row is useless data in both lists - remove it
    assocData.pop(0)
    equipData.pop(0)
    stateColName = 'State'
    statesForAnalysisAbbrev = getStateAbbrevs(statesForAnalysis)
    isolateGensInStates(equipData, statesForAnalysisAbbrev, stateColName)
    return [assocData, equipData]


def import860data(dataRoot):
    """Imports 860 equipment and association data

    :param dataRoot: (string) path to directory with input data
    :return: (2d list) EIA860 cooling association and equipment data
    """
    dir860 = os.path.join(dataRoot, 'EIA8602014')

    assocName = '6_1_EnviroAssoc_Y2014_cooling.csv'
    equipName = '6_2_EnviroEquip_Y2014_cooling.csv'
    assocData = readCSVto2dList(os.path.join(dir860, assocName))
    equipData = readCSVto2dList(os.path.join(dir860, equipName))
    return [assocData, equipData]


def addEmissionsRates(baseGenFleet, statesForAnalysis, dataRoot):
    """Add emission rates from EGrid to generator fleet

    Adds a column with eGRID emissions rates to the 2d list with generator fleet. It modifies the 2d list
    `baseGenFleet`. This function does not return anything.

    :param baseGenFleet: (2d list) generator fleet (see :func:`.importNEEDSFleet`)
    :param statesForAnalysis: (1d list) states for analysis
    :param dataRoot: (string) path to directory with input data
    """
    (egridBoiler, egridPlant) = importeGridData(statesForAnalysis, dataRoot)
    emsHeadersToAdd = ["NOxEmRate(lb/MMBtu)", "SO2EmRate(lb/MMBtu)",
                       "CO2EmRate(lb/MMBtu)"]
    addHeaders(baseGenFleet, emsHeadersToAdd)
    addEmissionsRatesValues(baseGenFleet, egridBoiler, egridPlant)
    # write2dListToCSV(baseGenFleet,'genFleetWithoutFilledInEmsRates.csv')
    # fill in missing values w/ avg of similar fuel & plant type
    fillMissingEmissionsRates(baseGenFleet, emsHeadersToAdd)


def fillMissingEmissionsRates(baseGenFleet, emsHeadersToAdd):
    """Fills missing generator emission rates with average for gens that have the same fuel and plant type.

    This function does not return anything.

    :param baseGenFleet: (2d list) generator fleet (see :func:`.importNEEDSFleet`)
    :param emsHeadersToAdd: (1d list) emissions headers to add
    """
    # Get headers and columns
    headersToColsMapBase = mapHeadersToCols(baseGenFleet)
    plantTypeCol = headersToColsMapBase['PlantType']
    fuelTypeCol = headersToColsMapBase['Modeled Fuels']
    noxCol = headersToColsMapBase[emsHeadersToAdd[0]]
    so2Col = headersToColsMapBase[emsHeadersToAdd[1]]
    co2Col = headersToColsMapBase[emsHeadersToAdd[2]]
    # Find and fill missing emissions rates values
    for idx in range(1, len(baseGenFleet)):
        if baseGenFleet[idx][noxCol] == 'NA':
            (plantType, fuelType) = (baseGenFleet[idx][plantTypeCol],
                                     baseGenFleet[idx][fuelTypeCol])
            [nox, so2, co2] = getEmsRatesOfMatchingFuelAndPlantType(baseGenFleet, plantType,
                                                                    fuelType, emsHeadersToAdd)
            [avgnox, avgso2, avgco2] = [avgListVals(nox), avgListVals(so2), avgListVals(co2)]
            baseGenFleet[idx][noxCol] = avgnox
            baseGenFleet[idx][so2Col] = avgso2
            baseGenFleet[idx][co2Col] = avgco2


def getEmsRatesOfMatchingFuelAndPlantType(baseGenFleet, plantType, fuelType, emsHeadersToAdd):
    """ Gets emission rates of generators w/ given plant & fuel type

    :param baseGenFleet: (2d list) generator fleet (see :func:`.importNEEDSFleet`)
    :param plantType: (string) name of plant type
    :param fuelType: (string) name of fuel type
    :param emsHeadersToAdd: (1d list) emissions headers to add
    :return: (1d lists) lists with NOx, SO2 and CO2 emissions rates
    """
    # Get headers
    headersToColsMapBase = mapHeadersToCols(baseGenFleet)
    noxCol = headersToColsMapBase[emsHeadersToAdd[0]]
    so2Col = headersToColsMapBase[emsHeadersToAdd[1]]
    co2Col = headersToColsMapBase[emsHeadersToAdd[2]]
    # Get cols w/ matching fuel & plant type
    matchingRowIdxs = getMatchingRowsFuelAndPlantType(baseGenFleet, plantType, fuelType,
                                                      noxCol)
    # If can't find on fuel & plant type, try just fuel type
    if matchingRowIdxs == []:
        matchingRowIdxs = getMatchingRowsFuelType(baseGenFleet, fuelType, noxCol)
    # If still can't get emissions rate, then ues other plant & fuel type:
    # LFG - NGCT, MSW - biomass, gas & oil O/G Steam - gas O/G Steam, Non-fossil waste -
    if matchingRowIdxs == [] and fuelType == 'Landfill Gas':
        matchingRowIdxs = getMatchingRowsFuelAndPlantType(baseGenFleet, 'Combustion Turbine',
                                                          'Natural Gas', noxCol)
    elif matchingRowIdxs == [] and fuelType == 'MSW':
        matchingRowIdxs = getMatchingRowsFuelAndPlantType(baseGenFleet, 'Biomass',
                                                          'Biomass', noxCol)
    elif matchingRowIdxs == [] and fuelType == 'Natural Gas& Distillate Fuel Oil& Residual Fuel Oil':
        matchingRowIdxs = getMatchingRowsFuelAndPlantType(baseGenFleet, 'O/G Steam',
                                                          'Natural Gas', noxCol)
    elif matchingRowIdxs == [] and fuelType == 'Non-Fossil Waste':
        matchingRowIdxs = getMatchingRowsFuelAndPlantType(baseGenFleet, 'Biomass',
                                                          'Biomass', noxCol)
    # Get emissions rates of matching rows
    [nox, so2, co2] = [[], [], []]
    for rowIdx in matchingRowIdxs:
        row = baseGenFleet[rowIdx]
        nox.append(row[noxCol])
        so2.append(row[so2Col])
        co2.append(row[co2Col])
    return [nox, so2, co2]


def getMatchingRowsFuelAndPlantType(baseGenFleet, plantType, fuelType, noxCol):
    """Gets row indexes in generator fleet of generators that match given plant & fuel type, filtering out units w/ no emissions rate data.

    :param baseGenFleet: (2d list) generator fleet (see :func:`.importNEEDSFleet`)
    :param plantType: (string) name of plant type
    :param fuelType: (string) name of fuel type
    :param noxCol: (int) index of column with NOx emission rate
    :return: (1d list) indices of row with matching plant & fuel type
    """
    headersToColsMapBase = mapHeadersToCols(baseGenFleet)
    plantTypeCol = headersToColsMapBase['PlantType']
    fuelTypeCol = headersToColsMapBase['Modeled Fuels']
    matchingRowIdxs = []
    for idx in range(len(baseGenFleet)):
        row = baseGenFleet[idx]
        if row[plantTypeCol] == plantType and row[fuelTypeCol] == fuelType:
            if row[noxCol] != 'NA':  # make sure has data!
                matchingRowIdxs.append(idx)
    return matchingRowIdxs


def getMatchingRowsFuelType(baseGenFleet, fuelType, noxCol):
    """Gets row indexes in generator fleet of gens with same fuel type, filtering out units with no emissions rate data.


    :param baseGenFleet: (2d list) generator fleet (see :func:`.importNEEDSFleet`)
    :param fuelType: (string) name of fuel type
    :param noxCol: (int) index of column with NOx emission rate
    :return: (1d list) indices of row with matching fuel type
    """
    headersToColsMapBase = mapHeadersToCols(baseGenFleet)
    fuelTypeCol = headersToColsMapBase['Modeled Fuels']
    matchingRowIdxs = []
    for idx in range(len(baseGenFleet)):
        row = baseGenFleet[idx]
        if row[fuelTypeCol] == fuelType:
            if row[noxCol] != 'NA':  # make sure has data!
                matchingRowIdxs.append(idx)
    return matchingRowIdxs


def addEmissionsRatesValues(baseFleet, egridBoiler, egridPlant):
    """Add eGRID emissions rates values to fleet

    This function adds eGRID emissions rates values to fleet, either using boiler specific data for coal & o/g steam
    units or plant level average data. Adds emission rate in order of nox, so2, and co2, as set by ems headers in
    addEmissionsRates. This function does not return anything.

    :param baseFleet: (2d list) generator fleet (see :func:`.importNEEDSFleet`)
    :param egridBoiler: (2d list) eGRID boiler data
    :param egridPlant: (2d list) eGRID plant data
    """
    headersToColsMapBase = mapHeadersToCols(baseFleet)
    headersToColsMapEgridBlr = mapHeadersToCols(egridBoiler)
    headersToColsMapEgridPlnt = mapHeadersToCols(egridPlant)
    basePlantTypeCol = headersToColsMapBase['PlantType']
    noEmissionPlantTypes = ['hydro', 'solar pv', 'wind', 'geothermal',
                            'solar thermal', 'pumped storage', 'nuclear']
    for idx in range(1, len(baseFleet)):
        plantType = baseFleet[idx][basePlantTypeCol].lower()
        if plantType == 'coal steam':
            [nox, so2, co2] = getBlrEmRates(baseFleet, idx, egridBoiler)
        elif plantType == 'o/g steam':
            [nox, so2, co2] = getBlrEmRates(baseFleet, idx, egridBoiler)
            if nox == 'NA':  # just test on nox, but all would be na
                [nox, so2, co2] = getPlantEmRates(baseFleet, idx, egridPlant)
        elif plantType in noEmissionPlantTypes:
            [nox, so2, co2] = [0, 0, 0]
        else:
            [nox, so2, co2] = getPlantEmRates(baseFleet, idx, egridPlant)
        # Some plants have no emissions info, so end up w/ zero emission values -
        # fill in 'NA' if so.
        if [nox, so2, co2] == [0, 0, 0] and plantType not in noEmissionPlantTypes:
            [nox, so2, co2] = ['NA', 'NA', 'NA']
        baseFleet[idx].extend([nox, so2, co2])


def getBlrEmRates(baseFleet, idx, egridBoiler):
    """Get boiler emission rates

    Look for boiler-level match of given gen in gen fleet to eGRID data, and return emissions rates if find match.

    :param baseFleet: (2d list) generator fleet (see :func:`.importNEEDSFleet`)
    :param idx: (int) index of row in `baseFleet`
    :param egridBoiler: (2d list) eGRID boiler data
    :return: (1d list) boiler-level nox, so2 & co2 ems rates
    """
    # Setup necessary data
    headersToColsMapBase = mapHeadersToCols(baseFleet)
    headersToColsMapEgridBlr = mapHeadersToCols(egridBoiler)
    (baseOrisCol, baseUnitCol) = (headersToColsMapBase["ORIS Plant Code"],
                                  headersToColsMapBase["Unit ID"])
    (egridOrisCol, egridBlrCol) = (headersToColsMapEgridBlr["DOE/EIA ORIS plant or facility code"],
                                   headersToColsMapEgridBlr["Boiler ID"])
    (egridBlrORISIDs, egridBlrIDs) = (colTo1dList(egridBoiler, egridOrisCol),
                                      colTo1dList(egridBoiler, egridBlrCol))
    # eGrid ORIS IDs are given w/ .0 @ end (e.g., 5834.0). So convert to int and back to str.
    removeTrailingDecimalFromEgridORIS(egridBlrORISIDs)
    # Do mapping
    (baseOrisID, baseUnitID) = (baseFleet[idx][baseOrisCol], baseFleet[idx][baseUnitCol])
    try:
        egridBlrRow = search2Lists(egridBlrORISIDs, egridBlrIDs, baseOrisID, baseUnitID)
        [nox, so2, co2] = calculateEmissionsRatesBlr(egridBoiler, egridBlrRow)
    except:
        # print('No matching boiler for: ORIS' + str(baseOrisID) + ' Blr' + str(baseUnitID))
        [nox, so2, co2] = ['NA', 'NA', 'NA']
    return [nox, so2, co2]


def getPlantEmRates(baseFleet, idx, egridPlant):
    """Get plant emission rates

    Looks for plant-level match of given unit in gen fleet to eGRID plant data, and returns plant-level ems rate of
    matching plant if found.

    :param baseFleet: (2d list) generator fleet (see :func:`.importNEEDSFleet`)
    :param idx: (int) index of row in `baseFleet`
    :param egridPlant: (2d list) eGRID plant data
    :return: (1d list) plant-level nox, so2 & co2 ems rate
    """
    # Setup necessary data
    headersToColsMapBase = mapHeadersToCols(baseFleet)
    headersToColsMapEgridPlnt = mapHeadersToCols(egridPlant)
    baseOrisCol = headersToColsMapBase["ORIS Plant Code"]
    egridOrisCol = headersToColsMapEgridPlnt["DOE/EIA ORIS plant or facility code"]
    egridORISIDs = colTo1dList(egridPlant, egridOrisCol)
    # eGrid ORIS IDs are given w/ .0 @ end (e.g., 5834.0). So convert to int and back to str.
    # removeTrailingDecimalFromEgridORIS(egridORISIDs)
    # Do mapping
    baseOrisID = baseFleet[idx][baseOrisCol]
    try:
        egridPlantRow = egridORISIDs.index(baseOrisID)
        [nox, so2, co2] = calculateEmissionsRatesPlnt(egridPlant, egridPlantRow)
    except:
        # print('No matching plant for: ORIS' + str(baseOrisID))
        [nox, so2, co2] = ['NA', 'NA', 'NA']
    return [nox, so2, co2]


def calculateEmissionsRatesBlr(egridBoiler, egridBoilerRow):
    """Calculate boiler-level emissions rates.

    :param egridBoiler: (2d list) eGRID boiler data
    :param egridBoilerRow: (int) row in boiler data
    :return: (1d list) boiler-level emission rates [lb/mmbtu]
    """
    scaleTonsToLbs = 2000
    # Define headers
    htInputHeader = 'Boiler unadjusted annual best heat input (MMBtu)'
    noxHeader = 'Boiler unadjusted annual best NOx emissions (tons)'
    so2Header = 'Boiler unadjusted annual best SO2 emissions (tons)'
    co2Header = 'Boiler unadjusted annual best CO2 emissions (tons)'
    # Calculate values
    headersToColsMap = mapHeadersToCols(egridBoiler)
    (htinputCol, noxCol, so2Col, co2Col) = (headersToColsMap[htInputHeader],
                                            headersToColsMap[noxHeader],
                                            headersToColsMap[so2Header],
                                            headersToColsMap[co2Header])
    blrData = egridBoiler[egridBoilerRow]
    (htInput, noxEms, so2Ems, co2Ems) = (blrData[htinputCol], blrData[noxCol],
                                         blrData[so2Col], blrData[co2Col])
    # Str nums have commas in them - use helper function to turn into numbers
    (htInput, noxEms, so2Ems, co2Ems) = (toNum(htInput), toNum(noxEms), toNum(so2Ems),
                                         toNum(co2Ems))
    (noxEmsRate, so2EmsRate, co2EmsRate) = (noxEms / htInput * scaleTonsToLbs,
                                            so2Ems / htInput * scaleTonsToLbs,
                                            co2Ems / htInput * scaleTonsToLbs)
    return [noxEmsRate, so2EmsRate, co2EmsRate]


def calculateEmissionsRatesPlnt(egridPlant, egridPlantRow):
    """Calculate plant-level emissions rates.

    :param egridPlant: (2d list) eGRID plant data
    :param egridPlantRow: (int) row in plant data
    :return: (1d list) plant-level nox, so2 and co2 emission rates [lb/mmbtu]
    """
    # Define headers
    noxEmsRateHeader = 'Plant annual NOx input emission rate (lb/MMBtu)'
    so2EmsRateHeader = 'Plant annual SO2 input emission rate (lb/MMBtu)'
    co2EmsRateHeader = 'Plant annual CO2 input emission rate (lb/MMBtu)'
    # Get values
    headersToColsMap = mapHeadersToCols(egridPlant)
    (noxCol, so2Col, co2Col) = (headersToColsMap[noxEmsRateHeader],
                                headersToColsMap[so2EmsRateHeader],
                                headersToColsMap[co2EmsRateHeader])
    plantData = egridPlant[egridPlantRow]
    (noxEmsRate, so2EmsRate, co2EmsRate) = [plantData[noxCol], plantData[so2Col], plantData[co2Col]]
    # Ems rate nums have commas - use helper func to turn into numbers
    (noxEmsRate, so2EmsRate, co2EmsRate) = (toNum(noxEmsRate),
                                            toNum(so2EmsRate),
                                            toNum(co2EmsRate))
    return [noxEmsRate, so2EmsRate, co2EmsRate]

#
################################################################################
#
################################################################################
#


def addLatLong(baseGenFleet, statesForAnalysis, dataRoot):
    """Add lat/long data to fleet

    This function gets lat/long data from egrid and adds columsn with this data to the generator fleet. This function
    does not return anything

    :param baseGenFleet: (2d list) generator fleet (see :func:`.importNEEDSFleet`)
    :param statesForAnalysis: (1-d list) list with states included in the analysis
    :param dataRoot: (string) path of directory wit input data
    """
    (egridBoiler, egridPlant) = importeGridData(statesForAnalysis, dataRoot)
    latLongHeadersToAdd = ["Latitude", "Longitude"]
    addHeaders(baseGenFleet, latLongHeadersToAdd)
    addLatLongValues(baseGenFleet, egridPlant)


def addLatLongValues(baseFleet, egridPlant):
    """ Worker function to add lat/long data to fleet

    This function adds lat/long values to base fleet using eGRID plant data. This function does not return anything

    :param baseFleet: (2d list) generator fleet (see :func:`.importNEEDSFleet`)
    :param egridPlant: (2d list) eGRID plant data
    """
    headersToColsMapBase = mapHeadersToCols(baseFleet)
    headersToColsMapEgrid = mapHeadersToCols(egridPlant)
    baseOrisCol = headersToColsMapBase['ORIS Plant Code']
    egridOrisCol = headersToColsMapEgrid['DOE/EIA ORIS plant or facility code']
    (egridLatCol, egridLongCol) = (headersToColsMapEgrid['Plant latitude'],
                                   headersToColsMapEgrid['Plant longitude'])
    egridORISIDs = colTo1dList(egridPlant, egridOrisCol)
    for idx in range(1, len(baseFleet)):
        genRow = baseFleet[idx]
        genORIS = genRow[baseOrisCol]
        if genORIS in egridORISIDs:
            egridRow = egridORISIDs.index(genORIS)
            (genLat, genLong) = (egridPlant[egridRow][egridLatCol],
                                 egridPlant[egridRow][egridLongCol])
        else:
            (genLat, genLong) = ('NA', 'NA')
        baseFleet[idx].extend([genLat, genLong])


################################################################################


def performFleetCompression(genFleet, ipmZones, plantTypesCurtailed):
    """Compress fleet by combining small units

    This function simplifies the generator fleet by combining small units (< 200 MW) of same type into a single plant.
    Units are combined according to plant type, subregion (ipm zone), age, among other factors. Thermal units that will
    have capacity deratings will not be combined (since capacity deratings are location-dependent).

    :param genFleet: (2d list) generator fleet (see :func:`.importNEEDSFleet`)
    :param ipmZones: (1-d list) list with names of ipm zones
    :param plantTypesCurtailed: (1-d list) list with names of plant types that are included in capacity derating simulations
    :return: (2-d list) updated compressed generator fleet with small units combined
    """
    fuelAndPlantTypeToCompress = [('Landfill Gas', 'Landfill Gas'),
                                  ('Distillate Fuel Oil', 'Combustion Turbine'), ('MSW', 'Municipal Solid Waste'),
                                  ('Natural Gas', 'Combustion Turbine'), ('Biomass', 'Biomass'),
                                  ('Natural Gas& Distillate Fuel Oil', 'O/G Steam'), ('Natural Gas', 'O/G Steam'),
                                  ('Natural Gas', 'Combined Cycle')]
    for (fuel, plant) in fuelAndPlantTypeToCompress:
        compressFuelAndPlantType(genFleet, fuel, plant, ipmZones, plantTypesCurtailed)

    return genFleet


def compressFuelAndPlantType(genFleet, fuel, plant, ipmZones, plantTypesCurtailed):
    """Compress fleet by combining small units of given plant type and fuel type

    This function modifies the 2-d list `genFleet` passed as parameter. It does not return anything

    :param genFleet: (2d list) generator fleet (see :func:`.importNEEDSFleet`)
    :param fuel: (string) name of fuel type
    :param plant: (string) name of plant type
    :param ipmZones: (1-d list) list with names of ipm zones
    :param plantTypesCurtailed: (1-d list) list with names of plant types that are included in capacity derating simulations
    """
    maxSizeToCombine = 200
    head = genFleet[0]
    (plantCol, fuelCol, capacCol) = (head.index('PlantType'), head.index('Modeled Fuels'),
                                     head.index('Capacity (MW)'))
    (coolTechCol, coolSourceCol) = (head.index('Cooling Tech'), head.index('Cooling Source'))
    idxsToRemoveAndCombine = []
    startFleetLength = len(genFleet)

    for idx in range(1, startFleetLength):
        rowFuel = isolateFirstFuelType(genFleet[idx][fuelCol])
        # Only combine plants that meet criteria AND don't have a cooling tech or source listed
        # if the plant will be curtailed.
        if rowFuel == fuel and genFleet[idx][plantCol] == plant and float(genFleet[idx][capacCol]) < maxSizeToCombine:
            if plant in plantTypesCurtailed:
                if genFleet[idx][coolTechCol] == 'NoMatch' and genFleet[idx][coolSourceCol] == 'NoMatch':
                    idxsToRemoveAndCombine.append(idx)
            else:
                idxsToRemoveAndCombine.append(idx)

    combineGenerators(genFleet, idxsToRemoveAndCombine, fuel, plant, startFleetLength, ipmZones)

    for idx in reversed(idxsToRemoveAndCombine):
        genFleet.pop(idx)


def combineGenerators(genFleet, idxsToRemoveAndCombine, fuel, plant, startFleetLength, ipmZones):
    """Combine generators based on when they came online

    Small generators of same fuel/plant-type that came online in the same decade will be combined. This function
    modifies the 2-d list `genFleet` passed as paraemter. It does not return anything

    :param genFleet: (2d list) generator fleet (see :func:`.importNEEDSFleet`)
    :param idxsToRemoveAndCombine: (1-d list) list with indexes of rows of units that are going to be combined (and then removed from fleet)
    :param fuel: (string) name of fuel type
    :param plant: (string) name of plant type
    :param startFleetLength: (int) length of original uncompressed fleet
    :param ipmZones: (1-d list) list with names of ipm zones
    """
    onlineYearCol = genFleet[0].index('On Line Year')
    zoneCol = genFleet[0].index('Region Name')
    onlineYears = [int(genFleet[idx][onlineYearCol]) for idx in idxsToRemoveAndCombine]
    (firstYr, lastYr, stepYr) = (1975, 2026, 10)
    yearIntervals = [yr for yr in range(firstYr, lastYr, stepYr)]

    for zone in ipmZones:
        for endingYear in yearIntervals:
            if endingYear == firstYr:
                beginningYear = 0
            else:
                beginningYear = endingYear - stepYr
            idxsInInterval = [idx for idx in idxsToRemoveAndCombine
                              if (beginningYear < int(genFleet[idx][onlineYearCol]) <= endingYear)]
            idxsInIntervalAndZone = [idx for idx in idxsInInterval if genFleet[idx][zoneCol] == zone]

            if len(idxsInInterval) > 0:
                combineGeneratorsInDecade(genFleet, idxsInIntervalAndZone, endingYear - stepYr // 2, fuel, plant,
                                          startFleetLength, zone)


def combineGeneratorsInDecade(genFleet, idxsInIntervalAndZone, medianYearInInterval, fuel, plant,
                              startFleetLength, zone):
    """Combine generators in same decade

    This function does the actual aggregation of generators in the same decade. It combines generators up to 500 MW
    of combined size. This function modifies the 2-d list `genFleet` passed as parameter. It does not return anything

    :param genFleet: (2d list) generator fleet (see :func:`.importNEEDSFleet`)
    :param idxsInIntervalAndZone: (1-d list) list with units in the time interval and ipm zone that can be combined
    :param medianYearInInterval: (int) median online year of plants being combined
    :param fuel: (string) name of fuel type
    :param plant: (string) name of plant type
    :param startFleetLength: (int) length of original uncompressed fleet
    :param zone: (string) name of ipm zone
    """
    maxCombinedSize = 500
    (runningCombinedSize, idxsToCombine) = (0, [])
    capacCol = genFleet[0].index('Capacity (MW)')
    for idx in idxsInIntervalAndZone:
        if runningCombinedSize + float(genFleet[idx][capacCol]) > maxCombinedSize:
            addCombinedIdxsToFleet(genFleet, idxsToCombine, runningCombinedSize, fuel, plant,
                                   medianYearInInterval, zone)
            runningCombinedSize = float(genFleet[idx][capacCol])
            idxsToCombine = [idx]
        else:
            runningCombinedSize += float(genFleet[idx][capacCol])
            idxsToCombine.append(idx)

    if len(idxsToCombine) > 0:  # combine remaining units
        addCombinedIdxsToFleet(genFleet, idxsToCombine, runningCombinedSize, fuel, plant,
                               medianYearInInterval, zone)


def addCombinedIdxsToFleet(genFleet, idxsToCombine, combinedCapac, fuel, plant, medianYearInInterval, zone):
    """Combine generators in given indexes

    This function does the actual aggregation of generators at given rows in the fleet's 2-d list. It also appends
    the row with the new generator to the fleet. This function modifies the 2-d list `genFleet` passed as parameter.
    It does not return anything

    :param genFleet: (2d list) generator fleet (see :func:`.importNEEDSFleet`)
    :param idxsToCombine: (1-d list) list with indexes of rows of generators to be combined
    :param combinedCapac: (float) total capacity of combined power plant
    :param fuel: (string) name of fuel type
    :param plant: (string) name of plant type
    :param medianYearInInterval: (int) median online year of plants being combined
    :param zone: (string) name of ipm zone
    """
    headers = genFleet[0]
    newRow = [''] * len(headers)
    rowsToCombine = [genFleet[idx] for idx in idxsToCombine]
    capacCol = headers.index('Capacity (MW)')
    capacWts = [float(row[capacCol]) / combinedCapac for row in rowsToCombine]
    parametersToCombine = ['NOxEmRate(lb/MMBtu)', 'SO2EmRate(lb/MMBtu)',
                           'CO2EmRate(lb/MMBtu)', 'Heat Rate (Btu/kWh)']
    for param in parametersToCombine:
        colNum = headers.index(param)
        paramVals = [float(row[colNum]) for row in rowsToCombine]
        newRow[colNum] = sum(list(map(operator.mul, capacWts, paramVals)))
    newRow[capacCol] = combinedCapac
    addStateZoneOrisFuelOnlineYearAndPlantType(genFleet, newRow, fuel, plant, zone, medianYearInInterval)

    genFleet.append(newRow)


def addStateZoneOrisFuelOnlineYearAndPlantType(genFleet, newRow, fuel, plant, zone, *onlineYear):
    """Add additional info to new combined generator

    This function adds State, Zone, Oris, Fuel, OnlineYear, and PlantType to the new generator resulting from the
    combination of small units. This function modifies the 1-d list `newrow` passed as parameter. It does not return
    anything

    :param genFleet: (2d list) generator fleet (see :func:`.importNEEDSFleet`)
    :param newRow: (1-d list) list with data of new combined generator that will be added to the fleet
    :param fuel: (string) name of fuel type
    :param plant: (string) name of plant type
    :param zone: (string) name of ipm zone
    :param onlineYear: (int) online year of new combined generator that will be added to the fleet
    """
    (stateCol, orisCol) = (genFleet[0].index('State Name'), genFleet[0].index('ORIS Plant Code'))
    (unitCol, fuelCol) = (genFleet[0].index('Unit ID'), genFleet[0].index('Modeled Fuels'))
    (onlineYearCol, ipmRetireCol) = (genFleet[0].index('On Line Year'), genFleet[0].index('Retirement Year'))
    zoneCol = genFleet[0].index('Region Name')
    plantCol = genFleet[0].index('PlantType')
    maxOris = max([int(row[orisCol]) for row in genFleet[1:]])
    (newRow[orisCol], newRow[stateCol], newRow[unitCol]) = (maxOris + 1, '', '1')
    (newRow[fuelCol], newRow[plantCol], newRow[zoneCol]) = (fuel, plant, zone)
    newRow[ipmRetireCol] = 9999

    if len(onlineYear) > 0: newRow[onlineYearCol] = onlineYear[0]


################################################################################
#
################################################################################


def addVOMandFOM(baseGenFleet, vomAndFomData):
    """Add columns with variable and fixed O&M data to generator fleet

    This function modifies the 2-d list `baseGenFleet` passed as parameter. It does not return anything

    :param baseGenFleet: (2d list) generator fleet (see :func:`.importNEEDSFleet`)
    :param vomAndFomData: (2-d list) 2-d list with variable of fixed O&M data for each plant type
    """
    vomAndFomHeader = ['VOM($/MWh)', 'FOM($/MW/yr)']
    addHeaders(baseGenFleet, vomAndFomHeader)
    addVomAndFomValues(baseGenFleet, vomAndFomData)


def addVomAndFomValues(baseGenFleet, vomAndFomData):
    """Add variable and fixed O&M data to generator fleet

    This function modifies the 2-d list `baseGenFleet` passed as parameter. It does not return anything

    :param baseGenFleet: (2d list) generator fleet (see :func:`.importNEEDSFleet`)
    :param vomAndFomData: (2-d list) 2-d list with variable of fixed O&M data for each plant type
    """
    plantTypeColFleet = baseGenFleet[0].index('PlantType')
    (vomColFleet, fomColFleet) = (baseGenFleet[0].index('VOM($/MWh)'), baseGenFleet[0].index('FOM($/MW/yr)'))
    plantTypeColData = vomAndFomData[0].index('PlantType')
    plantTypesData = [row[plantTypeColData] for row in vomAndFomData]
    (vomColData, fomColData) = (vomAndFomData[0].index('VOM(2012$/MWh)'), vomAndFomData[0].index('FOM(2012$/MW/yr)'))
    for row in baseGenFleet[1:]:
        dataRow = plantTypesData.index(row[plantTypeColFleet])
        vomValue = convertCostToTgtYr('vom', float(vomAndFomData[dataRow][vomColData]))
        fomValue = convertCostToTgtYr('fom', float(vomAndFomData[dataRow][fomColData]))
        if vomColFleet >= len(row):
            row.extend([vomValue, fomValue])  # haven't added VOM or FOM values yet
        else:
            (row[vomColFleet], row[fomColFleet]) = (
            vomValue, fomValue)  # filling in values in row that already has blank space for VOM & FOM


def importVomAndFomData(dataRoot):
    """Import variable and fixed O&M data from input files

    Variable and fixed O&M data are located in file `VOMandFOMValuesExistingPlants4Aug2016.csv` inside the directory
    `[dataRoot]/NewPlantData/`

    :param dataRoot: (string) path of directory wit input data
    :return: (2-d list) 2-d list with variable of fixed O&M data for each plant type
    """
    dirName = os.path.join(dataRoot, 'NewPlantData')

    fileName = 'VOMandFOMValuesExistingPlants4Aug2016.csv'
    fullFileName = os.path.join(dirName, fileName)
    return readCSVto2dList(fullFileName)


################################################################################

def addUnitCommitmentParameters(baseGenFleet, phorumData):
    """Add unit commitment parameters to generator fleet

    unit commitment parameters are based on fuel and plant type (data from PHORUM). Columns with UC are added to
    the `baseGenFleet` 2-d list.

    :param baseGenFleet: (2d list) generator fleet (see :func:`.importNEEDSFleet`)
    :param phorumData: (2d list) table with phorum parameters (see :func:`.importPhorumData`)
    """
    ucHeaders = ['MinDownTime(hrs)', 'RampRate(MW/hr)', 'MinLoad(MW)', 'StartCost($)']
    addHeaders(baseGenFleet, ucHeaders)
    addUCValues(baseGenFleet, ucHeaders, phorumData)


def addUCValues(baseGenFleet, ucHeaders, phorumData):
    """Worker function to ddd unit commitment parameters to generator fleet

    :param baseGenFleet: (2d list) generator fleet (see :func:`.importNEEDSFleet`)
    :param ucHeaders:  (1d list) list with names of UC parameters that will be added to header of generator fleet
    :param phorumData: (2d list) table with phorum parameters (see :func:`.importPhorumData`)
    """
    capacCol = baseGenFleet[0].index('Capacity (MW)')
    fuelCol = baseGenFleet[0].index('Modeled Fuels')
    plantTypeCol = baseGenFleet[0].index('PlantType')
    for ucHeader in ucHeaders:
        currCol = baseGenFleet[0].index(ucHeader)
        ucHeaderToPhorumParamName = mapHeadersToPhorumParamNames()
        phorumParamName = ucHeaderToPhorumParamName[ucHeader]
        for row in baseGenFleet[1:]:
            (fuel, plantType, size) = (row[fuelCol], row[plantTypeCol], float(row[capacCol]))

            # some plants have multiple modeled fuels divided by &
            fuel = isolateFirstFuelType(fuel)

            phorumValue = getMatchingPhorumValue(phorumData, fuel, plantType, size, phorumParamName)

            if ucHeader == 'MinDownTime(hrs)':
                valToAdd = phorumValue
            else:
                valToAdd = phorumValue * size
                if ucHeader == 'StartCost($)': valToAdd = convertCostToTgtYr('startup', valToAdd)
            if currCol >= len(row):
                row.append(valToAdd)  # for when first adding values
            else:
                row[currCol] = valToAdd  # for filling in values after have already added


def getMatchingPhorumValue(phorumData, fuel, plantType, size, paramName):
    """Get value of UC parameter from Phorum database

    :param phorumData: (2d list) table with phorum parameters (see :func:`.importPhorumData`)
    :param fuel: (string) fuel name
    :param plantType: (string) plant type name
    :param size: (float) installed capacity of plant
    :param paramName: (string) name of UC parameter
    :return: (float) value of UC parameter for this plant type & fuel in the phorum database
    """
    if plantType == 'Fuel Cell': plantType = 'Combustion Turbine'
    fuel = mapFleetFuelToPhorumFuels(fuel)
    phorumPropertyNameCol = phorumData[0].index('PropertyName')
    phorumFuelCol = phorumData[0].index('Fuel')
    phorumPlantTypeCol = phorumData[0].index('PlantType')
    phorumLowerSizeCol = phorumData[0].index('LowerPlantSizeLimit')
    phorumUpperSizeCol = phorumData[0].index('UpperPlantSizeLimit')
    phorumValueCol = phorumData[0].index('PropertyValue')
    phorumProperties = [row[phorumPropertyNameCol] for row in phorumData[1:]]
    phorumFuels = [row[phorumFuelCol] for row in phorumData[1:]]
    phorumPlantTypes = [row[phorumPlantTypeCol] for row in phorumData[1:]]
    phorumLowerSizes = [int(row[phorumLowerSizeCol]) for row in phorumData[1:]]
    phorumUpperSizes = [int(row[phorumUpperSizeCol]) for row in phorumData[1:]]
    phorumValues = [float(row[phorumValueCol]) for row in phorumData[1:]]
    for idx in range(len(phorumProperties)):
        if (phorumProperties[idx] == paramName and phorumFuels[idx] == fuel and
                (phorumPlantTypes[idx] == plantType or phorumPlantTypes[idx] == 'All') and
                (phorumLowerSizes[idx] <= size and phorumUpperSizes[idx] > size)):
            return float(phorumValues[idx])


def mapFleetFuelToPhorumFuels(fleetFuel):
    """Map name of fuel in fleet to Phorum data

    :param fleetFuel: (string) name if fuel in fleet
    :return: (string) name of fuel in Phorum
    """
    fleetFuelToPhorumFuelMap = {'Bituminous': 'Coal', 'Petroleum Coke': 'Pet. Coke',
                                'Subbituminous': 'Coal', 'Lignite': 'Coal', 'Natural Gas': 'NaturalGas',
                                'Distillate Fuel Oil': 'Oil', 'Hydro': 'Hydro', 'Landfill Gas': 'LF Gas',
                                'Biomass': 'Biomass', 'Solar': 'Solar', 'Non-Fossil Waste': 'Non-Fossil',
                                'MSW': 'MSW', 'Pumped Storage': 'Hydro', 'Residual Fuel Oil': 'Oil',
                                'Wind': 'Wind', 'Nuclear Fuel': 'Nuclear', 'Coal': 'Coal'}
    return fleetFuelToPhorumFuelMap[fleetFuel]


def mapHeadersToPhorumParamNames():
    """Get map name of UC parameter in fleet to Phorum

    This function returns a dictionary mapping the names of UC parameters used in the fleet to the ones
    used in Phorum.

    :return: dictionary {'name in fleet': 'name in phorum'}
    """
    return {'MinDownTime(hrs)': 'Min Down Time', 'RampRate(MW/hr)': 'Ramp Rate',
            'MinLoad(MW)': 'Min Stable Level', 'StartCost($)': 'Start Cost'}


################################################################################

################################################################################
# ADD FUEL PRICES
def addFuelPrices(baseGenFleet, currYear, fuelPriceTimeSeries):
    """Add column with fuel prices to fleet

    This function adds a column of fuel prices to the 2-d list with fleet data

    :param baseGenFleet: (2-d list) generator fleet (see :func:`.importNEEDSFleet`)
    :param currYear: (int) current year
    :param fuelPriceTimeSeries: (1-d list) list with time series of fuel prices
    """
    fuelPriceHeader = ['FuelPrice($/MMBtu)']
    addHeaders(baseGenFleet, fuelPriceHeader)
    addFuelPriceValues(baseGenFleet, fuelPriceTimeSeries, currYear)


def addFuelPriceValues(baseGenFleet, fuelPriceTimeSeries, currYear):
    """Add fuel values to column of fuel prices

    :param baseGenFleet: (2-d list) generator fleet (see :func:`.importNEEDSFleet`)
    :param fuelPriceTimeSeries: (1-d list) list with time series of fuel prices
    :param currYear: (int) current year
    """
    fuelCol = baseGenFleet[0].index('Modeled Fuels')
    for row in baseGenFleet[1:]: row.append(getFuelPrice(row, fuelCol, fuelPriceTimeSeries, currYear))


def getFuelPrice(fleetRow, fleetFuelCol, fuelPriceTimeSeries, currYear):
    """Get price of fuel of a specific plant in fleet for current year (adjusted for inflation)

    :param fleetRow: (1-d list) row in generator fleet table with data for a specific plant
    :param fleetFuelCol: (int) index of column with fuel prices in gen fleet
    :param fuelPriceTimeSeries: (1-d list) list with time series of fuel prices
    :param currYear: (int) current year
    :return: value of price of  fuel adjusted to inflation
    """
    fuel = fleetRow[fleetFuelCol]
    fuel = isolateFirstFuelType(fuel)
    fuelPriceDollarUnadjusted = getFuelPriceForFuelType(fuel, fuelPriceTimeSeries, currYear)
    return convertCostToTgtYr('fuel', fuelPriceDollarUnadjusted)


def getFuelPriceForFuelType(fuel, fuelPriceTimeSeries, currYear):
    """Get price value of given fuel for current year (not adjusted for inflation)

    :param fuel: (string) name of fuel
    :param fuelPriceTimeSeries: (1-d list) list with time series of fuel prices
    :param currYear: (int) current year
    :return: value of price of fuel (not adjusted to inflation)
    """
    fuel = mapFleetFuelToPhorumFuels(fuel)
    fuelPriceFuelCol = fuelPriceTimeSeries[0].index('FuelPrices($/MMBtu)')
    fuelPriceFuels = [row[fuelPriceFuelCol] for row in fuelPriceTimeSeries]
    fuelPriceYears = [int(yr) for yr in fuelPriceTimeSeries[0][1:]]
    fuelPriceRow = fuelPriceFuels.index(fuel)
    fuelPricePrices = [float(price) for price in fuelPriceTimeSeries[fuelPriceRow][1:]]
    if currYear in fuelPriceYears:
        yearCol = fuelPriceYears.index(currYear)
    elif currYear > max(fuelPriceYears):
        yearCol = fuelPriceYears.index(max(fuelPriceYears))
    return fuelPricePrices[yearCol]


################################################################################

################################################################################

def importNEEDSFleet(dataRoot):
    """Import base generator fleet from NEEDS

    Reads NEEDS fleet data from a CSV file. This CSV file must be located in the `[dataRoot]/NEEDS` folder

    The CSV file `needs_nocommas.csv` is originated from the NEEDS excel file. Commas in the excel file must be
    converted to `&`

    :param dataRoot: (string) path to directory with input data
    :return: (2-d list) 2-d list with headers containing generator fleet data
    """
    dirName = os.path.join(dataRoot, 'NEEDS')

    fileName = 'needs_v515_nocommas.csv'
    fullFileName = os.path.join(dirName, fileName)
    return readCSVto2dList(fullFileName)


def importeGridData(statesForAnalysis, dataRoot):
    """Import eGRID boiler and plant level data, then isolate plants and boilers in state

    :param statesForAnalysis: (1d list) states for analysis
    :param dataRoot: (string) path of directory wit input data
    :return: (2d lists) eGRID boiler and plant data
    """
    dirName = os.path.join(dataRoot, 'eGRID2015')

    egridBoiler = importeGridBoilerData(dirName)
    egridPlant = importeGridPlantData(dirName)
    egridStateColName = 'Plant state abbreviation'
    statesForAnalysisAbbrev = getStateAbbrevs(statesForAnalysis)
    isolateGensInStates(egridBoiler, statesForAnalysisAbbrev, egridStateColName)
    isolateGensInStates(egridPlant, statesForAnalysisAbbrev, egridStateColName)
    return (egridBoiler, egridPlant)


def importeGridBoilerData(dirName):
    """Import eGRID boiler data and remove extra headers

    The CSV file `egrid_data_boiler.csv` is originated from the eGrid excel file. Commas in the excel file must be
    converted to `&`

    :param dirName: (string) directory with egrid data
    :return: (2d list) boiler data
    """
    fileName = 'egrid2012_data_boiler.csv'
    fullFileName = os.path.join(dirName, fileName)
    boilerData = readCSVto2dList(fullFileName)
    boilerDataSlim = elimExtraneousHeaderInfo(boilerData, 'eGRID2012 file boiler sequence number')
    return boilerDataSlim


def importeGridPlantData(dirName):
    """Import eGRID plant data and remove extra headers

    The CSV file `egrid_data_plant.csv` is originated from the eGrid excel file. Commas in the excel file must be
    converted to `&`

    :param dirName: (string) directory with egrid data
    :return: (2d list) plant data
    """
    fileName = 'egrid2012_data_plant.csv'
    fullFileName = os.path.join(dirName, fileName)
    plantData = readCSVto2dList(fullFileName)
    plantDataSlim = elimExtraneousHeaderInfo(plantData, 'eGRID2012 file plant sequence number')
    return plantDataSlim


def elimExtraneousHeaderInfo(egridFleet, valueInFirstValidRow):
    """Eliminates first several rows in egrid CSV that has no useful info

    :param egridFleet: (2d list) eGRID fleet
    :param valueInFirstValidRow: (string) value in col 0 in first row w/ valid data that want to save
    :return: (2d list) eGRID fleet
    """
    for idx in range(len(egridFleet)):
        if egridFleet[idx][0] == valueInFirstValidRow:
            egridFleetSlim = copy.deepcopy(egridFleet[idx:])
    return egridFleetSlim


def removeRetiredUnits(baseGenFleet, retirementYearScreen):
    """Removes retired units from fleet based on input year

    :param baseGenFleet: (2d list) gen fleet (see :func:`.importNEEDSFleet`)
    :param retirementYearScreen: (int) year below which retired units should be removed from fleet
    """
    colName = "Retirement Year"
    colNum = baseGenFleet[0].index(colName)
    rowsToRemove = []
    for rowIdx in range(1, len(baseGenFleet)):
        retireYear = baseGenFleet[rowIdx][colNum]
        if int(retireYear) < retirementYearScreen: rowsToRemove.append(rowIdx)
    if rowsToRemove != []: removeRows(baseGenFleet, rowsToRemove)


def isolateGensInStates(baseGenFleet, statesForAnalysis, colName):
    """Isolates fleet to generators in states of interest

    :param baseGenFleet: (2d list) gen fleet data (see :func:`.importNEEDSFleet`)
    :param statesForAnalysis: (1d list) states for analysis
    :param colName: (string) name of column name with state data
    :return: (2d list) updated gen fleet data
    """
    rowsToRemove = identifyRowsToRemove(baseGenFleet, statesForAnalysis, colName)
    removeRows(baseGenFleet, rowsToRemove)
    return baseGenFleet


def isolateGensInPowerSystem(baseGenFleet, ipmZones):
    """Isolates fleet to generators in in the ipm regions of interest

    :param baseGenFleet: (2d list) gen fleet (see :func:`.importNEEDSFleet`)
    :param ipmZones: (1d list) list with names of ipm zones included in analysis
    :return: (2d list) updated gen fleet
    """
    colName = "Region Name"
    rowsToRemove = identifyRowsToRemove(baseGenFleet, ipmZones, colName)
    removeRows(baseGenFleet, rowsToRemove)
    return baseGenFleet


def importPhorumData(dataRoot):
    """Import PHORUM data (VOM + UC parameters)

    Phorum data is in a CSV file `PHORUMUCParameters.csv`.

    :param dataRoot: (string) path to directory with input data
    :return: (2d list) phorum data
    """
    dirName = os.path.join(dataRoot, 'PHORUM')

    fileName = 'PHORUMUCParameters.csv'

    return readCSVto2dList(os.path.join(dirName, fileName))


################################################################################

################################################################################


def addRandomOpCostAdder(baseGenFleet, ocAdderMin, ocAdderMax):
    """Add random operation cost adder to fleet in new column

    Add to all fuel types. Use value of 0.05 - makes up ~0.03% on average of fleet. Max addition to op cost of gen in fleet is 0.19%.

    :param baseGenFleet: (2d list) gen fleet (see :func:`.importNEEDSFleet`)
    :param ocAdderMin: (float) min value
    :param ocAdderMax: (float) max value
    """
    randValHeader = 'RandOpCostAdder($/MWh)'
    addHeaders(baseGenFleet, [randValHeader])
    randValCol = baseGenFleet[0].index(randValHeader)

    #random.seed()

    for row in baseGenFleet[1:]: row.append(random.uniform(ocAdderMin, ocAdderMax))


################################################################################

################################################################################


def addRegResOfferAndElig(baseGenFleet, regupCostCoeffs):
    """Add regulated reserve offer costs and eligibility

    (Based on params in Denholm 2013, Val of energy sto for grid apps)

    :param baseGenFleet: (2d list) gen fleet (see :func:`.importNEEDSFleet`)
    :param regupCostCoeffs:
    """
    baseGenFleet[0].extend(['RegOfferCost($/MW)', 'RegOfferElig'])
    plantTypeCol = baseGenFleet[0].index('PlantType')
    for row in baseGenFleet[1:]:
        currPlantType = row[plantTypeCol]
        if currPlantType in regupCostCoeffs:
            regupCost, regElig = regupCostCoeffs[currPlantType], 1
        else:
            regupCost, regElig = 0, 0
        row.extend([regupCost, regElig])


################################################################################

################################################################################

def aggregatePlantTypeToORIS(genFleet, pt):
    """Aggregate hydro units to single plant according to ORIS code

    :param genFleet: (2d list) gen fleet (see :func:`.importNEEDSFleet`)
    :param pt: (string) name of plant type
    """
    orisCol, plantCol = genFleet[0].index('ORIS Plant Code'), genFleet[0].index('PlantType')
    capacCol = genFleet[0].index('Capacity (MW)')
    # Just sum capacities for each ORIS, keep track of which idxs to remove
    orisToCapac, idxToRemove = dict(), list()
    for idx in range(len(genFleet)):
        if genFleet[idx][plantCol] == pt:
            oris, capac = genFleet[idx][orisCol], float(genFleet[idx][capacCol])
            if oris in orisToCapac:
                orisToCapac[oris] += capac
                idxToRemove.append(idx)
            else:
                orisToCapac[oris] = capac
    # Remove rows
    if idxToRemove != []: removeRows(genFleet, idxToRemove)
    # Replace capacities of ORIS IDs
    for row in genFleet:
        if row[orisCol] in orisToCapac: row[capacCol] = orisToCapac[row[orisCol]]


################################################################################

################################################################################
# GENERAL UTILITY FUNCTIONS
def isolateFirstFuelType(fuel):
    """Helper function to split name of fuel

    Some plants have multiple modeled fuels divided by ``&``

    :param fuel: (string) original name of fuel
    :return: (string) simplified name of fuel
    """
    multiFuelDivider = '&'  # some plants have multiple modeled fuels divided by &

    fuel = fuel.split(multiFuelDivider)[0]

    return fuel


def getStateAbbrevs(statesForAnalysis):
    """Get abbreviations (which eGRID uses but NEEDS does not)

    :param statesForAnalysis: (1d list) states for analysis
    :return: (dict) map of states names to state abbreviations
    """
    stateAbbreviations = {'Virginia': 'VA', 'North Carolina': 'NC', 'South Carolina': 'SC',
                          'Georgia': 'GA', 'Mississippi': 'MS', 'Alabama': 'AL', 'Louisiana': 'LA',
                          'Missouri': 'MO', 'Arkansas': 'AR', 'Illinois': 'IL',
                          'Kentucky': 'KY', 'Tennessee': 'TN', 'Texas': 'TX'}
    statesForAnalysisAbbrev = []
    for state in statesForAnalysis:
        statesForAnalysisAbbrev.append(stateAbbreviations[state])
    return statesForAnalysisAbbrev


def identifyRowsToRemove(list2d, valuesToKeep, colName):
    """Returns a list of rows to remove for values in a given column that don't equal any value in valuesToKeep.

    :param list2d: any 2d list
    :param valuesToKeep: (1d list) values in specified column to keep
    :param colName: (string) col name
    :return: (1d list) row indices to remove
    """
    headersToColsMap = mapHeadersToCols(list2d)
    colNumber = headersToColsMap[colName]
    rowsToRemove = []
    for row in range(1, len(list2d)):
        if list2d[row][colNumber] not in valuesToKeep:
            rowsToRemove.append(row)
    return rowsToRemove


def removeRows(baseGenFleet, rowsToRemove):
    """ Remove rows of given indexes from 2d list

    :param baseGenFleet: (2d list) data (see :func:`.importNEEDSFleet`)
    :param rowsToRemove: (1d list) list with row indexes to remove
    """
    for row in reversed(rowsToRemove):
        baseGenFleet.pop(row)


def mapHeadersToCols(fleet):
    """Returns a dictionary mapping headers to column numbers

    :param fleet: (2d list) fleet data with header (see :func:`.importNEEDSFleet`)
    :return: (dict) map of header name to header
    """
    headers = fleet[0]
    headersToColsMap = dict()
    for colNum in range(len(headers)):
        header = headers[colNum]
        headersToColsMap[header] = colNum
    return headersToColsMap


def addHeaders(fleet, listOfHeaders):
    """

    :param fleet: (2d list) fleet data (see :func:`.importNEEDSFleet`)
    :param listOfHeaders: (1d list) headers to add to first row of data
    """
    for header in listOfHeaders:
        fleet[0].append(header)


def avgListVals(listOfVals):
    """Returns average of values in input 1d list

    :param listOfVals: (1-d list) list with numeric values
    :return: (float) average value
    """
    (total, count) = (0, 0)
    for val in listOfVals:
        total += float(val)
        count += 1
    return total / count


def removeTrailingDecimalFromEgridORIS(egridORISIDs):
    """Removes '.0' from end of ORIS IDs in eGRID

    :param egridORISIDs: (1-d list) list with ORIS IDs
    """
    for idx in range(1, len(egridORISIDs)):
        egridORISIDs[idx] = egridORISIDs[idx][:-2]


def toNum(s):
    """Converts a string w/ commas in it to a float

    :param s: (string) a number in string format
    :return: (float) numeric value
    """
    numSegments = s.split(',')
    result = ""
    for segment in numSegments:
        result += segment
    return float(result)



def search2Lists(list1, list2, data1, data2):
    """Return row indexes (or False) where list1=data1 and list2=data2

    :param list1: list
    :param list2: list
    :param data1: value
    :param data2: value
    :return: False if not match. If there is a match, index of match
    """
    if (data1 not in list1) or (data2 not in list2):
        return False

    for idx in range(len(list1)):
        if list1[idx] == data1 and list2[idx] == data2:
            return idx

    return False


def colTo1dList(data, colNum):
    """Convert specified column in 2d list to a 1-d list

    :param data: (2-d list) a generic 2-d list
    :param colNum:  (int) index of column
    :return: (1-d list) list with data from column
    """
    listWithColData = []
    for dataRow in data:
        listWithColData.append(dataRow[colNum])
    return listWithColData


################################################################################

################################################################################
# TEST FUNCTIONS
# Several utility functions used here are tested in other scripts. Other
# functions were tested manually using output test fleets.
def testAvgListVals():
    print('testing avgListVals')
    assert (avgListVals([5, 3]) == 4.0)
    assert (avgListVals(['5', '3']) == 4.0)
    assert (avgListVals([1, 1]) == 1.0)
    assert (avgListVals([10, 8, 6]) == 8.0)
    print('passed')


def testToNum():
    print('testing toNum')
    assert (toNum('50,000') == float(50000))
    assert (toNum('4,000') == float(4000))
    print('passed')


def testSearch2Lists():
    print('testing search2Lists')
    assert (search2Lists([5, 3, 2], [1, 4, 4], 3, 4) == 1)
    assert (search2Lists([5, 3, 2], [1, 4, 4], 3, 1) == False)
    assert (search2Lists([5, 3, 2], [1, 4, 4], 8, 8) == False)
    assert (search2Lists([5, 3, 2], [1, 4, 4], 5, 1) == 0)
    print('passed')


def testColto1dList():
    print('testing colTo1dList')
    assert (colTo1dList([[3, 2], [5, 4]], 0) == [3, 5])
    assert (colTo1dList([[3, 2], [5, 4]], 1) == [2, 4])
    print('pasesd')


def testAll():
    testAvgListVals()
    testToNum()
    testSearch2Lists()
    testColto1dList()

    # testAll()
    ################################################################################

    # setupGeneratorFleet()
