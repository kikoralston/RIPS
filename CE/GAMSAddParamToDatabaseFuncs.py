# Michael Craig
# October 4, 2016
# Functions for adding parameters to GAMS database. Used for CE & UC models.

import sys
import copy, math
import random
from GAMSAuxFuncs import *
from CalculateOpCost import calcOpCostsTech, calcOpCosts
from AuxFuncs import convertCostToTgtYr, nested_dict_to_dict


################################################################################
##################### CE & UC PARAMETERS #######################################
################################################################################
def addDemandParam(db, demandCEZonal, zoneSet, hourSet, gcmSet, hoursForCE, ipmZones, ipmZoneNums, scaleMWtoGW):
    """
    ADD HOURLY DEMAND parameters (either dict of gcm:zone:hourly demand or dict of zone:demand)
    for both CE and UC models

    :param db:
    :param demandCEZonal:
    :param zoneSet:
    :param hourSet:
    :param gcmSet:
    :param hourSymbols:
    :param ipmZones:
    :param ipmZoneNums:
    :param scaleMWtoGW:
    """

    demandDict = dict()

    if gcmSet is not None:
        for gcm in demandCEZonal:
            demandCENumZoneSymbs = createDictIndexedByZone(demandCEZonal[gcm], ipmZones, ipmZoneNums)
            demandDict[gcm] = getHourly2dParamDict(demandCENumZoneSymbs, hoursForCE[gcm], 1 / scaleMWtoGW)
        list_domain = [gcmSet, zoneSet, hourSet]
    else:
        demandCENumZoneSymbs = createDictIndexedByZone(demandCEZonal, ipmZones, ipmZoneNums)
        demandDict = getHourly2dParamDict(demandCENumZoneSymbs, hoursForCE, 1 / scaleMWtoGW)
        list_domain = [zoneSet, hourSet]

    demandDict = nested_dict_to_dict(demandDict)

    (demandName, demandDescrip) = ('pDemand', 'hourly zonal demand (GWh)')
    demandParam = add_NdParam(db, demandDict, list_domain, demandName, demandDescrip)


def createDictIndexedByZone(dataDict, ipmZones, ipmZoneNums, *args):
    """

    Inputs: dictionary of zone:list, optional scalar. Outputs: dictionary of zone symbol:list.

    :param dataDict:
    :param ipmZones:
    :param ipmZoneNums:
    :param args:
    :return:
    """
    zoneDict = dict()
    for zone in ipmZones:
        if len(args) > 0:
            zoneDict[createZoneSymbol(ipmZoneNums[ipmZones.index(zone)])] = dataDict[zone] * args[0]
        else:
            zoneDict[createZoneSymbol(ipmZoneNums[ipmZones.index(zone)])] = dataDict[zone]
    return zoneDict


def addEguParams(db, genFleet, genSet, genSymbols, ipmZones, ipmZoneNums, scaleLbToShortTon, scaleMWtoGW):
    """ ADD EXISTING GENERATOR PARAMETERS

    :param db: data base object
    :param genFleet: 2d list with gen fleet
    :param genSet: gams set variable with existing generators
    :param genSymbols: symbols of generators
    :param ipmZones: symbols of ipm zones (symbols follow gams model nomenclature)
    :param ipmZoneNums: indexes of ipm zones
    :param scaleLbToShortTon: (float) conversion factor (pounds to short tons)
    :param scaleMWtoGW: (float) conversions factor (MW to GW)
    """
    # Heat rate
    #scalarHrToMmbtuPerMwh = 1 / 1000
    #hrDict = getEguParamDict(genFleet, 'Heat Rate (Btu/kWh)', scalarHrToMmbtuPerMwh * scaleMWtoGW)
    #(hrName, hrDescrip) = ('pHr', 'heat rate (MMBtu/GWh)')
    #hrParam = add1dParam(db, hrDict, genSet, genSymbols, hrName, hrDescrip)

    # Emissions rate ----------------
    emRateDict = getEguParamDict(genFleet, 'CO2EmRate(ton/GWh)', 1)
    emRateRandErrDict = getEguParamDict(genFleet, 'CO2EmRandErr(ton/GWh)', 1)

    # update Emission rate with random error value and divide by 1e3
    for egu in emRateDict:
        if emRateDict[egu] > 0:
            emRateDict[egu] = (emRateDict[egu] + emRateRandErrDict[egu]) / 1e3

    (emRateName, emRateDescrip) = ('pCO2emrate', 'emissions rate (10^3 short ton/GWh)')
    emRateParam = add1dParam(db, emRateDict, genSet, genSymbols, emRateName, emRateDescrip)

    # Zone ----------------
    zoneDict = getEguParamZoneDict(genFleet, 'Region Name', ipmZones, ipmZoneNums)
    (zoneName, zoneDesc) = ('pEguzones', 'zone for each egu')
    zoneParam = add1dParam(db, zoneDict, genSet, genSymbols, zoneName, zoneDesc)


def getEguParamZoneDict(genFleet, zoneCol, ipmZones, ipmZoneNums):
    """ Return dict of {genSymbol:zone num}

    :param genFleet: 2d list with gen fleet
    :param zoneCol:
    :param ipmZones:
    :param ipmZoneNums:
    :return:
    """
    zoneCol = genFleet[0].index(zoneCol)
    zoneDict = dict()
    for row in genFleet[1:]:
        zoneDict[createGenSymbol(row, genFleet[0])] = ipmZoneNums[ipmZones.index(row[zoneCol])]
    return zoneDict


def addEguOpCostParam(db, genFleet, genSet, genSymbols, scaleLbToShortTon, scaleMWtoGW, scaleDollarsToThousands,
                      *co2Price):
    """Add op cost parameter for existing gens

    :param db: data base object
    :param genFleet: 2d list with gen fleet
    :param genSet:
    :param genSymbols:
    :param scaleLbToShortTon:
    :param scaleMWtoGW:
    :param scaleDollarsToThousands:
    :param co2Price:
    """
    ocDict = getEguOpCostDict(genFleet, scaleLbToShortTon, scaleMWtoGW, scaleDollarsToThousands, co2Price)
    (ocName, ocDescrip) = ('pOpcost', 'op cost (thousand$/GWh)')
    ocParam = add1dParam(db, ocDict, genSet, genSymbols, ocName, ocDescrip)


def addEguHourlyParams(db, hourlyCapacsCE, gcmSet, genSet, hourSet, hoursForCE, scaleMWtoGW):
    """ ADD EXISTING GENERATOR HOURLY CAPACITIES

    Add hourly HR & capac param

    :param db:
    :param hourlyCapacsCE:
    :param gcmSet:
    :param genSet:
    :param hourSet:
    :param hoursForCE:
    :param scaleMWtoGW:
    """

    capacDict = dict()

    if gcmSet is not None:
        for gcm in hourlyCapacsCE:
            capacDict[gcm] = getHourly2dParamDict(hourlyCapacsCE[gcm], hoursForCE[gcm], 1 / scaleMWtoGW)
        list_domain = [gcmSet, genSet, hourSet]

    else:
        capacDict = getHourly2dParamDict(hourlyCapacsCE, hoursForCE, 1 / scaleMWtoGW)
        list_domain = [genSet, hourSet]

    capacDict = nested_dict_to_dict(capacDict)

    (capacName, capacDescrip) = ('pCapac', 'hourly capacity (GW)')
    capacParam = add_NdParam(db, capacDict, list_domain, capacName, capacDescrip)


def getHourly2dParamDict(hourlyParam, hoursForCE, scalar):
    """Creates dictionary of (key,hourSymbol):param   (for capac or HR, key=gen; for demand or wind & solar gen,
    key=zone)

    :param hourlyParam:
    :param hoursForCE:
    :param scalar:
    :return:
    """
    paramDict = dict()
    for key in hourlyParam:
        for idx, h in enumerate(hoursForCE):
            hourSymbol = createHourSymbol(h)
            paramDict[(key, hourSymbol)] = hourlyParam[key][idx] * scalar
    return paramDict


def addExistingRenewableMaxGenParams(db, gcmSet, zoneSet, ipmZones, ipmZoneNums, hourSet, hoursForCE,
                                     hourlySolarGenCEZonal, hourlyWindGenCEZonal, scaleMWtoGW):
    """ADD EXISTING RENEWABLE COMBINED MAXIMUM GENERATION VALUES

    Converts 1d list of param vals to hour-indexed dicts, then adds dicts to GAMS db

    :param db:
    :param zoneSet:
    :param ipmZones:
    :param ipmZoneNums:
    :param hourSet:
    :param hourSymbols:
    :param hourlySolarGenCEZonal:
    :param hourlyWindGenCEZonal:
    :param scaleMWtoGW:
    """

    maxSolarGenDict = dict()
    maxWindGenDict = dict()

    if gcmSet is not None:
        for gcm in hourlySolarGenCEZonal:
            hourlySolarGenCEZonalSymbs = createDictIndexedByZone(hourlySolarGenCEZonal[gcm], ipmZones, ipmZoneNums)
            hourlyWindGenCEZonalSymbs = createDictIndexedByZone(hourlyWindGenCEZonal[gcm], ipmZones, ipmZoneNums)
            maxSolarGenDict[gcm] = getHourly2dParamDict(hourlySolarGenCEZonalSymbs, hoursForCE[gcm], 1 / scaleMWtoGW)
            maxWindGenDict[gcm] = getHourly2dParamDict(hourlyWindGenCEZonalSymbs, hoursForCE[gcm], 1 / scaleMWtoGW)
            list_domain = [gcmSet, zoneSet, hourSet]

    else:
        hourlySolarGenCEZonalSymbs = createDictIndexedByZone(hourlySolarGenCEZonal, ipmZones, ipmZoneNums)
        hourlyWindGenCEZonalSymbs = createDictIndexedByZone(hourlyWindGenCEZonal, ipmZones, ipmZoneNums)
        maxSolarGenDict = getHourly2dParamDict(hourlySolarGenCEZonalSymbs, hoursForCE, 1 / scaleMWtoGW)
        maxWindGenDict = getHourly2dParamDict(hourlyWindGenCEZonalSymbs, hoursForCE, 1 / scaleMWtoGW)
        list_domain = [zoneSet, hourSet]


    maxSolarGenDict = nested_dict_to_dict(maxSolarGenDict)
    maxWindGenDict = nested_dict_to_dict(maxWindGenDict)

    (maxSolarGenName, maxSolarGenDescrip) = ('pMaxgensolar', 'max combined gen by existing solar')
    solarParam = add_NdParam(db, maxSolarGenDict, list_domain, maxSolarGenName, maxSolarGenDescrip)

    (maxWindGenName, maxWindGenDescrip) = ('pMaxgenwind', 'max combined gen by existing wind')
    windParam = add_NdParam(db, maxWindGenDict, list_domain, maxWindGenName, maxWindGenDescrip)


# Stores set of values into dictionary keyed by hour
# Inputs: set of param values (1d list), hour symbols (1d list), optional scalar
# Outputs: dictionary of (hour symbol:param val)
def getParamIndexedByHourDict(paramVals, hourSymbols, *scalar):
    paramIndexedByHourDict = dict()
    for idx in range(len(hourSymbols)): paramIndexedByHourDict[hourSymbols[idx]] = paramVals[idx] * scalar[0]
    return paramIndexedByHourDict


##### ADD ZONE AND LINE CONSTRAINTS
# Add parameter mapping lines to zone sources & sinks (pLinesources(l),pLinesinks(l))
def addLineSourceAndSink(db, lineSet, lines, ipmZones, ipmZoneNums):
    lineSources, lineSinks = dict(), dict()
    for line in lines:
        lineSource, lineSink = line.split('_to_')
        lineSources[line] = ipmZoneNums[ipmZones.index(lineSource)]  # convert zone from name to 1+ number
        lineSinks[line] = ipmZoneNums[ipmZones.index(lineSink)]
    sourceName, sourceDesc = 'pLinesources', 'source zone for each line'
    sinkName, sinkDesc = 'pLinesinks', 'sink zone for each line'
    sourceParam = add1dParam(db, lineSources, lineSet, lines, sourceName, sourceDesc)
    sinkParam = add1dParam(db, lineSinks, lineSet, lines, sinkName, sinkDesc)


def addLineCapacs(db, lineCapacs, lineSet, lines, scaleMWtoGW):
    """Add line capacities (pLinecapacs(l)). lineCapacs already a dict of line:capac;

    convert to GW

    :param db:
    :param lineCapacs:
    :param lineSet:
    :param lines:
    :param scaleMWtoGW:
    """
    lineCapacsGW = dict()
    for key in lineCapacs: lineCapacsGW[key] = lineCapacs[key] * 1 / scaleMWtoGW
    name, desc = 'pLinecapacs', 'line capacs (GW)'
    param = add1dParam(db, lineCapacsGW, lineSet, lines, name, desc)


##### ADD PUMPED HYDRO PARAMETERS
def addPumpHydroParams(db, genFleetForCE, phEff, phMaxSoc, phInitSoc, pumpHydroGenSet, pumpHydroGenSymbols,
                       scaleMWtoGW):
    effDict = dict()
    for symb in pumpHydroGenSymbols: effDict[symb] = phEff
    (effname, effdesc) = ('pEfficiency', 'efficiency')
    effParam = add1dParam(db, effDict, pumpHydroGenSet, pumpHydroGenSymbols, effname, effdesc)
    socDict = getSocDict(genFleetForCE, phMaxSoc, pumpHydroGenSymbols, scaleMWtoGW)
    (socname, socdesc) = ('pMaxsoc', 'max state of charge (GWh)')
    socparam = add1dParam(db, socDict, pumpHydroGenSet, pumpHydroGenSymbols, socname, socdesc)
    initSocDict = dict()
    for gen in socDict: initSocDict[gen] = socDict[gen] * phInitSoc  # given as fraction of max SOC
    (socname, socdesc) = ('pInitsoc', 'initial state of charge (GWh)')
    initsocparam = add1dParam(db, initSocDict, pumpHydroGenSet, pumpHydroGenSymbols, socname, socdesc)


# Get state of charge for pumped hydro. phMaxSoc equals multiple of capacity
def getSocDict(genFleetForCE, phMaxSoc, pumpHydroGenSymbols, scaleMWtoGW):
    capacCol = genFleetForCE[0].index('Capacity (MW)')
    socDict = dict()
    for row in genFleetForCE[1:]:
        if createGenSymbol(row, genFleetForCE[0]) in pumpHydroGenSymbols:
            socDict[createGenSymbol(row, genFleetForCE[0])] = float(row[capacCol]) / scaleMWtoGW * phMaxSoc
    return socDict


################################################################################
################################################################################
################################################################################

################################################################################
##################### CAPACITY EXPANSION PARAMETERS ############################
################################################################################
##### ADD NEW TECH PARAMS FOR CE
def addTechParams(db, newTechsCE, scaleMWtoGW, scaleDollarsToThousands, scaleLbToShortTon, ptCurtailed):

    # get sets written to data base
    cellSet = db.get_set('c')
    zoneSet = db.get_set('z')
    techSet = db.get_set('tech')
    techCurtailedSet = db.get_set('techcurtailed')
    techNotCurtailedSet = db.get_set('technotcurtailed')
    renewTechSet = db.get_set('techrenew')
    hourSet = db.get_set('h')

    # get lists with symbols for each set
    cellSymbols = [c.keys[0] for c in cellSet.__iter__()]
    zoneSymbols = [z.keys[0] for z in zoneSet.__iter__()]
    techSymbols = [t.keys[0] for t in techSet.__iter__()]
    techCurtailedSymbols = [t.keys[0] for t in techCurtailedSet.__iter__()]
    techNotCurtailedSymbols = [t.keys[0] for t in techNotCurtailedSet.__iter__()]
    renewTechSymbols = [t.keys[0] for t in renewTechSet.__iter__()]
    hourSymbols = [h.keys[0] for h in hourSet.__iter__()]

    # Nameplate capacity (for cost calculations)
    capacDict = getTechParamDict(newTechsCE, techSymbols, 'Capacity(MW)', ptCurtailed,
                                 1 / scaleMWtoGW)
    (capacName, capacDescrip) = ('pCapactech', 'capacity (GW) of techs')
    techCapacParam = add1dParam(db, capacDict, techSet, techSymbols, capacName, capacDescrip)

    # Op cost ---------------
    # write dictionaries of op. costs for new techs (differentiate costs by zone/cell)

    # 1d list or dictionary with operating costs for each plant type (without adding random error)
    opCosts = getTechOpCostDict(newTechsCE, ptCurtailed, 1)

    # techs that can be curtailed
    dict_techcurt_ocError = dict()
    dict_techcurt_oc = dict()
    for c in cellSymbols:
        for t in techCurtailedSymbols:
            dict_techcurt_ocError[(c, t)] = random.uniform(0, 0.05)
            dict_techcurt_oc[(c, t)] = opCosts[t] + dict_techcurt_ocError[(c, t)]
            dict_techcurt_oc[(c, t)] = dict_techcurt_oc[(c, t)]*(scaleMWtoGW / scaleDollarsToThousands)

    (ocName, ocDescrip) = ('pOpcosttechcurt', 'op cost for curt. tech (thousand$/GWh)')
    ocParam = add_NdParam(db, dict_techcurt_oc, [cellSet, techCurtailedSet], ocName, ocDescrip)

    # techs that cannot be curtailed
    dict_technotcurt_ocError = dict()
    dict_technotcurt_oc = dict()
    for z in zoneSymbols:
        for t in techNotCurtailedSymbols:
            dict_technotcurt_ocError[(z, t)] = random.uniform(0, 0.05)
            dict_technotcurt_oc[(z, t)] = opCosts[t] + dict_technotcurt_ocError[(z, t)]
            dict_technotcurt_oc[(z, t)] = dict_technotcurt_oc[(z, t)]*(scaleMWtoGW / scaleDollarsToThousands)

    (ocName, ocDescrip) = ('pOpcosttechnotcurt', 'op cost for techs not curt. (thousand$/GWh)')
    ocParam = add_NdParam(db, dict_technotcurt_oc, [zoneSet, techNotCurtailedSet], ocName, ocDescrip)

    # renewable techs
    dict_techrenew_ocError = dict()
    dict_techrenew_oc = dict()
    for z in zoneSymbols:
        for t in renewTechSymbols:
            dict_techrenew_ocError[(z, t)] = random.uniform(0, 0.05)
            dict_techrenew_oc[(z, t)] = opCosts[t] + dict_techrenew_ocError[(z, t)]
            dict_techrenew_oc[(z, t)] = dict_techrenew_oc[(z, t)]*(scaleMWtoGW / scaleDollarsToThousands)

    (ocName, ocDescrip) = ('pOpcosttechrenew', 'op cost for renew tech (thousand$/GWh)')
    ocParam = add_NdParam(db, dict_techrenew_oc, [zoneSet, renewTechSet], ocName, ocDescrip)

    # Fixed O&M ---------------
    fixedomDict = getTechParamDict(newTechsCE, techSymbols, 'FOM(2012$/MW/yr)', ptCurtailed,
                                   scaleMWtoGW * 1 / scaleDollarsToThousands)
    for tech in fixedomDict: fixedomDict[tech] = convertCostToTgtYr('fom', fixedomDict[tech])
    (fixedomName, fixedomDescrip) = ('pFom', 'fixed O&M (thousand$/GW/yr)')
    techFixedomParam = add1dParam(db, fixedomDict, techSet, techSymbols, fixedomName, fixedomDescrip)

    # Overnight capital cost ---------------
    occDict = getTechParamDict(newTechsCE, techSymbols, 'CAPEX(2012$/MW)', ptCurtailed,
                               scaleMWtoGW * 1 / scaleDollarsToThousands)
    for tech in occDict: occDict[tech] = convertCostToTgtYr('occ', occDict[tech])
    (occName, occDescrip) = ('pOcc', 'overnight capital cost (thousand$/GW)')
    techOccParam = add1dParam(db, occDict, techSet, techSymbols, occName, occDescrip)

    # Emissions rate ---------------
    # write dictionaries of co2 emission rates for new techs (differentiate rates by zone/cell)
    co2EmRates = getTechParamDict(newTechsCE, techSymbols, 'CO2EmRate(ton/GWh)', ptCurtailed)

    # techs that can be curtailed
    dict_techcurt_emRateError = dict()
    dict_techcurt_emRate = dict()
    for c in cellSymbols:
        for t in techCurtailedSymbols:
            dict_techcurt_emRateError[(c, t)] = random.uniform(0, 0.05)
            if co2EmRates[t] == 0:
                dict_techcurt_emRate[(c, t)] = 0
            else:
                dict_techcurt_emRate[(c, t)] = (co2EmRates[t] + dict_techcurt_emRateError[(c, t)]) / 1e3

    (emRateName, emRateDescrip) = ('pCO2emratetechcurt',
                                   'co2 emissions rate for curtailed techs (10^3 short ton/GWh)')
    techEmRateParam = add_NdParam(db, dict_techcurt_emRate, [cellSet, techCurtailedSet], emRateName, emRateDescrip)

    # techs that cannot be curtailed
    dict_technotcurt_emRateError = dict()
    dict_technotcurt_emRate = dict()
    for z in zoneSymbols:
        for t in techNotCurtailedSymbols:
            dict_technotcurt_emRateError[(z, t)] = random.uniform(0, 0.05)
            if co2EmRates[t] == 0:
                dict_technotcurt_emRate[(z, t)] = 0
            else:
                dict_technotcurt_emRate[(z, t)] = (co2EmRates[t] + dict_technotcurt_emRateError[(z, t)])/1e3

    (emRateName, emRateDescrip) = ('pCO2emratetechnotcurt',
                                   'co2 emissions rate for techs not curtailed (10^3 short ton/GWh)')
    techEmRateParam = add_NdParam(db, dict_technotcurt_emRate, [zoneSet, techNotCurtailedSet], emRateName, emRateDescrip)

    # renewable techs
    dict_techrenew_emRateError = dict()
    dict_techrenew_emRate = dict()
    for z in zoneSymbols:
        for t in renewTechSymbols:
            dict_techrenew_emRateError[(z, t)] = random.uniform(0, 0.05)
            if co2EmRates[t] == 0:
                dict_techrenew_emRate[(z, t)] = 0
            else:
                dict_techrenew_emRate[(z, t)] = (co2EmRates[t] + dict_techrenew_emRateError[(z, t)]) / 1e3

    (emRateName, emRateDescrip) = ('pCO2emratetechrenew',
                                   'co2 emissions rate for renewable techs (10^3 short ton/GWh)')
    techEmRateParam = add_NdParam(db, dict_techrenew_emRate, [zoneSet, renewTechSet], emRateName, emRateDescrip)

    # Lifetime ---------------------
    lifetimeDict = getTechParamDict(newTechsCE, techSymbols, 'Lifetime(years)', ptCurtailed)
    (lifetimeName, lifetimeDescrip) = ('pLife', 'years')
    techLifetimeParam = add1dParam(db, lifetimeDict, techSet, techSymbols, lifetimeName, lifetimeDescrip)


def getTechParamDict(newTechsCE, techSymbols, paramColName, ptCurtailed, *scalar):
    """Creates dict of (techSymbol:paramVal) for given parameter name

    :param newTechsCE:
    :param techSymbols:
    :param paramColName:
    :param ptCurtailed:
    :param scalar:
    :return:
    """
    techCol = newTechsCE[0].index('TechnologyType')
    paramCol = newTechsCE[0].index(paramColName)
    techSymbolsInNewTechsCE = [createTechSymbol(row, newTechsCE[0], ptCurtailed) for row in newTechsCE]
    paramDict = dict()

    for techSymbol in techSymbols:
        rowIdx = techSymbolsInNewTechsCE.index(techSymbol)
        if len(scalar) > 0:
            paramDict[techSymbol] = float(newTechsCE[rowIdx][paramCol]) * scalar[0]
        else:
            paramDict[techSymbol] = float(newTechsCE[rowIdx][paramCol])

    return paramDict


def getTechOpCostDict(newTechs, ptCurtailed, scalar):
    """Takes in techs and returns dictionary of (tech:opCost)

    :param newTechs:
    :param ptCurtailed:
    :param scalar:
    :return:
    """
    opCosts = calcOpCostsTech(newTechs)
    paramDict = dict()
    for idx in range(1, len(newTechs)):
        # op costs = 1d list of vals, so offset by 1
        paramDict[createTechSymbol(newTechs[idx], newTechs[0], ptCurtailed)] = opCosts[idx-1]*scalar
    return paramDict


def addTechCurtailedHourlyCapac(db, hourlyCurtailedTechCapacsCE, gcmSet, cellSet, techCurtailedSet, hourSet, hoursForCE,
                                scaleMWtoGW):
    """Add hourly curtailed tech capac.

    Input capacs: (tech,loc):[capacs]. Output dict added to db: idxed by (loc,tech,hr)

    :param db:
    :param hourlyCurtailedTechCapacsCE:
    :param cellSet:
    :param techCurtailedSet:
    :param hourSet:
    :param hoursForCE:
    :param scaleMWtoGW:
    """
    capacDict = dict()
    for gcm in hourlyCurtailedTechCapacsCE:
        capacDict[gcm] = getHourlyTechParamDict(hourlyCurtailedTechCapacsCE[gcm], hoursForCE[gcm], 1 / scaleMWtoGW)

    capacDict = nested_dict_to_dict(capacDict)

    (cName, cDes) = ('pCapactechcurtailed', 'hourly curtailed capac (GW)')
    capParam = add_NdParam(db, capacDict, [gcmSet, cellSet, techCurtailedSet, hourSet], cName, cDes)


def getHourlyTechParamDict(hourlyParam, hoursForCE, scalar):
    """Creates dictionary of (gcm, locSymbol, techSymbol, hourSymbol):capac or HR


    :param hourlyParam:
    :param hoursForCE:
    :param scalar:
    :return:
    """
    paramDict = dict()
    for (techSymbol, locSymbol) in hourlyParam:
        for idx, h in enumerate(hoursForCE):
            hourSymbol = createHourSymbol(h)
            paramDict[(locSymbol, techSymbol, hourSymbol)] = hourlyParam[(techSymbol, locSymbol)][idx] * scalar
    return paramDict


##### ADD MAP FROM CELLS TO ZONES
def addCellsToZones(db, cellSet, cellsToZones, ipmZones, ipmZoneNums):
    cellsToZoneNums, cellSymbols = dict(), list()
    for cell in cellsToZones:
        cellsToZoneNums[cell] = ipmZoneNums[ipmZones.index(cellsToZones[cell])]
        cellSymbols.append(cell)
    (name, desc) = ('pCellzones', 'zone each cell is in')
    cToZParam = add1dParam(db, cellsToZoneNums, cellSet, cellSymbols, name, desc)


##### ADD PLANNING RESERVE MARGIN FRACTION PARAMETER
# Add zonal planning reserve
def addPlanningReserveParam(db, planningReserveZonal, ipmZones, ipmZoneNums, zoneSet, zoneSymbols, scaleMWtoGW):
    reserveDict = createDictIndexedByZone(planningReserveZonal, ipmZones, ipmZoneNums, 1 / scaleMWtoGW)
    planName, planDesc = 'pPlanningreserve', 'planning reserve'
    add1dParam(db, reserveDict, zoneSet, zoneSymbols, planName, planDesc)


# Add map of peak hour to zone
def addPeakHourToZoneParam(db, peakDemandHourZonal, peakHourSet, peakHrSymbols, ipmZones, ipmZoneNums):
    peakDict = dict()
    for zone in peakDemandHourZonal:
        peakDict[createHourSymbol(peakDemandHourZonal[zone])] = ipmZoneNums[ipmZones.index(zone)]
    (peakName, peakDesc) = ('pPeakhtozone', 'map peak hours to zones')
    zoneParam = add1dParam(db, peakDict, peakHourSet, peakHrSymbols, peakName, peakDesc)


##### ADD DISCOUNT RATE PARAMETER
def addDiscountRateParam(db, discountRate):
    add0dParam(db, 'pR', 'discount rate', discountRate)


##### ADD FIRM FRACTION FOR EXISTING GENERATORS
# Firm fraction goes towards meeting planning reserve margin
def addExistingPlantFirmFractions(db, genFleet, genSet, genSymbols, firmCapacityCreditsExistingGens):
    firmCreditDict = getFirmCreditExistingGenDict(genFleet, firmCapacityCreditsExistingGens)
    (firmCreditName, firmCreditDescrip) = ('pFirmcapacfractionegu', 'firm capacity fraction')
    firmCreditExistingGenParam = add1dParam(db, firmCreditDict, genSet, genSymbols, firmCreditName, firmCreditDescrip)


# Returns dict of (genSymbol:capacCredit) based on plant type of each generator
def getFirmCreditExistingGenDict(genFleet, firmCapacityCreditsExistingGens):
    plantTypeCol = genFleet[0].index('PlantType')
    firmCapacityCreditsExistingGensDict = dict()
    for row in genFleet[1:]:
        capacCredit = firmCapacityCreditsExistingGens[row[plantTypeCol]]
        firmCapacityCreditsExistingGensDict[createGenSymbol(row, genFleet[0])] = capacCredit
    return firmCapacityCreditsExistingGensDict


def addRenewTechCFParams(db, renewTechSet, renewTechSymbols, gcmSet, zoneSet, hourSet, hoursForCE, newWindCFsCEZonal,
                         newSolarCFsCEZonal, ipmZones, ipmZoneNums):
    """ADD HOURLY CAPACITY FACTORS FOR NEW RENEWABLE TECHS

    Add CFs for new renew builds. Input: CFs as zone:[CF]

    Output: dict added to GAMS file as (z,tech,h):[CFs]

    :param db:
    :param renewTechSet:
    :param renewTechSymbols:
    :param zoneSet:
    :param hourSet:
    :param hourSymbols:
    :param newWindCFsCEZonal:
    :param newSolarCFsCEZonal:
    :param ipmZones:
    :param ipmZoneNums:
    """
    renewtechCfDict = dict()
    for gcm in newWindCFsCEZonal:
        for renewtech in renewTechSymbols:
            if renewtech == 'Wind':
                relevantCfs = copy.deepcopy(newWindCFsCEZonal[gcm])
            elif renewtech == 'Solar PV':
                relevantCfs = copy.deepcopy(newSolarCFsCEZonal[gcm])

            for zone in relevantCfs:
                for idx, h in enumerate(hoursForCE[gcm]):
                    renewtechCfDict[(gcm, createZoneSymbol(ipmZoneNums[ipmZones.index(zone)]),
                                     renewtech, createHourSymbol(h))] = relevantCfs[zone][idx]

    (renewtechCFName, renewtechCFDescrip) = ('pCf', 'capacity factors for new wind and solar')

    renewtechCfParam = add_NdParam(db, renewtechCfDict, [gcmSet, zoneSet, renewTechSet, hourSet], renewtechCFName,
                                   renewtechCFDescrip)


def addCppEmissionsCap(db, co2CppSercCurrYearLimit):
    """ADD CO2 EMISSIONS CAP FOR CE MODEL

    :param db: data base object
    :param co2CppSercCurrYearLimit: annual co2 upper bound (in short tons)
    """
    add0dParam(db, 'pCO2emcap', 'CPP co2 emissions cap [10^3 short tons]', co2CppSercCurrYearLimit/1e3)


def addSeasonDemandWeights(db, seasonDemandWeights):
    """ADD WEIGHTS TO SCALE REPRESENTATIVE SEASONAL DEMAND UP TO COMPLETE SEASON PERIOD

    :param db: data base object
    :param seasonDemandWeights: dictionary with season weights
    """
    for season in seasonDemandWeights:
        add0dParam(db, 'pWeight' + season, 'weight on rep. seasonal demand', seasonDemandWeights[season])


def addMaxNumNewBuilds(db, newTechsCE, zoneSet, ipmZones, ipmZoneNums, typeSet, maxAddedZonalCapacPerTech, ptCurtailed):
    """ ADD LIMIT ON MAX NUMBER OF NEW BUILDS PER PLANT TYPE BY ZONE

    :param db:
    :param newTechsCE:
    :param zoneSet:
    :param ipmZones:
    :param ipmZoneNums:
    :param typeSet:
    :param maxAddedZonalCapacPerTech:
    :param ptCurtailed:
    """
    capacCol = newTechsCE[0].index('Capacity(MW)')
    techCol = newTechsCE[0].index('TechnologyType')

    techMaxNewBuildsDict = dict()

    listTypes = [x.get_keys()[0] for x in typeSet]
    d = None

    if isinstance(maxAddedZonalCapacPerTech, int) or isinstance(maxAddedZonalCapacPerTech, float):
        d = dict()
        for ty in listTypes:
            d[ty] = maxAddedZonalCapacPerTech
    elif isinstance(maxAddedZonalCapacPerTech, dict):
        d = maxAddedZonalCapacPerTech
    else:
        print('------------------------------------------------------------------------------------------')
        print('ERROR!!!')
        print('Parameter maxAddedZonalCapacPerTech must be either a number or a dictionary')
        print('------------------------------------------------------------------------------------------')
        sys.exit()

    for ty in listTypes:
        # assumes ALL cooling techs in the same plant type have the same capacity size
        # iterate over list of new techs, find those that match plant type and return only the
        # first occurrence
        capacType = [int(row[capacCol]) for row in newTechsCE[1:] if row[techCol] == ty][0]

        for zone in ipmZones:
            ubound = d[ty]

            if ubound == float('inf'):
                ubound = sys.maxsize

            techMaxNewBuildsDict[(createZoneSymbol(ipmZoneNums[ipmZones.index(zone)]), ty)] = \
                math.ceil(ubound / capacType)

    (maxBuildName, maxBuildDescrip) = ('pNmax', 'max num builds per type of plant')
    maxBuildParam = add2dParam(db, techMaxNewBuildsDict, zoneSet, typeSet, maxBuildName, maxBuildDescrip)


def addHydroMaxGenPerSeason(db, hydroGenSet, gcmSet, hydroPotPerSeason, scaleMWtoGW):
    """ ADD MAX HYDRO GEN PER TIME BLOCK

    hydroPotPerSeason is a dict of season:genSymbol:gen (MWh), so just need to scale it.
    Note that sesaon can be 'special'!

    :param db:
    :param hydroGenSet:
    :param hydroGenSymbols:
    :param hydroPotPerSeason:
    :param scaleMWtoGW:
    """

    seasonList = hydroPotPerSeason[list(hydroPotPerSeason.keys())[0]].keys()

    for season in seasonList:
        maxGenDict = dict()
        for gcm in hydroPotPerSeason:
            auxDict = dict()
            for symb, pot in hydroPotPerSeason[gcm][season].items():
                auxDict[symb] = pot / scaleMWtoGW
            maxGenDict[gcm] = auxDict

        maxGenDict = nested_dict_to_dict(maxGenDict)

        name, desc = 'pMaxhydrogen' + season[:3], 'max gen by each hydro unit'  # ex: pMaxhydrogensum
        maxGenParam = add_NdParam(db, maxGenDict, [gcmSet, hydroGenSet], name, desc)


def addHydroMaxGenUC(db, hydroGenSet, hydroPot, scaleMWtoGW):
    """ ADD MAX HYDRO GEN FOR DAYS IN UC RUN

    hydroPot is a dict of genSymbol:gen (MWh), so just need to scale it.

    :param db:
    :param hydroGenSet:
    :param hydroPot:
    :param scaleMWtoGW:
    """

    maxGenDict = {genId: hydroPot[genId]/scaleMWtoGW for genId in hydroPot}

    name, desc = 'pMaxgenhydro', 'max gen by each hydro unit'  # ex: pMaxhydrogensum
    maxGenParam = add_NdParam(db, maxGenDict, [hydroGenSet], name, desc)


################################################################################
##################### UNIT COMMITMENT PARAMETERS ###############################
################################################################################
# Add UC parameters
def addEguUCParams(db, fleetUC, genSet, genSymbols, scaleMWtoGW, scaleDollarsToThousands):
    # Min load
    minLoadDict = getEguParamDict(fleetUC, 'MinLoad(MW)', 1 / scaleMWtoGW)
    (minLoadName, minLoadDescrip) = ('pMinload', 'min load (GW)')
    minLoadParam = add1dParam(db, minLoadDict, genSet, genSymbols, minLoadName, minLoadDescrip)
    # Ramp rate
    rampDict = getEguParamDict(fleetUC, 'RampRate(MW/hr)', 1 / scaleMWtoGW)
    (rampName, rampDescrip) = ('pRamprate', 'ramp rate (GW/hr)')
    rampParam = add1dParam(db, rampDict, genSet, genSymbols, rampName, rampDescrip)
    # Start up fixed cost
    startCostDict = getEguParamDict(fleetUC, 'StartCost($)', 1 / scaleDollarsToThousands)
    (startName, startDescrip) = ('pStartupfixedcost', 'startup fixed cost (thousand$)')
    startCostParam = add1dParam(db, startCostDict, genSet, genSymbols, startName, startDescrip)
    # Min down time
    minDownDict = getEguParamDict(fleetUC, 'MinDownTime(hrs)', 1)
    (minDownName, minDownDescrip) = ('pMindowntime', 'min down time (hrs)')
    minDownParam = add1dParam(db, minDownDict, genSet, genSymbols, minDownName, minDownDescrip)


# Add non-time-varying capacity param for generators
# def addEguCapacParam(db,genFleet,genSet,genSymbols,scaleMWtoGW):
#     capacDict = getEguParamDict(genFleet,'Capacity (MW)',1/scaleMWtoGW)
#     (capacName,capacDescrip) = ('pCapac','capacity (GW)')
#     capacParam = add1dParam(db,capacDict,genSet,genSymbols,capacName,capacDescrip)

##### RESERVE PARAMETERS
# Add reg reserve parameters
def addRegReserveParameters(db, regUp, regDown, rrToRegTime, hourSet, hourSymbols, zoneSet, modelName, genparam):

    rampToRegParam = db.add_parameter('pRampratetoregreservescalar', 0, 'convert ramp rate to reg timeframe')
    rampToRegParam.add_record().value = rrToRegTime

    # Add hourly reg reserves; in CE model, increases w/ built wind, hence diff
    # name than UC model.
    regUpNumZoneSymbs = createDictIndexedByZone(regUp, genparam.ipmZones, genparam.ipmZoneNums)
    regUpDict = getHourly2dParamDict(regUpNumZoneSymbs, hourSymbols, 1 / genparam.scaleMWtoGW)
    #regUpDict = getParamIndexedByHourDict(regUp, hourSymbols, 1 / scaleMWtoGW)

    if modelName == 'UC':
        regParamName = 'pRegupreserves'
    elif modelName == 'CE':
        regParamName = 'pRegupreserveinitial'
    (regUpName, regUpDescr) = (regParamName, 'hourly reg up reserves (GWh)')
    regParam = add_NdParam(db, regUpDict, [zoneSet, hourSet], regUpName, regUpDescr)

    regDownNumZoneSymbs = createDictIndexedByZone(regDown, genparam.ipmZones, genparam.ipmZoneNums)
    regDownDict = getHourly2dParamDict(regDownNumZoneSymbs, hourSymbols, 1 / genparam.scaleMWtoGW)
    #regDownDict = getParamIndexedByHourDict(regDown, hourSymbols, 1 / scaleMWtoGW)
    if modelName == 'UC':
        regParamName = 'pRegdownreserves'
    elif modelName == 'CE':
        regParamName = 'pRegdownreserveinitial'  # not used in current project

    (regDownName, regDownDescr) = (regParamName, 'hourly reg down reserves (GWh)')
    regParam = add_NdParam(db, regDownDict, [zoneSet, hourSet], regDownName, regDownDescr)


# Add reserve parameter quantities
def addFlexReserveParameters(db, flexRes, rrToFlexTime, hourSet, hourSymbols, zoneSet, modelName, genparam):

    rampToFlexParam = db.add_parameter('pRampratetoflexreservescalar', 0, 'convert ramp rate to flex timeframe')
    rampToFlexParam.add_record().value = rrToFlexTime

    flexNumZoneSymbs = createDictIndexedByZone(flexRes, genparam.ipmZones, genparam.ipmZoneNums)
    flexDict = getHourly2dParamDict(flexNumZoneSymbs, hourSymbols, 1 / genparam.scaleMWtoGW)
    #flexDict = getParamIndexedByHourDict(flexRes, hourSymbols, 1 / scaleMWtoGW)
    if modelName == 'UC':
        regParamName = 'pFlexreserves'
    elif modelName == 'CE':
        regParamName = 'pFlexreserveinitial'
    (flexName, flexDesc) = (regParamName, 'hourly flex reserves (GWh)')
    flexParam = add_NdParam(db, flexDict, [zoneSet, hourSet], flexName, flexDesc)


def addContReserveParameters(db, contRes, rrToContTime, hourSet, hourSymbols, zoneSet, genparam):

    rampToContParam = db.add_parameter('pRampratetocontreservescalar', 0, 'convert ramp rate to cont timeframe')
    rampToContParam.add_record().value = rrToContTime

    contNumZoneSymbs = createDictIndexedByZone(contRes, genparam.ipmZones, genparam.ipmZoneNums)
    contDict = getHourly2dParamDict(contNumZoneSymbs, hourSymbols, 1 / genparam.scaleMWtoGW)
    #contDict = getParamIndexedByHourDict(contRes, hourSymbols, 1 / scaleMWtoGW)

    (contName, contDesc) = ('pContreserves', 'hourly cont reserves (GWh)')
    contParam = add_NdParam(db, contDict, [zoneSet, hourSet], contName, contDesc)


# Add initial conditions
def addEguInitialConditions(db, genSet, genSymbols, fleetUC, onOffInitial, genAboveMinInitial,
                            mdtCarriedInitial, scaleMWtoGW):
    onOffInitialDict = getInitialCondsDict(fleetUC, onOffInitial, 1)
    (onOffInitialName, onOffInitialDescrip) = (
    'pOnoroffinitial', 'whether initially on (1) or off (0) based on last UC')
    onOffInitialParam = add1dParam(db, onOffInitialDict, genSet, genSymbols, onOffInitialName, onOffInitialDescrip)
    mdtCarryDict = getInitialCondsDict(fleetUC, mdtCarriedInitial, 1)
    (mdtCarryName, mdtCarryDescrip) = ('pMdtcarriedhours', 'remaining min down time hrs from last UC (hrs))')
    mdtCarryParam = add1dParam(db, mdtCarryDict, genSet, genSymbols, mdtCarryName, mdtCarryDescrip)
    genAboveMinDict = getInitialCondsDict(fleetUC, genAboveMinInitial, scaleMWtoGW)
    (genAboveMinName, genAboveMinDescrip) = ('pGenabovemininitial', 'initial gen above min load based on last UC (GW)')
    genAboveMinParam = add1dParam(db, genAboveMinDict, genSet, genSymbols, genAboveMinName, genAboveMinDescrip)


def getInitialCondsDict(fleetUC, initialCondValues, *scalar):
    initCondsDict = dict()
    for rowNum in range(1, len(fleetUC)):
        initCondsDict[createGenSymbol(fleetUC[rowNum], fleetUC[0])] = initialCondValues[rowNum - 1] * scalar[0]
    return initCondsDict


##### WHICH GENERATORS ARE ELIGIBLE TO PROVIDE RESERVES
# Add parameter for which existing generators can provide flex, cont, or reg reserves
def addEguEligibleToProvideRes(db, fleetUC, genSet, genSymbols, *stoMarket):

    fleetOrTechsFlag = 'fleet'

    eligibleFlexDict = getEligibleSpinDict(fleetUC, fleetOrTechsFlag, stoMarket)
    (eligFlexName, eligFlexDesc) = ('pFlexeligible', 'egu eligible to provide flex (1) or not (0)')
    eligFlexParam = add1dParam(db, eligibleFlexDict, genSet, genSymbols, eligFlexName, eligFlexDesc)

    eligContDict = getEligibleSpinDict(fleetUC, fleetOrTechsFlag, stoMarket)
    (eligContName, eligContDesc) = ('pConteligible', 'egu eligible to provide cont (1) or not (0)')
    eligContParam = add1dParam(db, eligContDict, genSet, genSymbols, eligContName, eligContDesc)

    eligibleRegDict = getEguParamDict(fleetUC, 'RegOfferElig', 1)
    (eligRegName, eligRegDescrip) = ('pRegeligible', 'egu eligible to provide reg res (1) or not (0)')
    eligibleRegParam = add1dParam(db, eligibleRegDict, genSet, genSymbols, eligRegName, eligRegDescrip)


# Returns dict of whether units can provide spin reserves or not based on the plant type
def getEligibleSpinDict(fleetOrTechsData, fleetOrTechsFlag, *stoMktOrPtCurt):

    (windPlantType, solarPlantType) = getWindAndSolarPlantTypes()
    plantTypesNotProvideRes = {windPlantType, solarPlantType}

    if fleetOrTechsFlag == 'fleet':
        if len(stoMktOrPtCurt) > 0:
            if len(stoMktOrPtCurt[0]) > 0 and stoMktOrPtCurt[0][0] == 'energy':
                plantTypesNotProvideRes.add('Storage')
        plantTypeCol = fleetOrTechsData[0].index('PlantType')

    elif fleetOrTechsFlag == 'techs':
        if len(stoMktOrPtCurt) > 0 and len(stoMktOrPtCurt[0]) > 0:
            ptCurtailed = stoMktOrPtCurt[0][0]
            print('check this eligible spin dict in GAMSAddParam!')
        plantTypeCol = fleetOrTechsData[0].index('TechnologyType')

    eligibleSpinDict = dict()
    for rowNum in range(1, len(fleetOrTechsData)):
        plantType = fleetOrTechsData[rowNum][plantTypeCol]

        if plantType in plantTypesNotProvideRes:
            provideSpin = 0
        else:
            provideSpin = 1

        if fleetOrTechsFlag == 'fleet':
            symbol = createGenSymbol(fleetOrTechsData[rowNum], fleetOrTechsData[0])
        elif fleetOrTechsFlag == 'techs':
            symbol = createTechSymbol(fleetOrTechsData[rowNum], fleetOrTechsData[0],
                                      ptCurtailed)

        eligibleSpinDict[symbol] = provideSpin

    return eligibleSpinDict


def getWindAndSolarPlantTypes():
    return ('Wind', 'Solar PV')


# Add cost of CNSE
def addCostNonservedEnergy(db, cnse, scaleMWtoGW, scaleDollarsToThousands):
    add0dParam(db, 'pCnse', 'cost of non-served energy (thousand$/GWh)',
               cnse * scaleMWtoGW * 1 / scaleDollarsToThousands)


# Add CO2 price
def addCo2Price(db, co2Price, scaleDollarsToThousands):
    add0dParam(db, 'pCO2price', 'co2 emissions price (thousand$/short ton)',
               co2Price * 1 / scaleDollarsToThousands)


################################################################################
################################################################################
################################################################################

################################################################################
############ GENERIC FUNCTIONS TO ADD PARAMS TO GAMS DB ########################
################################################################################
def add0dParam(db, paramName, paramDescrip, paramValue):
    addedParam = db.add_parameter(paramName, 0, paramDescrip)
    addedParam.add_record().value = paramValue


def add1dParam(db, paramDict, idxSet, setSymbols, paramName, paramDescrip):
    addedParam = db.add_parameter_dc(paramName, [idxSet], paramDescrip)
    for idx in setSymbols:
        addedParam.add_record(idx).value = paramDict[idx]
    return addedParam


def add2dParam(db, param2dDict, idxSet1, idxSet2, paramName, paramDescrip):
    addedParam = db.add_parameter_dc(paramName, [idxSet1, idxSet2], paramDescrip)
    for k, v in iter(param2dDict.items()):
        addedParam.add_record(k).value = v
    return addedParam


def add_NdParam(db, paramDict, list_idxSet, paramName, paramDescrip):
    """ Generic function to add a N-dimensional parameter to database

    paramDict is a simple dictionary (it MUST NOT BE a nested dictionary). See function 'nested_dict_to_dict()'
    to convert a nested dictionary into a simple dictionary where the key is a tuple with combination of nested keys

    :param db: database object
    :param paramDict: simple dictionary with data. keys must be a tuple with N values
    :param list_idxSet: 1d list of size N with GAMS sets for domain of parameters (must be in the correct orders)
    :param paramName: (string)
    :param paramDescrip: (string)
    :return:
    """
    addedParam = db.add_parameter_dc(paramName, list_idxSet, paramDescrip)
    for k, v in iter(paramDict.items()):
        addedParam.add_record(k).value = v
    return addedParam


def add3dParam(db, param3dDict, idxSet1, idxSet2, idxSet3, paramName, paramDescrip):
    addedParam = db.add_parameter_dc(paramName, [idxSet1, idxSet2, idxSet3], paramDescrip)
    for k, v in iter(param3dDict.items()):
        addedParam.add_record(k).value = v
    return addedParam


# Takes in gen fleet and param col name, and returns a dictionary of (genSymbol:paramVal)
# for each row.
def getEguParamDict(genFleet, paramColName, *scalar):
    paramCol = genFleet[0].index(paramColName)
    paramDict = dict()

    for row in genFleet[1:]:
        paramDict[createGenSymbol(row, genFleet[0])] = float(row[paramCol]) * scalar[0]

    return paramDict


# Takes in gen fleet and returns dictionary of (genSymbol:opCost)
def getEguOpCostDict(genFleet, scaleLbToShortTon, scaleMWtoGW, scaleDollarsToThousands, *co2Price):
    if len(co2Price[0]) > 0:
        (opCosts, hrs) = calcOpCosts(genFleet, scaleLbToShortTon, co2Price[0][0])  # thousand $/GWh
    else:
        (opCosts, hrs) = calcOpCosts(genFleet, scaleLbToShortTon)  # thousand $/GWh
    paramDict = dict()
    for idx in range(1, len(genFleet)):
        genSymb = createGenSymbol(genFleet[idx], genFleet[0])
        paramDict[genSymb] = opCosts[
                                 idx - 1] * scaleMWtoGW / scaleDollarsToThousands  # op costs = 1d list of vals, so offset by 1
    return paramDict
################################################################################
################################################################################
################################################################################

# Add hourly HR params to tech.
# Input HRs: (tech,loc):[HR]. Output dict added to db: idxed by (loc,tech,hr)
# def addTechHourlyHR(db,hourlyTechHrsCE,cellSet,techSet,hourSet,hourSymbols,scaleMWtoGW):
#     hrDict = getHourlyTechParamDict(hourlyTechHrsCE,hourSymbols,1)
#     (hrName,hrDescrip) = ('pHrtech','heat rate (MMBtu/GWh)')
#     techHrParam = add3dParam(db,hrDict,cellSet,techSet,hourSet,hrName,hrDescrip)
