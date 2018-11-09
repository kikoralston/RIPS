# Michael Craig
# Jan 5, 2017

from AuxFuncs import *
from GAMSAuxFuncs import *
from SetupResultLists import (setupHourlyGenByPlant, setupHourlySystemResultsWithHourSymbols,
                              setupHourlyLineflow, setupEmptyHourly2dList)
import copy, csv
import pandas as pd


def saveCapacExpOperationalData(ceModel):
    """Saves operational CE results: operations by plants & new techs, and MCs on demand and reserve constraints.

    :param ceModel:
    :return:
    """

    (genByPlant, genByCTechAndCell, genByRETechAndZone, genByNCTechAndZone) = saveGeneratorSpecificResults(ceModel)

    (chargeByPH, socByPH) = savePumpedHydroResults(ceModel)

    flowByLine = saveLineFlows(ceModel)

    #sysResults = saveSystemResults(ceModel, hoursForCESymbols, ipmZoneNums)
    #co2EmAndCostResults = saveCo2EmAndCostResults(ceModel)

    return genByPlant, genByCTechAndCell, genByRETechAndZone, genByNCTechAndZone, flowByLine, chargeByPH, socByPH


def saveGeneratorSpecificResults(ceModel):
    """SAVE GENERATOR SPECIFIC RESULTS

    All result arrays need to have same exact format & row/col idxs for gens/hours

    :param ceModel:
    :param genFleetForCE:
    :param curtailedTechs:
    :param renewTechs:
    :param notCurtailedTechs:
    :param cellsForTechs:
    :param ipmZoneNums:
    :param hoursForCESymbols:
    :return:
    """

    d = dict()
    for rec in ceModel.out_db['vPegu']:
        d[tuple(rec.keys)] = rec.get_level()

    genByPlant = dict_tuples_to_list2d(d)

    # use pandas to sort resulting list
    genByPlant = pd.DataFrame(genByPlant, columns=['gcm', 'genId', 'hour', 'value'])
    genByPlant = genByPlant.sort_values(by=['gcm', 'hour', 'genId']).reset_index(drop=True)
    genByPlant = genByPlant.values.tolist()

    # for curtailed candidate techs, the size of the list may be too long and most of the candidates will
    # not be chosen. So filter only to the chosen candidates.

    # get build decisions for candidate techs
    d = dict()
    for rec in ceModel.out_db['vNcurtailed']:
        d[tuple(rec.keys)] = rec.get_level()

    # list with (cell, tech) of curtailed techs that were built
    build_techs = [k for k in d if d[k] > 0]

    # populate dictionary with generation by built candidates
    daux = dict()
    for rec in ceModel.out_db['vPtechcurtailed']:
        if (rec.key(1), rec.key(2)) in build_techs:
            daux[tuple(rec.keys)] = rec.get_level()

    genByCTechAndCell = dict_tuples_to_list2d(daux)

    # use pandas to sort resulting list
    genByCTechAndCell = pd.DataFrame(genByCTechAndCell, columns=['gcm', 'cell', 'tech', 'hour', 'value'])
    genByCTechAndCell = genByCTechAndCell.sort_values(by=['gcm', 'hour', 'tech', 'cell']).reset_index(drop=True)
    genByCTechAndCell = genByCTechAndCell.values.tolist()

    d = dict()
    for rec in ceModel.out_db['vPtechrenew']:
        d[tuple(rec.keys)] = rec.get_level()

    genByRETechAndZone = dict_tuples_to_list2d(d)

    # use pandas to sort resulting list
    genByRETechAndZone = pd.DataFrame(genByRETechAndZone, columns=['gcm', 'zone', 'tech', 'hour', 'value'])
    genByRETechAndZone = genByRETechAndZone.sort_values(by=['gcm', 'hour', 'tech', 'zone']).reset_index(drop=True)
    genByRETechAndZone = genByRETechAndZone.values.tolist()

    d = dict()
    for rec in ceModel.out_db['vPtechnotcurtailed']:
        d[tuple(rec.keys)] = rec.get_level()

    genByNCTechAndZone = dict_tuples_to_list2d(d)

    # use pandas to sort resulting list
    genByNCTechAndZone = pd.DataFrame(genByNCTechAndZone, columns=['gcm', 'zone', 'tech', 'hour', 'value'])
    genByNCTechAndZone = genByNCTechAndZone.sort_values(by=['gcm', 'hour', 'tech', 'zone']).reset_index(drop=True)
    genByNCTechAndZone = genByNCTechAndZone.values.tolist()

    return genByPlant, genByCTechAndCell, genByRETechAndZone, genByNCTechAndZone


def setupCETechResultsLists(curtailedTechs, renewTechs, notCurtailedTechs, cellsForTechs,
                            ipmZoneNums, hoursForCESymbols):
    """Setup same-formatted results lists as above for by-tech results

    :param curtailedTechs:
    :param renewTechs:
    :param notCurtailedTechs:
    :param cellsForTechs:
    :param ipmZoneNums:
    :param hoursForCESymbols:
    :return:
    """
    (genByCTechAndCell, cTechAndCellToRow, hourToCol) = setupHourlyGenByTechAndLocs(cellsForTechs, 'cells',
                                                                                    curtailedTechs, hoursForCESymbols)
    (genByRETechAndZone, reTechAndZoneToRow, hourToCol) = setupHourlyGenByTechAndLocs(ipmZoneNums, 'zones',
                                                                                      renewTechs, hoursForCESymbols)
    (genByNCTechAndZone, ncTechAndZoneToRow, hourToCol) = setupHourlyGenByTechAndLocs(ipmZoneNums, 'zones',
                                                                                      notCurtailedTechs,
                                                                                      hoursForCESymbols)
    return (genByCTechAndCell, genByRETechAndZone, genByNCTechAndZone, cTechAndCellToRow,
            reTechAndZoneToRow, ncTechAndZoneToRow, hourToCol)


# Setup results lists for techs indexed by cell (i.e., curtailed techs) or zone
# (i.e., renew or not curtailed techs).
# Output: 2d list of [[tech_loc,hour1,hour2,...],[tech_loc1,genvals...],...]
def setupHourlyGenByTechAndLocs(locsForTechs, cellOrZones, curtailedTechs, hoursForCESymbols):
    (techAndLocToRow, hourToCol) = (dict(), dict())
    # Create empty 2d list
    numRows = len(curtailedTechs) * len(locsForTechs) + 1  # +1 for header in new 2d list
    hourlyGenByTechAndLoc = []
    for idx in range(numRows): hourlyGenByTechAndLoc.append([''] * (1 + len(hoursForCESymbols)))
    # Add hours as first row, starting at col 1 since first col is gen IDs
    genIDLabel = 'genID_cell'
    hourlyGenByTechAndLoc[0] = [genIDLabel] + hoursForCESymbols
    # Create dict mapping hours to col #s
    for idx in range(1, len(hourlyGenByTechAndLoc[0])): hourToCol[hourlyGenByTechAndLoc[0][idx]] = idx
    # Add tech_cell as first col, starting at row 1 since first row is hours
    idx = 1
    for tech in curtailedTechs:
        if cellOrZones == 'cells':
            techAndCells = [createTechAndLocLabel(tech, loc) for loc in locsForTechs]
        elif cellOrZones == 'zones':
            techAndCells = [createTechAndLocLabel(tech, createZoneSymbol(loc)) for loc in locsForTechs]
        for techAndLoc in techAndCells:
            hourlyGenByTechAndLoc[idx][0] = techAndLoc
            techAndLocToRow[techAndLoc] = idx
            idx += 1
    assert (idx == len(hourlyGenByTechAndLoc))
    return (hourlyGenByTechAndLoc, techAndLocToRow, hourToCol)


############ SAVE GENERATOR RESULTS
# Save gen-level CE results
def saveCEResultsByPlantVar(genByPlant, genToRow, hourToCol, ceModel):
    saveHourByPlantVarCE(genByPlant, genToRow, hourToCol, ceModel, 'vPegu')


def saveCEResultsByTechVar(genByCTechAndCell, genByRETechAndZone, genByNCTechAndZone, cTechAndCellToRow,
                           reTechAndZoneToRow, ncTechAndZoneToRow, hourToCol, ceModel):

    saveLocByTechByHourVarCE(genByCTechAndCell, cTechAndCellToRow, hourToCol, ceModel, 'vPtechcurtailed')
    saveLocByTechByHourVarCE(genByRETechAndZone, reTechAndZoneToRow, hourToCol, ceModel, 'vPtechrenew')
    saveLocByTechByHourVarCE(genByNCTechAndZone, ncTechAndZoneToRow, hourToCol, ceModel, 'vPtechnotcurtailed')


# Extract results from CE GAMS model and add to proper list
def saveHourByPlantVarCE(varHourByPlantList, genToRow, hourToCol, ceModel, varName):
    for rec in ceModel.out_db[varName]:
        (rowIdx, colIdx) = (genToRow[rec.key(0)], hourToCol[rec.key(1)])  # vars indexed as egu,h or tech,h
        varHourByPlantList[rowIdx][colIdx] = rec.level


# Vars indexed as c/z,tech,h
def saveLocByTechByHourVarCE(result, genByTechAndLocToRow, hourToCol, ceModel, varName):
    for rec in ceModel.out_db[varName]:
        loc, tech, hr = rec.key(0), rec.key(1), rec.key(2)
        rowIdx = genByTechAndLocToRow[createTechAndLocLabel(tech, loc)]
        colIdx = hourToCol[hr]  # vars indexed as egu,h or tech,h
        result[rowIdx][colIdx] = rec.level
    ################################################################################


def savePumpedHydroResults(ceModel):
    """SAVE PUMPED HYDRO RESULTS

    :param ceModel:
    :param genFleetForCE:
    :param pumpedHydroSymbs:
    :param hoursForCESymbols:
    :return:
    """

    d = dict()
    for rec in ceModel.out_db['vCharge']:
        d[tuple(rec.keys)] = rec.get_level()

    chargeByPH = dict_tuples_to_list2d(d)

    # use pandas to sort resulting list
    chargeByPH = pd.DataFrame(chargeByPH, columns=['gcm', 'genId', 'hour', 'value'])
    chargeByPH = chargeByPH.sort_values(by=['gcm', 'hour', 'genId']).reset_index(drop=True)
    chargeByPH = chargeByPH.values.tolist()

    d = dict()
    for rec in ceModel.out_db['vSoc']:
        d[tuple(rec.keys)] = rec.get_level()

    socByPH = dict_tuples_to_list2d(d)

    # use pandas to sort resulting list
    socByPH = pd.DataFrame(socByPH, columns=['gcm', 'genId', 'hour', 'value'])
    socByPH = socByPH.sort_values(by=['gcm', 'hour', 'genId']).reset_index(drop=True)
    socByPH = socByPH.values.tolist()

    return chargeByPH, socByPH


def saveLineFlows(ceModel):
    """SAVE LINE FLOWS BETWEEN ZONES

    :param ceModel:
    :return:
    """
    d = dict()
    for rec in ceModel.out_db['vLineflow']:
        d[tuple(rec.keys)] = rec.get_level()

    flowByLine = dict_tuples_to_list2d(d)

    # use pandas to sort resulting list
    flowByLine = pd.DataFrame(flowByLine, columns=['gcm', 'line', 'hour', 'value'])
    flowByLine = flowByLine.sort_values(by=['gcm', 'hour', 'line']).reset_index(drop=True)
    flowByLine = flowByLine.values.tolist()


    return flowByLine


################################################################################

################### SAVE SYSTEM RESULTS ########################################
def saveSystemResults(ceModel, hoursForCESymbols, ipmZoneNums):
    resultLabels = ['mcGen']
    sysResults, resultToRow, hourToCol = setupHourlyZonalSysResults(hoursForCESymbols, resultLabels, ipmZoneNums)
    saveCEResultsByZonalSysVar(sysResults, resultToRow, hourToCol, ceModel)

    # Commented code is for hourly system results not indexed by zone
    # sysResults,resultToRow,hourToCol = setupHourlySystemResultsWithHourSymbols(hoursForCESymbols,resultLabels)
    # saveCEResultsBySysVar(sysResults,resultToRow,hourToCol,ceModel)

    return sysResults


########## SETUP EMPTY ZONAL LIST
def setupHourlyZonalSysResults(hourSymbols, resultLabels, ipmZoneNums):
    hourToColSys, resultToRow, sysResults, numRows = dict(), dict(), [], len(resultLabels) * len(ipmZoneNums) + 1
    # Initialize empty 2d list
    for idx in range(numRows): sysResults.append([''] * (1 + len(hourSymbols)))
    # Add hour labels across top
    sysResults[0] = ['Hour'] + hourSymbols
    for idx in range(1, len(sysResults[0])): hourToColSys[sysResults[0][idx]] = idx
    # Add row labels
    idx = 1
    for result in resultLabels:
        resultAndZones = [createResultAndZoneLabel(result, createZoneSymbol(zone)) for zone in ipmZoneNums]
        for resultAndZone in resultAndZones:
            sysResults[idx][0] = resultAndZone
            resultToRow[resultAndZone] = idx
            idx += 1
    return (sysResults, resultToRow, hourToColSys)


def createResultAndZoneLabel(result, zone):
    return result + '-' + zone


########### SAVE SYSTEM RESULTS
# Note that setupHourlySystemResultsWitHourSymbols also produces a row w/ nse,
# but nse is not in CE, so don't save those values here. That's also why
# need to do "for result in resultLabel..." rather than in resultToRow.
def saveCEResultsByZonalSysVar(sysResults, resultToRow, hourToColSys, ceModel):
    resultLabelToEqnName = {'mcGen': 'meetdemand'}
    for result in resultLabelToEqnName:
        varName = resultLabelToEqnName[result]
        for rec in ceModel.out_db[varName]:
            zone, hr = rec.key(0), rec.key(1)
            rowLabel = createResultAndZoneLabel(result, zone)
            (rowIdx, colIdx) = (resultToRow[rowLabel], hourToColSys[hr])
            if 'mc' in result:
                sysResults[rowIdx][colIdx] = rec.marginal
            else:
                sysResults[rowIdx][colIdx] = rec.level

        # #Note that setupHourlySystemResultsWitHourSymbols also produces a row w/ nse,


# #but nse is not in CE, so don't save those values here. That's also why
# #need to do "for result in resultLabel..." rather than in resultToRow.
# def saveCEResultsBySysVar(sysResults,resultToRow,hourToColSys,ceModel):
#     resultLabelToEqnName = {'mcGen':'meetdemand'}
#     for result in resultLabelToEqnName:
#         varName = resultLabelToEqnName[result]
#         for rec in ceModel.out_db[varName]:
#             (rowIdx,colIdx) = (resultToRow[result],hourToColSys[rec.key(0)])
#             if 'mc' in result: sysResults[rowIdx][colIdx] = rec.marginal
#             else: sysResults[rowIdx][colIdx] = rec.level 
################################################################################

################### SAVE ANNUAL COST AND CO2 EMISSIONS #########################
def saveCo2EmAndCostResults(ceModel):
    co2Ems = extract0dVarResultsFromGAMSModel(ceModel, 'vCO2emsannual')
    totalCost = extract0dVarResultsFromGAMSModel(ceModel, 'vZ')
    capCost = extract0dVarResultsFromGAMSModel(ceModel, 'vIc')
    varCost = extract0dVarResultsFromGAMSModel(ceModel, 'vVc')
    return [['co2ems tons', co2Ems], ['totalcost thousand$', totalCost],
            ['capcost thousand$', capCost], ['varcost thousand$', varCost]]
################################################################################
