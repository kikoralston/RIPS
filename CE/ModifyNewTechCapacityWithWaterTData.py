# Michael Craig
# August 16, 2016

# Determine hourly capacity factors for candidate new plants in
# capacity expansion, depending on where they are assigned.

import os, csv, copy
from AuxFuncs import *
from GAMSAuxFuncs import createTechSymbol
from AssignCellsToIPMZones import mapCellToIPMZone
from ModifyGeneratorCapacityWithWaterTData import loadWaterAndMetData, read_netcdf_full
from CurtailmentRegressions import (calcCurtailmentForGenOrTech, loadRegCoeffs, getKeyCurtailParamsNewTechs,
                                    getCoeffsForGenOrTech)
from AssignCellsToStates import getStateOfPt
import numpy as np
import pandas as pd
import progressbar


################################################################################
####### MASTER FUNCTION ########################################################
################################################################################
# Returns dict of (plant+cooltype,cell):hrly tech curtailments
def determineHrlyCurtailmentsForNewTechs(eligibleCellWaterTs, newTechsCE, currYear, genparam, curtailparam,
                                         resultsDir, pbar=True):

    # Isolate water Ts to year of analysis
    eligibleCellWaterTsCurrYear = getWaterTsInCurrYear(currYear, eligibleCellWaterTs)
    cellWaterTsForNewTechs = selectCells(eligibleCellWaterTsCurrYear, genparam.cellNewTechCriteria,
                                         genparam.fipsToZones, genparam.fipsToPolys)
    # Do curtailments
    (hrlyCurtailmentsAllTechsInTgtYr, hrlyTechCurtailmentsList) = (dict(), [])
    regCoeffs = loadRegCoeffs(genparam.dataRoot, 'capacity.json')  # dict of cooling type: reg coeffs

    # read full meteo data for current year (full water data is already loaded)
    fname = curtailparam.basenamemeteo

    # TODO: find a way to combine different GCMs
    gcm = curtailparam.listgcms[0]

    fname = os.path.join(curtailparam.rbmDataDir, fname.format(gcm, currYear))
    meteodata = read_netcdf_full(currYear, fname, curtailparam)

    if pbar:
        ProgressBar = progressbar.ProgressBar
    else:
        ProgressBar = progressbar.NullBar

    bar = ProgressBar()

    for cell in bar(cellWaterTsForNewTechs):
        t0 = time.time()
        #print(cell)

        cellLat, cellLong = float(cell.split('_')[0]), float(cell.split('_')[1])
        state = getStateOfPt(genparam.statePolys, cellLat, cellLong)  # need for regression

        #print('  Got state. ' + str_elapsedtime(t0))

        metAndWaterData = loadWaterAndMetData(currYear, cellLat, cellLong, genparam, curtailparam, metdatatot=meteodata,
                                              waterDatatot=cellWaterTsForNewTechs)

        #print('  Got met data. ' + str_elapsedtime(t0))

        for techRow in newTechsCE[1:]:
            (plantType, hr, fuelAndCoalType, coolType,
             fgdType, cap, coolDesignT) = getKeyCurtailParamsNewTechs(newTechsCE, techRow)

            #print('  Starting computation for plant ' + plantType)

            coeffs = getCoeffsForGenOrTech(plantType, coolType, genparam.ptCurtailed, regCoeffs, coolDesignT)

            if (coeffs is not None) and (plantType in genparam.ptCurtailedRegs):  # skip gens that aren't curtailed

                hrlyCurtailmentsGen = calcCurtailmentForGenOrTech(plantType, fuelAndCoalType, coolType, state, cap,
                                                                  coolDesignT, metAndWaterData, coeffs, genparam,
                                                                  curtailparam)

                #print('  computed curtailment. ' + str_elapsedtime(t0))
                hrlyCurtailmentsAllTechsInTgtYr[
                    (createTechSymbol(techRow, newTechsCE[0], genparam.ptCurtailedAll),
                     cell)] = hrlyCurtailmentsGen
                hrlyTechCurtailmentsList.append(
                    [createTechSymbol(techRow, newTechsCE[0], genparam.ptCurtailedAll), cell]
                    + list(hrlyCurtailmentsGen))

    return hrlyCurtailmentsAllTechsInTgtYr, hrlyTechCurtailmentsList


################################################################################
####### ISOLATE AVG WATER TS FOR CURRENT YEAR FOR EACH CELL ####################
################################################################################
# Returns dict of cell folder name : [[Datetime],[AverageWaterT(degC)]]
def getWaterTsInCurrYear(currYear, eligibleCellWaterTs):

    eligibleCellWaterTsCurrYear = dict()
    for cell in eligibleCellWaterTs:
        cellWaterTs = eligibleCellWaterTs[cell]
        dateCol = cellWaterTs[0].index('date')
        eligibleCellWaterTsCurrYear[cell] = [cellWaterTs[0]] + [row for row in
                                                                cellWaterTs[1:] if str(currYear) in str(row[dateCol])]
    return eligibleCellWaterTsCurrYear


################################################################################
################################################################################
################################################################################

################################################################################
####### PICK GRID CELL TO PUT NEW TECH IN ######################################
################################################################################
def selectCells(eligibleCellWaterTsCurrYear, cellNewTechCriteria, fipsToZones, fipsToPolys):
    if cellNewTechCriteria == 'all':
        cellWaterTsForNewTechs = copy.deepcopy(eligibleCellWaterTsCurrYear)
    elif cellNewTechCriteria == 'maxWaterT':
        cellsPerZone = getCellsInEachZone(eligibleCellWaterTsCurrYear, fipsToZones, fipsToPolys)
        cellWaterTsForNewTechs = placeTechInMaxWaterTCell(eligibleCellWaterTsCurrYear, cellsPerZone)
    else:
        print('Parameter \'cellWaterTsForNewTechs\' must be either \'all\' or \'maxWaterT\'.... Ending simulation!')
        sys.exit()

    return cellWaterTsForNewTechs


# Takes in dict of cell:waterTs, and returns dict of zone:all cells in that zone
def getCellsInEachZone(eligibleCellWaterTsCurrYear, fipsToZones, fipsToPolys):
    cellsPerZone = dict()
    for cell in eligibleCellWaterTsCurrYear:
        cellZone = mapCellToIPMZone(cell, fipsToZones, fipsToPolys)
        if cellZone in cellsPerZone:
            cellsPerZone[cellZone].append(cell)
        else:
            cellsPerZone[cellZone] = [cell]
    return cellsPerZone


# Choose cell w/ highest peak water T in each zone and return in dict of
# cell:water Ts
def placeTechInMaxWaterTCell(eligibleCellWaterTsCurrYear, cellsPerZone):
    maxPeakWaterTCellWaterTs = dict()
    for zone in cellsPerZone:  # for each zone, get cells
        (maxPeakWaterT, maxPeakWaterTCell) = (None, None)
        cellsInZone = cellsPerZone[zone]
        for cell in cellsInZone:  # across cells in zone, get cell w/ max water T
            # already know cell is in eligible... b/c that's where they came from (see getCellsInEachZone)
            cellWaterTData = eligibleCellWaterTsCurrYear[cell]
            waterTCol = cellWaterTData[0].index('AverageWaterT(degC)')
            cellWaterTs = [float(row[waterTCol]) for row in cellWaterTData[1:]]
            peakWaterT = max(cellWaterTs)
            if maxPeakWaterT == None or peakWaterT > maxPeakWaterT:
                (maxPeakWaterT, maxPeakWaterTCell) = (peakWaterT, cell)
        maxPeakWaterTCellWaterTs[maxPeakWaterTCell] = copy.deepcopy(eligibleCellWaterTsCurrYear[maxPeakWaterTCell])
    return maxPeakWaterTCellWaterTs


################################################################################
################################################################################
################################################################################

################################################################################
####### LOAD MET DATA FOR CELL AND SLIM WATER T DATA TO YEAR ###################
################################################################################
def loadMetAndAddWaterData(cellLat, cellLong, rbmOutputDir, locPrecision, currYear,
                           cellWaterTAndDates, resultsDir):
    metData = loadMetData(cellLat, cellLong, rbmOutputDir, locPrecision, currYear)
    cellWaterTs = [float(row[cellWaterTAndDates[0].index('AverageWaterT(degC)')])
                   for row in cellWaterTAndDates[1:]]
    hourlyWaterT = list(np.array([[val] * 24 for val in cellWaterTs]).flatten())
    assert (len(hourlyWaterT) == metData.shape[0])
    metData['waterC'] = pd.Series(hourlyWaterT)
    metData['waterF'] = metData.loc[:, 'waterC'] * 9 / 5 + 32
    metData['airF'] = metData.loc[:, 'tC'] * 9 / 5 + 32
    metData.to_csv(os.path.join(resultsDir, 'metAndWaterTech' + str(cellLat) + '_' + str(cellLong) + '.csv'))
    return metData
################################################################################
################################################################################
################################################################################
