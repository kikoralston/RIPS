# Michael Craig
# August 16, 2016

# Determine hourly capacity factors for candidate new plants in
# capacity expansion, depending on where they are assigned.

import os, copy
from AuxFuncs import *
from GAMSUtil.GAMSAuxFuncs import createTechSymbol
from AssignCellsToIPMZones import mapCellToIPMZone
from thermalderatings.AuxCurtailmentFuncs import loadWaterAndMetData, read_netcdf_full
from thermalderatings.CurtailmentRegressions import (calcCurtailmentForGenOrTech, loadRegCoeffs, getKeyCurtailParamsNewTechs,
                                                     getCoeffsForGenOrTech)
from AssignCellsToStates import getStateOfPt
import progressbar
import multiprocessing as mp


def determineHrlyCurtailmentsForNewTechs(eligibleCellWaterTs, newTechsCE, currYear, genparam, curtailparam, pbar=True):
    """ Main wrapper function to compute deratings of candidate thermal technologies for capacity expansion simulation.

    This function computes thermal deratings of candidate techs using parallel multiprocessing. It combines results from
    all different GCMs into a single nested dictionary.

    :param eligibleCellWaterTs:
    :param newTechsCE:
    :param currYear: (integer) current year of simulation
    :param genparam: object of class :mod:`Generalparameters`
    :param curtailparam: object of class :mod:`Curtailmentparameters`
    :param pbar: (boolean) If True show progress bar
    :return: nested dict of {gcm: {(plant+cooltype,cell): [hrly tech curtailments]}}
    """
    # Isolate water Ts to year of analysis
    eligibleCellWaterTsCurrYear = getWaterTsInCurrYear(currYear, eligibleCellWaterTs)
    cellWaterTsForNewTechs = selectCells(eligibleCellWaterTsCurrYear, genparam.cellNewTechCriteria,
                                         genparam.fipsToZones, genparam.fipsToPolys)

    # compile list of arguments to call multiprocessing
    args_list = [[cellWaterTsForNewTechs[gcm], newTechsCE, currYear, genparam, curtailparam, gcm]
                 for gcm in curtailparam.listgcms]

    ncores = min(len(curtailparam.listgcms), genparam.ncores_py)

    if ncores == 1:
        # do normal map instead of pool.map to avoid overhead
        list_curtailments = list(map(worker_new_tech_curtailments, args_list))
    else:
        # if ncores > 1, use multiprocessing
        with mp.Pool(processes=ncores) as pool:
            list_curtailments = pool.map(worker_new_tech_curtailments, args_list)

    # create nested dictionary
    hrlyCurtailmentsAllTechsInTgtYr = dict(zip(curtailparam.listgcms, list_curtailments))

    return hrlyCurtailmentsAllTechsInTgtYr


def worker_new_tech_curtailments(list_args):
    """ Worker function for using multiprocessing with function calculateGeneratorCurtailments.

    :param list_args: list with arguments for function
    :return: dict of {(plant+cooltype,cell): [hrly tech curtailments]} (see function calculateGeneratorCurtailments)
    """

    curtailments_out = calculateTechsCurtailments(cellWaterTsForNewTechs=list_args[0], newTechsCE=list_args[1],
                                                  currYear=list_args[2], genparam=list_args[3],
                                                  curtailparam=list_args[4], gcm=list_args[5])
    return curtailments_out


def calculateTechsCurtailments(cellWaterTsForNewTechs, newTechsCE, currYear, genparam, curtailparam, gcm, pbar=True):
    """ Computes deratings of candidate thermal technologies for capacity expansion simulation

    :param cellWaterTsForNewTechs:
    :param newTechsCE:
    :param currYear: (integer) current year of simulation
    :param genparam: object of class :mod:`Generalparameters`
    :param curtailparam: object of class :mod:`Curtailmentparameters`
    :param gcm: (string) name of GCM being simulated
    :param pbar: (boolean) show progress bar on standard output?
    :return: dict of {(plant+cooltype,cell): [hrly tech curtailments]}
    """

    # Do curtailments
    (hrlyCurtailmentsAllTechsInTgtYr, hrlyTechCurtailmentsList) = (dict(), [])
    regCoeffs = loadRegCoeffs(genparam.dataRoot, 'capacity.json')  # dict of cooling type: reg coeffs

    if genparam.incCurtailments or genparam.incRegs:
        # read full meteo data for current year (full water data is already loaded)
        fname = curtailparam.basenamemeteo
        name_gcm = gcm

        fname = os.path.join(curtailparam.rbmRootDir, fname.format(name_gcm, currYear))
        meteodata = read_netcdf_full(currYear, fname, curtailparam)
    else:
        meteodata = None

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

        if genparam.incCurtailments or genparam.incRegs:
            metAndWaterData = loadWaterAndMetData(currYear, cellLat, cellLong, genparam, curtailparam, gcm=gcm,
                                                  metdatatot=meteodata, waterDatatot=cellWaterTsForNewTechs)
        else:
            metAndWaterData = None

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
#                hrlyTechCurtailmentsList.append(
#                    [createTechSymbol(techRow, newTechsCE[0], genparam.ptCurtailedAll), cell]
#                    + list(hrlyCurtailmentsGen))

    return hrlyCurtailmentsAllTechsInTgtYr #, hrlyTechCurtailmentsList


def getWaterTsInCurrYear(currYear, eligibleCellWaterTs):
    """ISOLATE AVG WATER TS FOR CURRENT YEAR FOR EACH CELL

    :param currYear: (integer)
    :param eligibleCellWaterTs:
    :return: dict of cell folder name : [[Datetime],[AverageWaterT(degC)]]
    """
    eligibleCellWaterTsCurrYear = dict()

    for gcm in eligibleCellWaterTs:
        auxDict = {}
        for cell in eligibleCellWaterTs[gcm]:
            cellWaterTs = eligibleCellWaterTs[gcm][cell]
            if cellWaterTs is not None:
                dateCol = cellWaterTs[0].index('date')
                auxDict[cell] = [cellWaterTs[0]] + [row for row in cellWaterTs[1:] if str(currYear) in str(row[dateCol])]
            else:
                auxDict[cell] = None

        eligibleCellWaterTsCurrYear[gcm] = auxDict

    return eligibleCellWaterTsCurrYear


def selectCells(eligibleCellWaterTsCurrYear, cellNewTechCriteria, fipsToZones, fipsToPolys):
    """Filter list of grid cells that can host new techs

    This function filters the list of eligible cells and chooses only the ones with maximum water temperature in each zone
    if the criteria is 'maxWaterT'

    The parameter ``cellNewTechCriteria` is defined in the :mod:`Generalparameters` object. It accepts two options:

    * 'all': no filtering is applied. All eligible cells are considered
    * 'maxWaterT': only the cell with maximum water temeprature in each zone is considered

    :param eligibleCellWaterTsCurrYear:
    :param cellNewTechCriteria:
    :param fipsToZones:
    :param fipsToPolys:
    :return:
    """
    if cellNewTechCriteria == 'all':
        cellWaterTsForNewTechs = copy.deepcopy(eligibleCellWaterTsCurrYear)

    elif cellNewTechCriteria == 'maxWaterT':
        cellsPerZone = getCellsInEachZone(eligibleCellWaterTsCurrYear, fipsToZones, fipsToPolys)
        cellWaterTsForNewTechs = placeTechInMaxWaterTCell(eligibleCellWaterTsCurrYear, cellsPerZone)
    else:
        print('Parameter \'cellWaterTsForNewTechs\' must be either \'all\' or \'maxWaterT\'.... Ending simulation!')
        sys.exit()

    return cellWaterTsForNewTechs


def getCellsInEachZone(eligibleCellWaterTsCurrYear, fipsToZones, fipsToPolys):
    """Takes in dict of cell:waterTs, and returns dict of zone:all cells in that zone

    :param eligibleCellWaterTsCurrYear:
    :param fipsToZones:
    :param fipsToPolys:
    :return:
    """
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

