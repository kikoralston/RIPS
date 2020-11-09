# Michael Craig
# March 10, 2017
# Determine which coal plants should retire based on CF in PRIOR CE run. (May
# not retire after prior CE run to maintain planning margin.)

from GAMSUtil.GAMSAuxFuncs import *
from ProcessResults.ProcessCEResults import markRetiredUnitsFromCE
import pandas as pd
import numpy as np


def retireUnitsCFPriorCE(genFleet, genFleetPriorCE, retirementCFCutoff, priorCEout_db, priorHoursCE, scaleMWtoGW,
                         ptEligRetireCF, currYear):

    unitsRetireCF = markAndSaveRetiredUnitsFromPriorCE(retirementCFCutoff, genFleet, genFleetPriorCE,
                                                       priorCEout_db, priorHoursCE, scaleMWtoGW, ptEligRetireCF,
                                                       currYear)
    return unitsRetireCF


def markAndSaveRetiredUnitsFromPriorCE(retirementCFCutoff, genFleet, genFleetPriorCE, priorCEout_db,
                                       priorHoursCE, scaleMWtoGW, ptEligRetireCF, currYear):
    """RETIRE UNITS BASED ON CF IN PRIOR CE RUN

    Retires coal plants that don't meet CF threshold in PRIOR CE run. These units may not retire due to limiting
    retirements per planning margin.

    :param retirementCFCutoff: Retirement CF cutoff
    :param genFleet: generator fleet
    :param genFleetPriorCE:
    :param priorCEout_db: GAMS CE model output db from prior CE run
    :param priorHoursCE: hours included in CE
    :param scaleMWtoGW: scale MW to GW
    :param ptEligRetireCF: 1d list of plant types eligible to retire based on CFs
    :param currYear:
    :return: 1d list w/ units that retire due to CF
    """
    unitsRetireCF = getUnitsRetireByCF(retirementCFCutoff, genFleet, genFleetPriorCE, priorCEout_db,
                                       priorHoursCE, scaleMWtoGW, ptEligRetireCF)
    markRetiredUnitsFromCE(genFleet, unitsRetireCF, currYear)

    return unitsRetireCF


def getUnitsRetireByCF(retirementCFCutoff, genFleet, genFleetPriorCE, priorCEout_db, priorHoursCE,
                       scaleMWtoGW, ptEligRetireCF):
    """Determines which units retire due to CF in prior CE run

# Inputs: , dictionary (genID:CF) for generators eligible to retire based on CF,
# planning reserve, curr & end year, gen fleet w/ only online gens, plant types that
# can retire due to CF in CE.
# Outputs: 1d list of units to retire for economic reasons


    :param retirementCFCutoff: CF cutoff retirement
    :param genFleet:
    :param genFleetPriorCE:
    :param priorCEout_db:
    :param priorHoursCE:
    :param scaleMWtoGW:
    :param ptEligRetireCF:
    :return:
    """
    # (gcm,genID,hr):hourly gen [GW]
    hourlyGenByGens = dict()
    for rec in priorCEout_db['vPegu']:
        hourlyGenByGens[tuple(rec.keys)] = float(rec.level) * scaleMWtoGW

    # convert to list
    hourlyGenByGens = list(hourlyGenByGens.items())

    # convert to pd data frame
    hourlyGenByGensDf = pd.DataFrame({'gcm': [k[0][0] for k in hourlyGenByGens],
                                      'genId': [k[0][1] for k in hourlyGenByGens],
                                      'genValue': [i[1] for i in hourlyGenByGens]})

    # aggregate average annual generation for each plant (over all gcms, over all hours)
    hourlyGenByGensDf = hourlyGenByGensDf.groupby(['gcm', 'genId'], as_index=False).aggregate(np.sum)
    hourlyGenByGensDf = hourlyGenByGensDf.groupby('genId', as_index=False).aggregate(np.mean)

    # convert to dictionary (original format)
    keys = hourlyGenByGensDf['genId'].values
    values = hourlyGenByGensDf['genValue'].values
    ceHoursGenByGens = dict(zip(keys, values))

    # get longest hoursForCE (different gcms could have different lengths because of special hours)
    auxArray = np.array([[i[0], float(len(i[1]))] for i in priorHoursCE.items()], dtype=object)
    idxMax = np.argmax(auxArray[:, 1])
    longGcm = auxArray[idxMax, 0]

    priorHoursCEAux = priorHoursCE[longGcm]

    #hourlyGenByGens = extract2dVarResultsIntoDict(priorCapacExpModel, 'vPegu')  # (genID,hr):hourly gen [GW]
    #ceHoursGenByGens = sumHourlyGenByGensInCE(hourlyGenByGens, scaleMWtoGW)

    gensEligToRetireCFs = getGenCFsInCENotAlreadyRetired(ceHoursGenByGens, genFleet,
                                                         genFleetPriorCE, ptEligRetireCF, priorHoursCEAux)
    unitsRetireCF = []
    if len(gensEligToRetireCFs) > 0:
        minCF = min([gensEligToRetireCFs[gen] for gen in gensEligToRetireCFs])
        if minCF < retirementCFCutoff:  # if any plants eligible for retirement
            addAllUnitsWithCFBelowCutoff(gensEligToRetireCFs, retirementCFCutoff, unitsRetireCF)

    return unitsRetireCF


def getGenCFsInCENotAlreadyRetired(ceHoursGenByGens, genFleet, genFleetPriorCE, plantTypesEligibleForRetirementByCF,
                                   hoursForCE):
    """Determines which gens retire due to CF in prior CE run that didn't already retire after last CE run.

# Inputs: total gen by generators for CE hours (dict of genID:total gen), gen fleet
# only w/ online generators (2d list), list of plant types that can retire based on CF,
# 1d list of hours input to CE
# Outputs: dictionary (genID:CF) for generators eligible to retire based on CF

    :param ceHoursGenByGens:
    :param genFleet:
    :param genFleetPriorCE:
    :param plantTypesEligibleForRetirementByCF:
    :param hoursForCE:
    :return:
    """
    gensEligToRetireCFs = dict()
    (capacCol, plantTypeCol) = (genFleet[0].index('Capacity (MW)'),
                                genFleet[0].index('PlantType'))
    econRetCol = genFleet[0].index('YearRetiredByCE')
    genSymbolsFleetFull = [createGenSymbol(row, genFleet[0]) for row in genFleet]
    genSymbolsFleetPriorCE = [createGenSymbol(row, genFleetPriorCE[0]) for row in genFleetPriorCE]
    for gen in ceHoursGenByGens:
        # Need to screen out wind and solar plants in genFleetPriorCE, since these
        # were tacked on @ end of fleet and are not in genFleet. Consequently, if don't
        # have this if statement and don't build new plants, genSymbolsFleetFull.index(gen)
        # call will not find gen listed.
        if (genFleetPriorCE[genSymbolsFleetPriorCE.index(gen)][plantTypeCol] != 'Wind' and
                genFleetPriorCE[genSymbolsFleetPriorCE.index(gen)][plantTypeCol] != 'Solar PV'):
            genRow = genSymbolsFleetFull.index(gen)
            if genFleet[genRow][plantTypeCol] in plantTypesEligibleForRetirementByCF:
                if genFleet[genRow][econRetCol] == '':
                    genCapac = genFleet[genRow][capacCol]
                    genCF = ceHoursGenByGens[gen] / (float(genCapac) * len(hoursForCE))
                    gensEligToRetireCFs[gen] = genCF
    return gensEligToRetireCFs


def addAllUnitsWithCFBelowCutoff(gensEligToRetireCFs, retirementCFCutoff, unitsRetireCF):
    """Adds all units elig to retire w/ CF below cutoff to unitsToRetire list

# Inputs: dictionary (genID:CF) for gens elig to retire, retirement CF cutoff,
# empty list to which genIDs for units that should retire are added

    :param gensEligToRetireCFs:
    :param retirementCFCutoff:
    :param unitsRetireCF:
    """
    for gen in gensEligToRetireCFs:
        genCF = gensEligToRetireCFs[gen]
        if genCF < retirementCFCutoff: unitsRetireCF.append(gen)