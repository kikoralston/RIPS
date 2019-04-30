# Michael Craig, 31 May 2016
# Run capacity expansion and dispatch models w/ climate impacts.

import os, sys
import csv, operator, copy, time, random
import numpy as np
import datetime as dt
import multiprocessing as mp
import shutil
import gc
import platform
import pickle as pk

try:
    from gams import *
except ImportError:
    print('gams module not found! GAMS functions will not work.')

from AuxFuncs import *
from GAMSAuxFuncs import *
from SetupGeneratorFleet import setupGeneratorFleet, aggregatePlantTypeToORIS
# from ImportDemandProfile import importZonalHourlyDemand
from ForecastDemandWithRegression import forecastZonalDemandWithReg
from UpdateFuelPriceFuncs import *
from DemandFuncs import *
from DemandFuncsCE import *
from GetHydroMaxGenPotential import getHydroEPotential, importHydroPotentialGen, getDailyHydroPotentialsUC
from CO2CapCalculations import getCo2Cap, interpolateCO2Cap
from SetInitCondsUC import *
from ImportNewTechs import getNewTechs
from RetireUnitsCFPriorCE import retireUnitsCFPriorCE
from CreateFleetForCELoop import createFleetForCurrentCELoop, onlineAndNotRetired
from GetRenewableCFs import getRenewableCFs
from GetNewRenewableCFs import getNewWindAndSolarCFs, trimNewRECFsToCEHours
from AssignCellsToIPMZones import getIPMPolys, assignCellsToIPMZones
from AssignCellsToStates import getStatePolys
from LoadEligibleCellWaterTs import loadEligibleCellWaterTs
from GetHourlyCapacsForCE import getHourlyNonRECapacsForCE, getHourlyCurtailedTechCapacsForCE
from GetHourlyCapacsForUC import getHourlyCapacitiesForDays
from CalculateHourlyCapacsWithCurtailments import calculateHourlyCapacsWithCurtailments
from ProcessCEResults import *
from CombineWindAndSolarGensToSingleGen import combineWindAndSolarToSinglePlant
from TrimDemandREGenAndResForUC import getDemandAndREGenForUC, getResForUC
from GAMSAddSetToDatabaseFuncs import *
from GAMSAddParamToDatabaseFuncs import *
from ConvertCO2CapToPrice import convertCo2CapToPrice
from SetupResultLists import setupHourlyResultsByPlant, setupHourlySystemResults, setupHourlyPHResults
from SaveHourlyResults import saveHourlyResultsByPlant, saveHourlySystemResults, saveHourlyPumpedHydroResults
from WriteUCResults import writeHourlyResultsByPlant, writeHourlyStoResults
from ReservesWWSIS import calcWWSISReserves
from SaveCEOperationalResults import saveCapacExpOperationalData
from LoadCEFleet import loadCEFleet
from ModifyGeneratorCapacityWithWaterTData import determineHrlyCurtailmentsForExistingGens
from PreProcessRBM import processRBMDataIntoIndividualCellFiles
from ModifyNewTechCapacityWithWaterTData import determineHrlyCurtailmentsForNewTechs

from Parameters import *

sys.stdout.flush()


def masterFunction(genparam, reserveparam, curtailparam):
    """MASTER FUNCTION

    """
    # initialize random number generator with fixed value (for 'reproducibility')
    random.seed(a=27081979)

    if not os.path.exists(genparam.resultsDir):
        os.makedirs(genparam.resultsDir)

    # create new folder to store GAMS models for this instance
    if not os.path.exists(os.path.join(genparam.resultsDir, 'GAMS')):
        os.makedirs(os.path.join(genparam.resultsDir, 'GAMS'))

    # copy GAMS files to GAMS dir in results folder
    try:
        shutil.copy(os.path.join(genparam.dataRoot, 'GAMS', genparam.capacExpFilename),
                    os.path.join(genparam.resultsDir, 'GAMS', genparam.capacExpFilename))
        shutil.copy(os.path.join(genparam.dataRoot, 'GAMS', genparam.ucFilename),
                    os.path.join(genparam.resultsDir, 'GAMS', genparam.ucFilename))
        # eg. src and dest are the same file
    except shutil.Error as e:
        print('Error: {}'.format(e))
        # eg. source or destination doesn't exist
    except IOError as e:
        print('Error: {}'.format(e.strerror))

    write2dListToCSV([['Zone', 'ZoneNum']] + rotate([genparam.ipmZones, genparam.ipmZoneNums]),
                     os.path.join(genparam.resultsDir, 'zoneNamesToNumbers.csv'))

    create_description_file(genparam, curtailparam)

    print()

    genFleet = getInitialFleetAndDemand(genparam, reserveparam)

    if genparam.processRBMData:
        print('Processing raw RBM data')
        processRBMDataIntoIndividualCellFiles(curtailparam=curtailparam)
        print('Processed RBM data')

    (capacExpModelsEachYear, capacExpBuilds, capacExpGenByGens, capacExpRetiredUnitsByCE,
     capacExpRetiredUnitsByAge) = ([], [['TechnologyType']], [['ORIS+UnitID']], [], [])

    t_start = time.time()
    # begin loop
    for currYear in range(genparam.startYear, genparam.endYear, genparam.yearStepCE):

        if not genparam.referenceCase:

            if genparam.gcmranking is not None:

                file_demand = 'df_demand_{1:}_{0:4d}.pk'.format(currYear, genparam.rcp)
                file_curtailment = 'curtailments_{1}_{0:4d}.pk'.format(currYear, genparam.rcp)

                # read demand file
                with open(os.path.join(curtailparam.rbmDataDir, file_demand), 'rb') as f:
                    df_demand = pk.load(f)

                # compute total hourly system demand
                df_demand = df_demand.groupby(['gcm', 'hour']).agg({'demand': sum}).reset_index()

                df_demand = df_demand.groupby(['gcm']).agg({'demand': max}).reset_index().sort_values(
                    by=['demand']).reset_index(drop=True)

                # get names of GCms in according to given ranking in genparam
                gcms_chosen = list(df_demand.iloc[genparam.gcmranking,]['gcm'].astype('str'))
            else:
                # if no ranking was given, just use list of gcms in curtailparam
                gcms_chosen = list(curtailparam.listgcms)

            print('GCMs used in year {}:'.format(currYear))
            print(gcms_chosen)
            print()

            # make copy of curtail param and update list of GCMs
            curtparam_year = copy.deepcopy(curtailparam)
            curtparam_year.listgcms = gcms_chosen
        else:
            curtparam_year = copy.deepcopy(curtailparam)

        t_year = time.time()
        print('\n----------------------------------------------------------------\n')
        print('Starting loop for year {0:4d}\n'.format(currYear))
        currCo2Cap = interpolateCO2Cap(currYear, genparam)

        print('CO2 cap in year {1:4d}: {0:,.3f} million tons'.format(currCo2Cap / 1e6, currYear))

        zonalDemandProfile, zonalTempDfs = forecastZonalDemandWithReg(currYear, genparam, curtparam_year)

        if genparam.runCE:
            # run capacity expansion model

            if genparam.coldStart and currYear == genparam.startYear:
                # if it is cold start and curr year is start year, read results from previous run
                # (this was created to treat cases when simulation crashes in the middle)

                priorYearCE = genparam.startYear - genparam.yearStepCE

                genFleet = readCSVto2dList(os.path.join(genparam.resultsDir, 'CE',
                                                        'genFleetAfterCE{0}.csv'.format(priorYearCE)))
                genFleetNoRetiredUnits = readCSVto2dList(os.path.join(genparam.resultsDir, 'CEtoUC',
                                                                      'genFleetCEtoUC{0}.csv'.format(priorYearCE)))
                genFleetPriorCE = readCSVto2dList(os.path.join(genparam.resultsDir, 'CE',
                                                               'genFleetForCE{0}.csv'.format(priorYearCE)))

                priorCEout_db = GamsWorkspace().add_database_from_gdx(os.path.join(genparam.resultsDir, 'CE',
                                                                                   'gdxOutYear{}.gdx'.format(
                                                                                       priorYearCE)))

                capacExpBuilds = readCSVto2dList(os.path.join(genparam.resultsDir, 'CE',
                                                              'genAdditionsCE{0}.csv'.format(priorYearCE)))
                capacExpGenByGens = readCSVto2dList(os.path.join(genparam.resultsDir, 'CE',
                                                                 'genByGensCE{0}.csv'.format(priorYearCE)))
                capacExpRetiredUnitsByCE = readCSVto2dList(os.path.join(genparam.resultsDir, 'CE',
                                                                        'genRetirementsEconCE{0}.csv'.format(
                                                                            priorYearCE)))
                capacExpRetiredUnitsByAge = readCSVto2dList(os.path.join(genparam.resultsDir, 'CE',
                                                                         'genRetirementsAgeCE{0}.csv'.format(
                                                                             priorYearCE)))

                with open(os.path.join(genparam.resultsDir, 'CE', 'hoursCE_{0}.pkl'.format(priorYearCE)), 'rb') as f:
                    priorHoursCE = pk.load(f)

            elif (not genparam.coldStart) and (currYear == genparam.startYear + genparam.yearStepCE):
                # not cold start and first CE run
                priorCEout_db, priorHoursCE, genFleetPriorCE = None, None, None  # first CE run

            # only run model in genparam.startYear if it is a cold start
            if genparam.coldStart or currYear > genparam.startYear:
                (genFleet, genFleetNoRetiredUnits, genFleetPriorCE, priorCEout_db,
                 priorHoursCE) = runCapacityExpansion(genFleet, zonalDemandProfile, currYear, currCo2Cap,
                                                      capacExpModelsEachYear, capacExpBuilds, capacExpGenByGens,
                                                      capacExpRetiredUnitsByCE, capacExpRetiredUnitsByAge,
                                                      genFleetPriorCE, priorCEout_db, priorHoursCE,
                                                      genparam, reserveparam, curtparam_year)

                # write hours for CE to pickle file
                with open(os.path.join(genparam.resultsDir, 'CE', 'hoursCE_{0}.pkl'.format(currYear)), 'wb') as f:
                    pk.dump(priorHoursCE, f, pk.HIGHEST_PROTOCOL)

            # since model was executed successfully, set coldStart to False before running next year
            genparam.coldStart = False

        if genparam.runUC:
            # Either only runs 2015, or runs in all but 2015
            if ((currYear == genparam.startYear and genparam.runFirstUCYear) or
                    (currYear > genparam.startYear and not genparam.runFirstUCYear)):
                if currYear == genparam.startYear:
                    genFleetNoRetiredUnits = genFleet
                else:
                    genFleetNoRetiredUnits = loadCEFleet(currYear, genparam.resultsDir)

                (ucResultsByDay, hourlyGenerationByPlants) = runUnitCommitment(genFleetNoRetiredUnits,
                                                                               zonalDemandProfile, currYear, currCo2Cap,
                                                                               genparam, reserveparam, curtailparam)

        print()
        print('Elapsed Time: ' + str_elapsedtime(t_year))

        print('\n----------------------------------------------------------------')

    print()
    print('************ Total Elapsed Time: ' + str_elapsedtime(t_start) + '************')


def getInitialFleetAndDemand(genparam, reserveparam):
    """SET UP INITIAL FLEET

    Reads folders containing data of existing power plant data and compiles 2d list with generator fleet data. If the
    genparam.analysisArea == 'test', it reads a file 'genFleetInitial.csv' from genparam.dataRoot.

    :param genparam: object of class Generalparameters
    :param reserveparam: object of class Reserveparameters
    :return: 2d list with initial generator fleet data
    """

    if genparam.analysisArea == 'test':
        fileNameWithDir = os.path.join(genparam.dataRoot, 'genFleetInitial.csv')

        if os.path.isfile(fileNameWithDir):
            genFleet = readCSVto2dList(fileNameWithDir=fileNameWithDir)
        else:
            print('Error! analysisArea is set to \'test\'. There must be a file \'genFleetInitial.csv\' at the ' +
                  'root of the data folder')
            sys.exit()
    else:
        genFleet = setupGeneratorFleet(genparam.startYear, genparam, reserveparam)

        aggregatePlantTypeToORIS(genFleet, 'Hydro')
        aggregatePlantTypeToORIS(genFleet, 'Pumped Storage')
        # Add placeholder columns for additions & retirements by CE
        ceHeaders = ['YearAddedCE', 'YearRetiredByCE', 'YearRetiredByAge']
        genFleet[0].extend(ceHeaders)
        for row in genFleet[1:]: row.extend([''] * len(ceHeaders))

    write2dListToCSV(genFleet, os.path.join(genparam.resultsDir, 'genFleetInitial.csv'))

    return genFleet


def runCapacityExpansion(genFleet, zonalDemandProfile, currYear, currCo2Cap, capacExpModelsEachYear, capacExpBuilds,
                         capacExpGenByGens, capacExpRetiredUnitsByCE, capacExpRetiredUnitsByAge, genFleetPriorCE,
                         priorCEout_db, priorHoursCE, genparam, reserveparam, curtailparam):
    """RUN CAPACITY EXPANSION

    This function does all the reading and preprocessing before executing optimization

    :param genFleet: 2d list with original existing generator fleet (in the first year of the simulation)
    :param zonalDemandProfile: nested dictionary with hourly demand by ipm zone in each gcm case being simulated.
                               Dictionary is {gcm: {zone: [hourly demand]}}
    :param currYear:
    :param currCo2Cap:
    :param capacExpModelsEachYear:
    :param capacExpBuilds:
    :param capacExpGenByGens:
    :param capacExpRetiredUnitsByCE:
    :param capacExpRetiredUnitsByAge:
    :param genFleetPriorCE: 2d list with existing generator fleet in current year of simulation
                            (including prior additions)
    :param priorCEout_db:
    :param priorHoursCE:
    :param genparam: object of class Generalparameters
    :param reserveparam: object of class Reserveparameters
    :param curtailparam: object of class Curtailmentparameters
    :return:
    """

    resultsDir = os.path.join(genparam.resultsDir, 'CE')

    if not os.path.exists(resultsDir): os.makedirs(resultsDir)
    print('Entering CE loop for year {0:4d}'.format(currYear))

    write2dListToCSV([[currCo2Cap]], os.path.join(resultsDir, 'co2CapCE' + str(currYear) + '.csv'))

    newTechsCE = getNewTechs(currYear, genparam, reserveparam)

    updateFuelPrices(genFleet, newTechsCE, currYear, genparam.fuelPricesTimeSeries)

    write2dListToCSV(newTechsCE, os.path.join(resultsDir, 'newTechsCE' + str(currYear) + '.csv'))

    if priorCEout_db is not None:  # if not in first CE loop
        unitsRetireCFPriorCE = retireUnitsCFPriorCE(genFleet, genFleetPriorCE, genparam.retirementCFCutoff,
                                                    priorCEout_db, priorHoursCE, genparam.scaleMWtoGW,
                                                    genparam.ptEligRetCF, currYear)
        print('Units that retire due to econ from prior CE ' + str(currYear) + ':', unitsRetireCFPriorCE)
        write2dListToCSV([unitsRetireCFPriorCE],
                         os.path.join(resultsDir, 'genRetirementsEconCEPrior' + str(currYear) + '.csv'))

    genFleetForCE = createFleetForCurrentCELoop(genFleet, currYear, capacExpRetiredUnitsByAge,
                                                genparam.dataRoot, genparam.scenario)  # removes all retired units

    print('Num units that retire due to age in ' + str(currYear) + ':' + str(len(capacExpRetiredUnitsByAge[-1]) - 1))

    write2dListToCSV(genFleetForCE, os.path.join(resultsDir, 'genFleetForCEPreRECombine' + str(currYear) + '.csv'))
    combineWindAndSolarToSinglePlant(genFleetForCE, genparam.ipmZones, genparam.dataRoot)  # combines by zone
    write2dListToCSV(genFleetForCE, os.path.join(resultsDir, 'genFleetForCE' + str(currYear) + '.csv'))

    writeDictToCSV(zonalDemandProfile, os.path.join(resultsDir, 'demandFullYrZonalCE' + str(currYear) + '.csv'))

    (startWindCapacForCFs, startSolarCapacForCFs) = (0, 0)
    zoneCol = genFleetForCE[0].index('Region Name')
    zonalHourlyWindGen, zonalHourlySolarGen = dict(), dict()
    zonalNewWindCFs, zonalNewSolarCFs = dict(), dict()

    print('Computing CFs for Renewables...')
    for zone in genparam.ipmZones:
        print('Zone ', zone)

        start_time = time.time()
        # Check 5/18/17: next 2 lines are working properly!
        zonalGenFleet = [genFleetForCE[0]] + [row for row in genFleetForCE if row[zoneCol] == zone]

        (windCFs, windCfsDtHr, windCfsDtSubhr, ewdIdAndCapac, solarCFs, solarCfsDtHr, solarCfsDtSubhr,
         solarFilenameAndCapac) = getRenewableCFs(zonalGenFleet, startWindCapacForCFs, startSolarCapacForCFs,
                                                  genparam.tzAnalysis, genparam.dataRoot, genparam.windGenDataYr,
                                                  zone, genparam.fipsToZones, genparam.fipsToPolys,
                                                  ncores_py=genparam.ncores_py)
        print('Got RE CFs. Elapsed time: ' + str_elapsedtime(start_time))

        start_time = time.time()
        if windCFs is not None:
            write2dListToCSV(windCFs, os.path.join(resultsDir, 'windCFsFullYrCE' + zone + str(currYear) + '.csv'))
            write2dListToCSV(windCfsDtHr, os.path.join(resultsDir, 'windCFsDtFullYrCE' + zone + str(currYear) + '.csv'))
            # write2dListToCSV(windCfsDtSubhr,os.path.join(resultsDir,'windCFsDtSubhrFullYrCE' + zone + str(currYear) + '.csv'))
            write2dListToCSV(ewdIdAndCapac,
                             os.path.join(resultsDir, 'windIdAndCapacCE' + zone + str(currYear) + '.csv'))
        if solarCFs is not None:
            write2dListToCSV(solarCFs, os.path.join(resultsDir, 'solarCFsFullYrCE' + zone + str(currYear) + '.csv'))
            write2dListToCSV(solarCfsDtHr,
                             os.path.join(resultsDir, 'solarCFsDtFullYrCE' + zone + str(currYear) + '.csv'))
            # write2dListToCSV(solarCfsDtSubhr,os.path.join(resultsDir,'solarCFsDtSubhrFullYrCE' + zone + str(currYear) + '.csv'))
            write2dListToCSV(solarFilenameAndCapac,
                             os.path.join(resultsDir, 'solarIdAndCapacCE' + zone + str(currYear) + '.csv'))

        (newWindCFs, newWindCFsSubhr, newSolarCFs, newSolarCFsSubhr, newWindCfsDtHr, newWindCfsDtSubhr,
         newWindIdAndCapac, newSolarCfsDtHr, newSolarCfsDtSubhr, newSolarFilenameAndCapac,
         addedWindCapac, addedSolarCapac) = getNewWindAndSolarCFs(zonalGenFleet, currYear, 'CE', genparam.tzAnalysis,
                                                                  genparam.dataRoot, resultsDir, genparam.windGenDataYr,
                                                                  zone, genparam.fipsToZones, genparam.fipsToPolys,
                                                                  ncores_py=genparam.ncores_py)
        print('Got new RE CFs. Elapsed time: ' + str_elapsedtime(start_time))

        write2dListToCSV(newWindIdAndCapac,
                         os.path.join(resultsDir, 'windNewIdAndCapacCE' + zone + str(currYear) + '.csv'))
        write2dListToCSV(newSolarFilenameAndCapac,
                         os.path.join(resultsDir, 'solarNewIdAndCapacCE' + zone + str(currYear) + '.csv'))

        # get number of hours in year (assume that is the same over all gcms)
        gcm = curtailparam.listgcms[0]
        nhours_year = len(zonalDemandProfile[gcm][zone])

        zonalHourlyWindGen[zone], zonalHourlySolarGen[zone] = getAggregateSolarAndWind(windCFs, ewdIdAndCapac,
                                                                                       solarCFs, solarFilenameAndCapac,
                                                                                       nhours_year)
        zonalNewWindCFs[zone], zonalNewSolarCFs[zone] = newWindCFs, newSolarCFs

    zonalNetDemand = dict()
    for gcm in curtailparam.listgcms:
        auxDemandDict = dict()
        for zone in genparam.ipmZones:
            netDemand = getNetDemand(zonalDemandProfile[gcm][zone], zonalHourlyWindGen[zone], zonalHourlySolarGen[zone])
            auxDemandDict[zone] = netDemand
        zonalNetDemand[gcm] = auxDemandDict

    writeDictToCSV(zonalHourlyWindGen, os.path.join(resultsDir, 'windGenFullYrCE' + str(currYear) + '.csv'))
    writeDictToCSV(zonalHourlySolarGen, os.path.join(resultsDir, 'solarGenFullYrCE' + str(currYear) + '.csv'))
    writeDictToCSV(zonalNewWindCFs, os.path.join(resultsDir, 'windNewCFsFullYrCE' + str(currYear) + '.csv'))
    writeDictToCSV(zonalNewSolarCFs, os.path.join(resultsDir, 'solarNewCFsFullYrCE' + str(currYear) + '.csv'))
    writeDictToCSV(zonalNetDemand, os.path.join(resultsDir, 'demandNetFullYrCE' + str(currYear) + '.csv'))

    start_time = time.time()
    print('Computing hourly curtailments for existing generators... ', flush=True)

    hrlyCurtailmentsAllGensInTgtYr = importHourlyThermalCurtailments(genFleetForCE, currYear, 'CE', resultsDir,
                                                                     genparam, curtailparam)
    print('Done! Elapsed time: ' + str_elapsedtime(start_time))

    start_time = time.time()

    print('Selecting subset of hours for CE model... ', flush=True)
    (demandCE, hoursForCE, repHrsBySeason, specialHrs, regHrsBySeason, demandCEZonal,
     hourlyWindGenCEZonal, hourlySolarGenCEZonal, peakDemandHourZonal, planningReserveZonal,
     repAndSpeHoursDict) = selectWeeksForExpansion(zonalDemandProfile, zonalNetDemand, zonalHourlyWindGen,
                                                   zonalHourlySolarGen, genparam.daysPerSeason,
                                                   genparam.selectCurtailDays, hrlyCurtailmentsAllGensInTgtYr,
                                                   currYear, resultsDir, genparam.planningReserve)

    writeDictToCSV(hoursForCE, os.path.join(resultsDir, 'hoursCE' + str(currYear) + '.csv'))
    writeDictToCSV(demandCE, os.path.join(resultsDir, 'demandCE' + str(currYear) + '.csv'))
    writeDictToCSV(demandCEZonal, os.path.join(resultsDir, 'demandZonalCE' + str(currYear) + '.csv'))
    writeDictToCSV(hourlyWindGenCEZonal, os.path.join(resultsDir, 'windGenZonalCE' + str(currYear) + '.csv'))
    writeDictToCSV(hourlySolarGenCEZonal, os.path.join(resultsDir, 'solarGenZonalCE' + str(currYear) + '.csv'))
    writeDictToCSV(planningReserveZonal, os.path.join(resultsDir, 'planningReserveZonalCE' + str(currYear) + '.csv'))
    writeDictToCSV(peakDemandHourZonal, os.path.join(resultsDir, 'peakDemandHoursZonalCE' + str(currYear) + '.csv'))

    seasonDemandWeights, weightsList = calculateSeasonalWeights(zonalDemandProfile, repHrsBySeason, regHrsBySeason)

    write2dListToCSV(weightsList, os.path.join(resultsDir, 'seasonWeightsCE' + str(currYear) + '.csv'))

    (newWindCFsCEZonal, newSolarCFsCEZonal) = trimNewRECFsToCEHours(zonalNewWindCFs, zonalNewSolarCFs, hoursForCE)

    writeDictToCSV(newWindCFsCEZonal, os.path.join(resultsDir, 'windNewCFsZonalCE' + str(currYear) + '.csv'))
    writeDictToCSV(newSolarCFsCEZonal, os.path.join(resultsDir, 'solarNewCFsZonalCE' + str(currYear) + '.csv'))

    print('Done! Elapsed time: ' + str_elapsedtime(start_time))

    start_time = time.time()
    print('Importing hydro generation data... ', flush=True)
    # Import hydropower generation potential (dict of {gcm: {season: {ORIS ID: potential gen}}})
    hydroPotPerSeason = getHydroEPotential(genFleetForCE, zonalDemandProfile, repAndSpeHoursDict, currYear, genparam,
                                           curtailparam)

    writeDictToCSV(hydroPotPerSeason, os.path.join(resultsDir, 'hydroPot' + 'CE' + str(currYear) + '.csv'))

    print('Done! Elapsed time: ' + str_elapsedtime(start_time))

    hourlyCapacsCE = getHourlyNonRECapacsForCE(genFleetForCE, hrlyCurtailmentsAllGensInTgtYr, hoursForCE, currYear)
    writeDictToCSV(hourlyCapacsCE, os.path.join(resultsDir, 'curtailedHourlyCapacsCE' + str(currYear) + '.csv'))

    start_time = time.time()
    print('Computing hourly curtailments for new generators...', flush=True)
    eligibleCellWaterTs = loadEligibleCellWaterTs(genFleet, currYear, genparam, curtailparam)

    print('  Got eligible cells: ' + str_elapsedtime(start_time))

    # hrlyCurtailmentsAllTechsInTgtYr is a dict of (tech,cell):[hourlycapac]
    hrlyCurtailmentsAllTechsInTgtYr = \
        determineHrlyCurtailmentsForNewTechs(eligibleCellWaterTs, newTechsCE, currYear, genparam, curtailparam,
                                             resultsDir)

    writeDictToCSV(hrlyCurtailmentsAllTechsInTgtYr, os.path.join(resultsDir, 'curtailmentsTechHourly' +
                                                                 str(currYear) + '.csv'))
    hourlyCurtailedTechCapacsCE = getHourlyCurtailedTechCapacsForCE(newTechsCE, hrlyCurtailmentsAllTechsInTgtYr,
                                                                    hoursForCE, currYear, genparam.ptCurtailedAll)

    writeDictToCSV(hourlyCurtailedTechCapacsCE, os.path.join(resultsDir,
                                                             'curtailedTechHourlyCapacsCE' + str(currYear) + '.csv'))
    print('Done! Elapsed time: ' + str_elapsedtime(start_time))

    print('Got hourly generation profiles with curtailments for existing and new techs')
    # Extract set of cells for techs & assign to zones
    cellsForNewTechs = set()

    for (tech, cell) in hrlyCurtailmentsAllTechsInTgtYr[list(hrlyCurtailmentsAllTechsInTgtYr.keys())[0]]:
        cellsForNewTechs.add(cell)

    cellsToZones = assignCellsToIPMZones(cellsForNewTechs, genparam.fipsToZones, genparam.fipsToPolys)
    writeDictToCSV(cellsToZones, os.path.join(resultsDir, 'cellsToZonesCE' + str(currYear) + '.csv'))

    # free memory before calling GAMS/CPLEX!
    del hrlyCurtailmentsAllTechsInTgtYr
    del eligibleCellWaterTs
    del hrlyCurtailmentsAllGensInTgtYr
    del zonalHourlyWindGen
    del zonalHourlySolarGen
    del zonalNewSolarCFs
    del zonalNewWindCFs

    gc.collect()

    # Run CE
    print('Set inputs, running CE model...')
    t0 = time.time()
    (capacExpModel, ms, ss, curtailedTechs, renewTechs, notCurtailedTechs,
     pumpedHydroSymbs) = callCapacityExpansion(genFleetForCE, hourlyCapacsCE, hourlyCurtailedTechCapacsCE,
                                               hourlyWindGenCEZonal, hourlySolarGenCEZonal, demandCEZonal, newTechsCE,
                                               planningReserveZonal, hoursForCE, newWindCFsCEZonal, newSolarCFsCEZonal,
                                               currCo2Cap, seasonDemandWeights, repHrsBySeason, specialHrs,
                                               peakDemandHourZonal, cellsToZones, hydroPotPerSeason, genparam)

    write2dListToCSV([['ms', 'ss'], [ms, ss]], os.path.join(resultsDir, 'msAndSsCE' + str(currYear) + '.csv'))
    print('Time (secs) for CE year {0:4d}: {1:.2f}'.format(currYear, (time.time() - t0)))

    # copy input and output GDX files to output folder and rename them
    gamsFileDir = os.path.join(genparam.resultsDir, 'GAMS')
    try:
        shutil.copy(os.path.join(gamsFileDir, '_gams_py_gdb0.gdx'),
                    os.path.join(resultsDir, 'gdxInYear{}.gdx'.format(currYear)))
        shutil.copy(os.path.join(gamsFileDir, '_gams_py_gdb1.gdx'),
                    os.path.join(resultsDir, 'gdxOutYear{}.gdx'.format(currYear)))
    except IOError as e:
        print("Unable to copy file. %s" % e)
    except:
        print("Unexpected error:", sys.exc_info())

    # Write and Save resulting data
    print()
    print('Writing resulting decision variables to csv files... ', end='', flush=True)
    (genByPlant, genByCTechAndCell, genByRETechAndZone, genByNCTechAndZone, flowByLine, chargeByPH,
     socByPH) = saveCapacExpOperationalData(capacExpModel)

    write2dListToCSV(genByPlant, os.path.join(resultsDir, 'genByPlantCE' + str(currYear) + '.csv'))
    write2dListToCSV(genByCTechAndCell, os.path.join(resultsDir, 'genByCurtTechCE' + str(currYear) + '.csv'))
    write2dListToCSV(genByRETechAndZone, os.path.join(resultsDir, 'genByRETechCE' + str(currYear) + '.csv'))
    write2dListToCSV(genByNCTechAndZone, os.path.join(resultsDir, 'genByNotCurtTechCE' + str(currYear) + '.csv'))
    # write2dListToCSV(sysResults, os.path.join(resultsDir, 'systemResultsCE' + str(currYear) + '.csv'))
    # write2dListToCSV(co2EmAndCostResults, os.path.join(resultsDir, 'co2EmsAndCostsCE' + str(currYear) + '.csv'))
    write2dListToCSV(flowByLine, os.path.join(resultsDir, 'lineFlowsCE' + str(currYear) + '.csv'))
    write2dListToCSV(chargeByPH, os.path.join(resultsDir, 'chargeByPumpHydroCE' + str(currYear) + '.csv'))
    write2dListToCSV(socByPH, os.path.join(resultsDir, 'socByPumpHydroCE' + str(currYear) + '.csv'))
    print('Done!')
    print()

    capacExpModelsEachYear.append((currYear, capacExpModel))
    newCurtTech, newRETech, newNotCurtTech = saveCapacExpBuilds(capacExpBuilds, capacExpModel, currYear)

    genFleet = addNewGensToFleet(genFleet, newCurtTech, newRETech, newNotCurtTech, newTechsCE, currYear,
                                 genparam.ipmZones, genparam.ipmZoneNums, genparam.ocAdderMin, genparam.ocAdderMax,
                                 cellsToZones, genparam.ptCurtailedAll, genparam.statePolys)

    genFleet = selectAndMarkUnitsRetiredByCE(genFleet, genFleetForCE, genparam.retirementCFCutoff, capacExpModel,
                                             currYear, capacExpGenByGens, capacExpRetiredUnitsByCE,
                                             genparam.scaleMWtoGW, hoursForCE, planningReserveZonal, genparam.endYear,
                                             capacExpRetiredUnitsByAge, demandCEZonal, hourlyWindGenCEZonal,
                                             hourlySolarGenCEZonal, newWindCFsCEZonal, newSolarCFsCEZonal,
                                             genparam.ptEligRetCF, peakDemandHourZonal)

    # removes all retired units; [] is dummy list b/c not adding ret age units to list
    genFleetNoRetiredUnits = createFleetForCurrentCELoop(genFleet, currYear, [], genparam.dataRoot, genparam.scenario)

    writeCEInfoToCSVs(capacExpBuilds, capacExpGenByGens, capacExpRetiredUnitsByCE,
                      capacExpRetiredUnitsByAge, resultsDir, currYear)
    write2dListToCSV(genFleet, os.path.join(resultsDir, 'genFleetAfterCE' + str(currYear) + '.csv'))

    # Write gen fleet for UC to special folder for ease of transfer
    ceUCDir = os.path.join(genparam.resultsDir, 'CEtoUC')

    if not os.path.exists(ceUCDir): os.makedirs(ceUCDir)
    write2dListToCSV(genFleetNoRetiredUnits, os.path.join(ceUCDir, 'genFleetCEtoUC' + str(currYear) + '.csv'))

    return genFleet, genFleetNoRetiredUnits, genFleetForCE, capacExpModel.out_db, hoursForCE


def importHourlyThermalCurtailments(genFleet, currYear, modelName, resultsDir, genparam, curtailparam):
    """IMPORT HOURLY THERMAL CURTAILMENTS BY GENERATOR

    :param genFleet: 2d list (matrix) with data of gen fleet
    :param currYear: current year of analysis
    :param modelName: string with model name ('CE' or 'UC')
    :param resultsDir: string specifying folder where files will be saved (it is a subfolder of genparam.dataRoot)
    :param genparam: object of class General Parameters
    :param curtailparam: object of class General Parameters
    :return: nested dict mapping each gen to 1d numpy array of hourly net capacity curtailments(MW)
             dict is {gcm:{genId: 1d np.array()}}
    """
    (hrlyCurtailmentsAllGensInTgtYr,
     genToCellLatLongsList) = determineHrlyCurtailmentsForExistingGens(genFleet, currYear, genparam, curtailparam)

    write2dListToCSV(genToCellLatLongsList, os.path.join(resultsDir,
                                                         'mapGensToCells' + modelName + str(currYear) + '.csv'))

    # convert 1d numpy array to list in order to print to csv file
    auxDict = {}

    for gcm in hrlyCurtailmentsAllGensInTgtYr.keys():
        auxDict2 = {}
        for gen in hrlyCurtailmentsAllGensInTgtYr[gcm].keys():
            auxDict2[gen] = list(hrlyCurtailmentsAllGensInTgtYr[gcm][gen])
        auxDict[gcm] = auxDict2

    writeDictToCSV(auxDict, os.path.join(resultsDir, 'curtailmentsHourly' + str(currYear) + '.csv'))

    return hrlyCurtailmentsAllGensInTgtYr


def callCapacityExpansion(genFleetForCE, hourlyCapacsCE, hourlyCurtailedTechCapacsCE, hourlyWindGenCEZonal,
                          hourlySolarGenCEZonal, demandCEZonal, newTechsCE, planningReserveZonal, hoursForCE,
                          newWindCFsCEZonal, newSolarCFsCEZonal, currCo2Cap, seasonDemandWeights, repHrsBySeason,
                          specialHrs, peakDemandHourZonal, cellsToZones, hydroPotPerSeason, genparam):
    """CALL CAPACITY EXPANSION

    :param genFleetForCE:
    :param hourlyCapacsCE:
    :param hourlyCurtailedTechCapacsCE:
    :param hourlyWindGenCEZonal:
    :param hourlySolarGenCEZonal:
    :param demandCEZonal:
    :param newTechsCE:
    :param planningReserveZonal:
    :param hoursForCE:
    :param newWindCFsCEZonal:
    :param newSolarCFsCEZonal:
    :param currCo2Cap:
    :param seasonDemandWeights:
    :param repHrsBySeason:
    :param specialHrs:
    :param peakDemandHourZonal:
    :param cellsToZones:
    :param hydroPotPerSeason:
    :param genparam:
    :return:
    """

    currDir = os.getcwd()

    # temporarily assign fields of general parameters to variables (old way)
    (discountRate, scaleMWtoGW, scaleDollarsToThousands, capacExpFilename, scaleLbToShortTon, dataRoot,
     maxAddedZonalCapacPerTech, pathSysGams, ipmZones, lineList, lineCapacs, ipmZoneNums, ptCurtailedAll,
     phEff, phMaxSoc, phInitSoc) = (genparam.discountRate, genparam.scaleMWtoGW, genparam.scaleDollarsToThousands,
                                    genparam.capacExpFilename, genparam.scaleLbToShortTon, genparam.dataRoot,
                                    genparam.maxAddedZonalCapacPerTech, genparam.pathSysGams, genparam.ipmZones,
                                    genparam.lineList, genparam.lineCapacs, genparam.ipmZoneNums,
                                    genparam.ptCurtailedAll, genparam.phEff, genparam.phMaxSoc, genparam.phInitSoc)

    # folder with gams file for this instance (this is why this is set to resultsdir)
    gamsFileDir = os.path.join(genparam.resultsDir, 'GAMS')
    gamsSysDir = pathSysGams.strip()

    if not gamsSysDir == '':
        ws = GamsWorkspace(working_directory=gamsFileDir, system_directory=gamsSysDir)
    else:
        ws = GamsWorkspace(working_directory=gamsFileDir)

    db = ws.add_database()

    # Add sets and parameters to database
    (genSet, genSymbols, hourSet, hourSymbols, techSet, techSymbols, techCurtailedSet, techCurtailedSymbols,
     renewTechSet, renewTechSymbols, techNotCurtailedSymbols, hydroGenSet, hydroGenSymbols, zoneSet, zoneSymbols,
     lineSet, cellSet, peakHourSet, peakHrSymbols, pumpHydroGenSet, pumpHydroGenSymbols, typeSet, gcmSet) = \
        addSetsToDatabaseCE(db, genFleetForCE, hoursForCE, newTechsCE, repHrsBySeason, specialHrs, peakDemandHourZonal,
                            cellsToZones, genparam)

    addParametersToDatabaseCE(db, hourlyCapacsCE, hourlyWindGenCEZonal, hourlySolarGenCEZonal, demandCEZonal,
                              newTechsCE, genFleetForCE, genSet, genSymbols, hydroGenSet, hoursForCE, hourSet,
                              hourSymbols, techSet, techSymbols, techCurtailedSet, techCurtailedSymbols, renewTechSet,
                              renewTechSymbols, planningReserveZonal, discountRate, newWindCFsCEZonal,
                              newSolarCFsCEZonal, scaleMWtoGW, scaleDollarsToThousands, currCo2Cap,
                              seasonDemandWeights, hourlyCurtailedTechCapacsCE, scaleLbToShortTon,
                              maxAddedZonalCapacPerTech, zoneSet, zoneSymbols, lineSet, ipmZones, ipmZoneNums,
                              lineList, lineCapacs, cellsToZones, cellSet, peakDemandHourZonal, peakHourSet,
                              peakHrSymbols, hydroPotPerSeason, ptCurtailedAll, pumpHydroGenSet, pumpHydroGenSymbols,
                              phEff, phMaxSoc, phInitSoc, typeSet, gcmSet)
    # Load GAMS model
    capacExpFile = capacExpFilename
    capacExpModel = ws.add_job_from_file(capacExpFile)
    # Run GAMS model
    opt = GamsOptions(ws)
    opt.defines['gdxincname'] = db.name

    # get numbers of cores
    n_cores = genparam.ncores_gams

    # add CPLEX option file
    file = open(os.path.join(ws.working_directory, "cplex.opt"), "w")
    file.write("mipdisplay 4\n")
    file.write("mipinterval 10\n")
    file.write("threads {}\n".format(n_cores))
    file.close()
    opt.optfile = 1

    # Run GAMS/CPLEX model
    capacExpModel.run(opt, databases=db, output=sys.stdout)

    ms = int(capacExpModel.out_db['pModelstat'].find_record().value)
    ss = int(capacExpModel.out_db['pSolvestat'].find_record().value)

    # print status of solution
    print('End of CE Model!')

    ms_description = list(dict(ModelStat.__dict__).keys())[list(dict(ModelStat.__dict__).values()).index(ms)]
    ss_description = list(dict(SolveStat.__dict__).keys())[list(dict(SolveStat.__dict__).values()).index(ss)]

    print()
    print('Model Status:  {0:4d} ({1})'.format(ms, ms_description))
    print('Solver Status: {0:4d} ({1})'.format(ss, ss_description))
    print()

    path_gams = gamsFileDir
    str_target = 'SET z  zones'
    fname = os.path.join(path_gams, '_gams_py_gjo0.lst')

    with open(fname) as f:
        lines = [line.rstrip('\n') for line in f]

    idx_set_zone = [i for i in range(len(lines)) if str_target in lines[i]]

    if len(idx_set_zone) > 0:
        idx_set_zone = idx_set_zone[0]
        str_zones = lines[idx_set_zone + 2]

        print()
        print('-------- Order of Zone set used in GAMS: ------------')
        print(str_zones)
        print()

    else:
        print()
        print('---- DISPLAY OF ZONE SET NOT FOUND IN LST FILE ------')
        print()

    if ss != 1:
        print()
        print('-------- INVALID SOLUTION FOUND BY SOLVER ------------')
        print('Exiting simulation...')
        print()
        print()
        sys.exit()

    return capacExpModel, ms, ss, techCurtailedSymbols, renewTechSymbols, techNotCurtailedSymbols, pumpHydroGenSymbols


################################### ADD SETS
def addSetsToDatabaseCE(db, genFleetForCE, hoursForCE, newTechsCE, repHrsBySeason, specialHrs, peakDemandHourZonal,
                        cellsToZones, genparam):
    ipmZoneNums, lineList, ptCurtailedAll = genparam.ipmZoneNums, genparam.lineList, genparam.ptCurtailedAll

    zoneSet, zoneSymbols = addZoneSets(db, ipmZoneNums)  # create string for set IDs
    lineSet = addLineSets(db, lineList)
    cellSet = addCellSet(db, cellsToZones)

    # compile 1d list with all hours for CE model
    hoursForCETotal = get_list_all_hours_ce(hoursForCE)

    (hourSet, hourSymbols) = addHourSet(db, hoursForCETotal)

    (gcmSetName, gcmSetDescrip, gcmSetDim) = ('g', 'gcms scenarios', 1)
    gcmSet = db.add_set(gcmSetName, gcmSetDim, gcmSetDescrip)
    for gcm in hoursForCE:
        gcmSet.add_record(gcm)

    (genSet, genSymbols, hydroGenSet, hydroGenSymbols, pumpHydroGenSet,
     pumpHydroGenSymbols) = addGeneratorSets(db, genFleetForCE)

    # write 2d set with hours for CE in each gcm
    (hours2dSetName, hours2dSetDescrip, hours2dSetDim) = ('h2', '2d set mapping gcm to hours for CE', 2)
    hours2dSet = db.add_set(hours2dSetName, hours2dSetDim, hours2dSetDescrip)
    for gcm in hoursForCE:
        for h in hoursForCE[gcm]:
            hsymbol = createHourSymbol(h)
            hours2dSet.add_record([gcm, hsymbol])

    addHourSeasonSubsets(db, repHrsBySeason)

    addHourSpecialSubset(db, specialHrs)

    (techSet, techSymbols, techCurtailedSet, techCurtailedSymbols, renewTechSet, renewTechSymbols,
     techNotCurtailedSet, techNotCurtailedSymbols) = addNewTechsSets(db, newTechsCE, ptCurtailedAll)

    # create additional set with type of plant (no cooling info)
    typeSymbols = set(map(lambda x: x.split('+')[0], techSymbols))
    (typeSetName, typeSetDescrip, typeSetDim) = ('type', 'set with types of new plants (no cooling info)', 1)
    typeSet = addSet(db, typeSymbols, typeSetName, typeSetDescrip, typeSetDim)

    # create additional 2d set mapping type of plant to plant+cooling
    typeList = list(map(lambda x: x.split('+')[0], techSymbols))
    (tech2dSetName, tech2dSetDescrip, tech2dSetDim) = ('tech2d', '2d set mapping types to techs', 2)
    techSet2 = db.add_set(tech2dSetName, tech2dSetDim, tech2dSetDescrip)
    for t1 in typeSymbols:
        idx = [i for i, x in enumerate(typeList) if x == t1]
        list_techs = [techSymbols[i] for i in idx]
        for t2 in list_techs:
            techSet2.add_record([t1, t2])

    peakHourSet, peakHrSymbols = addPeakHourSubset(db, peakDemandHourZonal, genparam)

    return (genSet, genSymbols, hourSet, hourSymbols, techSet, techSymbols, techCurtailedSet,
            techCurtailedSymbols, renewTechSet, renewTechSymbols, techNotCurtailedSymbols,
            hydroGenSet, hydroGenSymbols, zoneSet, zoneSymbols, lineSet, cellSet, peakHourSet,
            peakHrSymbols, pumpHydroGenSet, pumpHydroGenSymbols, typeSet, gcmSet)


################################### ADD PARAMETERS
def addParametersToDatabaseCE(db, hourlyCapacsCE, hourlyWindGenCEZonal, hourlySolarGenCEZonal, demandCEZonal,
                              newTechsCE, genFleetForCE, genSet, genSymbols, hydroGenSet, hoursForCE, hourSet,
                              hourSymbols, techSet, techSymbols, techCurtailedSet, techCurtailedSymbols, renewTechSet,
                              renewTechSymbols, planningReserveZonal, discountRate, newWindCFsCEZonal,
                              newSolarCFsCEZonal, scaleMWtoGW, scaleDollarsToThousands, currCo2Cap, seasonDemandWeights,
                              hourlyCurtailedTechCapacsCE, scaleLbToShortTon, maxAddedZonalCapacPerTech, zoneSet,
                              zoneSymbols, lineSet, ipmZones, ipmZoneNums, lineList, lineCapacs, cellsToZones,
                              cellSet, peakDemandHourZonal, peakHourSet, peakHrSymbols, hydroPotPerSeason,
                              ptCurtailedAll, pumpHydroGenSet, pumpHydroGenSymbols, phEff, phMaxSoc, phInitSoc,
                              typeSet, gcmSet):
    addTechParams(db, newTechsCE, scaleMWtoGW, scaleDollarsToThousands, scaleLbToShortTon, ptCurtailedAll)

    addEguParams(db, genFleetForCE, genSet, genSymbols, ipmZones, ipmZoneNums, scaleLbToShortTon, scaleMWtoGW)
    addPumpHydroParams(db, genFleetForCE, phEff, phMaxSoc, phInitSoc, pumpHydroGenSet, pumpHydroGenSymbols, scaleMWtoGW)
    addDiscountRateParam(db, discountRate)
    addCppEmissionsCap(db, currCo2Cap)
    addSeasonDemandWeights(db, seasonDemandWeights)

    addMaxNumNewBuilds(db, newTechsCE, zoneSet, ipmZones, ipmZoneNums, typeSet, maxAddedZonalCapacPerTech,
                       ptCurtailedAll)

    # Add hydro max gen per time block
    addHydroMaxGenPerSeason(db, hydroGenSet, gcmSet, hydroPotPerSeason, scaleMWtoGW)

    # Add zonal demand, planning margin, and RE generation
    addDemandParam(db, demandCEZonal, zoneSet, hourSet, gcmSet, hoursForCE, ipmZones, ipmZoneNums, scaleMWtoGW)

    addPlanningReserveParam(db, planningReserveZonal, ipmZones, ipmZoneNums, zoneSet, zoneSymbols, scaleMWtoGW)
    #    addPeakHourToZoneParam(db, peakDemandHourZonal, peakHourSet, peakHrSymbols, ipmZones, ipmZoneNums)

    addExistingRenewableMaxGenParams(db, gcmSet, zoneSet, ipmZones, ipmZoneNums, hourSet, hoursForCE,
                                     hourlySolarGenCEZonal, hourlyWindGenCEZonal, scaleMWtoGW)

    addRenewTechCFParams(db, renewTechSet, renewTechSymbols, gcmSet, zoneSet, hourSet, hoursForCE,
                         newWindCFsCEZonal, newSolarCFsCEZonal, ipmZones, ipmZoneNums)

    # Add hourly capacity for curtailed techs
    addTechCurtailedHourlyCapac(db, hourlyCurtailedTechCapacsCE, gcmSet, cellSet, techCurtailedSet, hourSet, hoursForCE,
                                scaleMWtoGW)

    # Add tech heat rate and op cost params that reflect degradation
    addEguHourlyParams(db, hourlyCapacsCE, gcmSet, genSet, hourSet, hoursForCE, scaleMWtoGW)

    addEguOpCostParam(db, genFleetForCE, genSet, genSymbols, scaleLbToShortTon, scaleMWtoGW, scaleDollarsToThousands)

    # Add zone & line parameters
    addLineCapacs(db, lineCapacs, lineSet, lineList, scaleMWtoGW)
    addLineSourceAndSink(db, lineSet, lineList, ipmZones, ipmZoneNums)
    addCellsToZones(db, cellSet, cellsToZones, ipmZones, ipmZoneNums)


def runUnitCommitment(genFleet, zonalDemandProfile, ucYear, currCo2Cap, genparam, reserveparam, curtailparam):
    """RUN UNIT COMMITMENT

    :param genFleet:
    :param zonalDemandProfile:
    :param ucYear:
    :param currCo2Cap:
    :param genparam:
    :param reserveparam:
    :param curtailparam:
    :return:
    """

    daysForUC = list(range(genparam.ucDayInitial, (genparam.ucDayEnd + 1)))

    (startYear, fuelPricesTimeSeries, scaleMWtoGW, scaleDollarsToThousands, scaleLbToShortTon,
     daysOpt, daysLA, projectName, resultsDir, ocAdderMin, ocAdderMax,
     windGenDataYr) = (genparam.startYear, genparam.fuelPricesTimeSeries, genparam.scaleMWtoGW,
                       genparam.scaleDollarsToThousands, genparam.scaleLbToShortTon, genparam.daysOpt, genparam.daysLA,
                       genparam.projectName, genparam.resultsDir, genparam.ocAdderMin, genparam.ocAdderMax,
                       genparam.windGenDataYr)

    (regLoadFrac, contLoadFrac, regErrorPercentile, flexErrorPercentile, rrToRegTime, rrToFlexTime, rrToContTime,
     regUpCostCoeffs) = (reserveparam.regLoadFrac, reserveparam.contLoadFrac, reserveparam.regErrorPercentile,
                         reserveparam.flexErrorPercentile, reserveparam.rampRateToRegReserveScalar,
                         reserveparam.rampRateToFlexReserveScalar, reserveparam.rampRateToContReserveScalar,
                         copy.deepcopy(reserveparam.regUpCostCoeffs))

    resultsDir = os.path.join(resultsDir, 'UC')
    if not os.path.exists(resultsDir): os.makedirs(resultsDir)
    print('Entering UC loop for year ' + str(ucYear))

    fleetUC = copy.deepcopy(genFleet)

    (startWindCapacForCFs, startSolarCapacForCFs) = (0, 0)
    zoneCol = fleetUC[0].index('Region Name')
    zonalHourlyWindGen, zonalHourlySolarGen = dict(), dict()

    print('Computing CFs for Renewables...')
    for zone in genparam.ipmZones:
        print('Zone ', zone)

        start_time = time.time()
        zonalGenFleet = [fleetUC[0]] + [row for row in fleetUC if row[zoneCol] == zone]

        (windCFs, windCfsDtHr, windCfsDtSubhr, windIdAndCapac, solarCFs, solarCfsDtHr, solarCfsDtSubhr,
         solarFilenameAndCapac) = getRenewableCFs(zonalGenFleet, startWindCapacForCFs, startSolarCapacForCFs,
                                                  genparam.tzAnalysis, genparam.dataRoot, genparam.windGenDataYr,
                                                  zone, genparam.fipsToZones, genparam.fipsToPolys, subHour=True,
                                                  ncores_py=genparam.ncores_py)

        print('Got RE CFs. Elapsed time: ' + str_elapsedtime(start_time))

        start_time = time.time()
        if windCFs is not None:
            write2dListToCSV(windCFs, os.path.join(resultsDir, 'windCFsFullYrCE' + zone + str(ucYear) + '.csv'))
            write2dListToCSV(windCfsDtHr, os.path.join(resultsDir, 'windCFsDtFullYrCE' + zone + str(ucYear) + '.csv'))
            write2dListToCSV(windCfsDtSubhr,
                             os.path.join(resultsDir, 'windCFsDtSubhrFullYrCE' + zone + str(ucYear) + '.csv'))
            write2dListToCSV(windIdAndCapac, os.path.join(resultsDir, 'windIdAndCapacCE' + zone + str(ucYear) + '.csv'))
        if solarCFs is not None:
            write2dListToCSV(solarCFs, os.path.join(resultsDir, 'solarCFsFullYrCE' + zone + str(ucYear) + '.csv'))
            write2dListToCSV(solarCfsDtHr, os.path.join(resultsDir, 'solarCFsDtFullYrCE' + zone + str(ucYear) + '.csv'))
            write2dListToCSV(solarCfsDtSubhr,
                             os.path.join(resultsDir, 'solarCFsDtSubhrFullYrCE' + zone + str(ucYear) + '.csv'))
            write2dListToCSV(solarFilenameAndCapac,
                             os.path.join(resultsDir, 'solarIdAndCapacCE' + zone + str(ucYear) + '.csv'))

        # get number of hours in year (assume that is the same over all gcms)
        gcm = curtailparam.listgcms[0]
        nhours_year = len(zonalDemandProfile[gcm][zone])

        zonalHourlyWindGen[zone], zonalHourlySolarGen[zone] = getAggregateSolarAndWind(windCFs, windIdAndCapac,
                                                                                       solarCFs, solarFilenameAndCapac,
                                                                                       nhours_year)

    #    (windCFs, windCfsDtHr, windCfsDtSubhr, windIdAndCapac, solarCFs, solarCfsDtHr, solarCfsDtSubhr,
    #     solarFilenameAndCapac) = getRenewableCFs(fleetUC, startWindCapacForCFs, startSolarCapacForCFs, tzAnalysis, dataRoot,
    #                                              windGenDataYr, z)
    #     write2dListToCSV([zonalDemandProfile], os.path.join(resultsDir, 'demandUC' + str(ucYear) + '.csv'))
    #     write2dListToCSV(windCFs, os.path.join(resultsDir, 'windCFsUC' + str(ucYear) + '.csv'))
    #     write2dListToCSV(windCfsDtHr, os.path.join(resultsDir, 'windCFsDtUC' + str(ucYear) + '.csv'))
    #     write2dListToCSV(windCfsDtSubhr, os.path.join(resultsDir, 'windCFsDtSubhrUC' + str(ucYear) + '.csv'))
    #     write2dListToCSV(windIdAndCapac, os.path.join(resultsDir, 'windIdAndCapacUC' + str(ucYear) + '.csv'))
    #     write2dListToCSV(solarCFs, os.path.join(resultsDir, 'solarCFsUC' + str(ucYear) + '.csv'))
    #     write2dListToCSV(solarCfsDtHr, os.path.join(resultsDir, 'solarCFsDtUC' + str(ucYear) + '.csv'))
    #     write2dListToCSV(solarCfsDtSubhr, os.path.join(resultsDir, 'solarCFsDtSubhrUC' + str(ucYear) + '.csv'))
    #     write2dListToCSV(solarFilenameAndCapac, os.path.join(resultsDir, 'solarIdAndCapacUC' + str(ucYear) + '.csv'))

    (contResHourly, regUpHourly, regDownHourly, flexResHourly, allRes, regUpWind, regDownWind, regUpSolar,
     regDownSolar, flexWind, flexSolar) = (dict(), dict(), dict(), dict(), dict(), dict(), dict(), dict(), dict(),
                                           dict(), dict())

    for gcm in curtailparam.listgcms:

        (auxcontResHourly, auxregUpHourly, auxregDownHourly, auxflexResHourly, auxallRes, auxregUpWind, auxregDownWind,
         auxregUpSolar, auxregDownSolar, auxflexWind, auxflexSolar) = (dict(), dict(), dict(), dict(), dict(), dict(),
                                                                       dict(), dict(), dict(), dict(), dict())

        for zone in genparam.ipmZones:
            (zonalcontResHourly, zonalregUpHourly, zonalregDownHourly, zonalflexResHourly, zonalallRes, zonalregUpWind,
             zonalregDownWind, zonalregUpSolar, zonalregDownSolar, zonalflexWind,
             zonalflexSolar) = calcWWSISReserves(windCfsDtHr, windCfsDtSubhr, windIdAndCapac, solarCfsDtHr,
                                                 solarCfsDtSubhr, solarFilenameAndCapac, zonalDemandProfile[gcm][zone],
                                                 regLoadFrac, contLoadFrac, regErrorPercentile, flexErrorPercentile)
            auxcontResHourly[zone] = zonalcontResHourly
            auxregUpHourly[zone] = zonalregUpHourly
            auxregDownHourly[zone] = zonalregDownHourly
            auxflexResHourly[zone] = zonalflexResHourly
            auxallRes[zone] = zonalallRes
            auxregUpWind[zone] = zonalregUpWind
            auxregDownWind[zone] = zonalregDownWind
            auxregUpSolar[zone] = zonalregUpSolar
            auxregDownSolar[zone] = zonalregDownSolar
            auxflexWind[zone] = zonalflexWind
            auxflexSolar[zone] = zonalflexSolar

        contResHourly[gcm] = auxcontResHourly
        regUpHourly[gcm] = auxregUpHourly
        regDownHourly[gcm] = auxregDownHourly
        flexResHourly[gcm] = auxflexResHourly
        allRes[gcm] = auxallRes
        regUpWind[gcm] = auxregUpWind
        regDownWind[gcm] = auxregDownWind
        regUpSolar[gcm] = auxregUpSolar
        regDownSolar[gcm] = auxregDownSolar
        flexWind[gcm] = auxflexWind
        flexSolar[gcm] = auxflexSolar

    # convert allRes dict to data frame for printing
    writeDictToCSV(allRes, os.path.join(resultsDir, 'reservesUC' + str(ucYear) + '.csv'))

    writeDictToCSV(contResHourly, os.path.join(resultsDir, 'reservesContUC' + str(ucYear) + '.csv'))
    writeDictToCSV(regUpHourly, os.path.join(resultsDir, 'reservesRegUpUC' + str(ucYear) + '.csv'))
    writeDictToCSV(regDownHourly, os.path.join(resultsDir, 'reservesRegDownUC' + str(ucYear) + '.csv'))
    writeDictToCSV(flexResHourly, os.path.join(resultsDir, 'reservesFlexUC' + str(ucYear) + '.csv'))

    # (hourlyWindGen, hourlySolarGen) = getAggregateSolarAndWind(windCFs, ewdIdAndCapac, solarCFs, solarFilenameAndCapac)
    # netDemand = getNetDemand(zonalDemandProfile, hourlyWindGen, hourlySolarGen)

    # write2dListToCSV([[val] for val in hourlyWindGen], os.path.join(resultsDir, 'windGenUC' + str(ucYear) + '.csv'))
    # write2dListToCSV([[val] for val in hourlySolarGen], os.path.join(resultsDir, 'solarGenUC' + str(ucYear) + '.csv'))

    updateFuelPricesExistingGens(fleetUC, ucYear, fuelPricesTimeSeries)
    combineWindAndSolarToSinglePlant(fleetUC, genparam.ipmZones, genparam.dataRoot)

    hourlyCapacsCurtailedGens = importHourlyThermalCurtailments(fleetUC, ucYear, 'UC', resultsDir, genparam,
                                                                curtailparam)

    hourlyCapacsAllGens = calculateHourlyCapacsWithCurtailments(fleetUC, hourlyCapacsCurtailedGens, ucYear)

    if genparam.calculateCO2Price:
        co2Price = convertCo2CapToPrice(fleetUC, zonalHourlyWindGen, zonalHourlySolarGen, zonalDemandProfile,
                                        currCo2Cap, scaleMWtoGW, scaleDollarsToThousands, scaleLbToShortTon,
                                        genparam.dataRoot)
    else:
        co2Price = 0

    print('CO2 price:', co2Price, '$/ton')

    write2dListToCSV([[co2Price]], os.path.join(resultsDir, 'co2PriceUC' + str(ucYear) + '.csv'))
    write2dListToCSV(fleetUC, os.path.join(resultsDir, 'genFleetUC' + str(ucYear) + '.csv'))

    # Set up results lists
    ucResultsByDay = []  # list of UC GAMS models

    (genByPlant, regUpByPlant, flexByPlant, contByPlant, turnonByPlant, turnoffByPlant, regDownByPlant,
     onOffByPlant, genToRow, hourToColPlant) = setupHourlyResultsByPlant(daysForUC, fleetUC)

    pumphydroSoc, pumphydroCharge = setupHourlyPHResults(genByPlant, fleetUC)

    (sysResults, resultToRow, hourToColSys) = setupHourlySystemResults(daysForUC)

    msAndSs = [['day', 'ms', 'ss']]  # store modelstat & solvestat from GAMS

    # *******************************************************
    # ------- FOR NOW JUST USE A SINGLE GCM ----------------
    # *******************************************************
    gcm = curtailparam.listgcms[0]

    hydroPotentials = importHydroPotentialGen(ucYear, genparam, gcm)

    for dayIdx in range(0, len(daysForUC), daysOpt):

        day = daysForUC[dayIdx]

        print('---------------------------------------------------------------')
        print('Starting UC simulation of day {0:4d}'.format(day))

        daysForUCAux = list(range(day, day + daysOpt + daysLA))

        (demandUC, hourlyWindGenUC, hourlySolarGenUC, hoursForUC) = getDemandAndREGenForUC(day, daysOpt, daysLA,
                                                                                           zonalDemandProfile[gcm],
                                                                                           zonalHourlyWindGen,
                                                                                           zonalHourlySolarGen)

        (regUpUC, regDownUC, flexUC, contUC) = getResForUC(day, daysOpt, daysLA, regUpHourly[gcm], regDownHourly[gcm],
                                                           flexResHourly[gcm], contResHourly[gcm])

        hourlyCapacsUC = getHourlyCapacitiesForDays(fleetUC, hourlyCapacsAllGens[gcm], hoursForUC)

        hydroPotentialUC = getDailyHydroPotentialsUC(fleetUC, hydroPotentials, daysForUCAux, ucYear)

        if daysForUC[0] == day:  # first day, therefore no initial conditions defined. MW energy values
            (onOffInitial, genAboveMinInitial, mdtCarriedInitial, initSocDict) = setInitCondsFirstUC(fleetUC, genparam)
        else:  # other days. MW energy values
            (onOffInitial, genAboveMinInitial, mdtCarriedInitial, initSocDict) = setInitCondsPerPriorUC(ucModel, fleetUC,
                                                                                                        hoursForUC, daysOpt, daysLA,
                                                                                                        scaleMWtoGW)
        t0 = time.time()
        ucModel, ms, ss = callUnitCommitment(fleetUC, hourlyCapacsUC, hourlyWindGenUC, hourlySolarGenUC,
                                             hydroPotentialUC, demandUC, hoursForUC, onOffInitial, genAboveMinInitial,
                                             mdtCarriedInitial, initSocDict, co2Price, regUpUC, regDownUC, flexUC,
                                             contUC, genparam, reserveparam)

        print('Time (secs) for UC day ' + str(day) + ': ' + str(time.time() - t0))
        print('---------------------------------------------------------------')
        print()

        # just saves GAMS model of current day to list
        ucResultsByDay.append((day, ucModel))

        saveHourlyResultsByPlant(genByPlant, regUpByPlant, regDownByPlant, flexByPlant, contByPlant, turnonByPlant,
                                 turnoffByPlant, onOffByPlant, genToRow, hourToColPlant, ucModel, day, daysOpt)

        saveHourlyPumpedHydroResults(pumphydroSoc, pumphydroCharge, ucModel, day, daysOpt)

        saveHourlySystemResults(sysResults, resultToRow, hourToColSys, ucModel, day, daysOpt)
        msAndSs.append([day, ms, ss])

        # check if file with flag to save GDX files exists in GAMS folder
        gamsFileDir = os.path.join(genparam.resultsDir, 'GAMS')
        if os.path.exists(os.path.join(gamsFileDir, 'saveGDX.flag')):

            # read flag from file
            with open(os.path.join(gamsFileDir, 'saveGDX.flag'), 'r') as flagFile:
                saveGDX = (flagFile.readline().lower().strip() == 'true')

            if saveGDX:
                # copy input and output GDX files to output folder and rename them
                try:
                    shutil.copy(os.path.join(gamsFileDir, '_gams_py_gdb0.gdx'),
                                os.path.join(resultsDir, 'gdxUCInYear{0:4d}Day{1:04d}.gdx'.format(ucYear, day)))
                    shutil.copy(os.path.join(gamsFileDir, '_gams_py_gdb1.gdx'),
                                os.path.join(resultsDir, 'gdxUCOutYear{0:4d}Day{1:04d}.gdx'.format(ucYear, day)))
                except IOError as e:
                    print("Unable to copy file. %s" % e)
                except:
                    print("Unexpected error:", sys.exc_info())

        if day % 5 == 0:
            # write files with current results every 5 days
            writeHourlyResultsByPlant(genByPlant, regUpByPlant, regDownByPlant, flexByPlant, contByPlant,
                                      turnonByPlant, turnoffByPlant, onOffByPlant, resultsDir, ucYear, 'UC', 'Plant')
            writeHourlyStoResults(pumphydroCharge, pumphydroSoc, resultsDir, ucYear)

            write2dListToCSV(sysResults, os.path.join(resultsDir, 'systemResultsUC' + str(ucYear) + '.csv'))
            write2dListToCSV(msAndSs, os.path.join(resultsDir, 'msAndSsUC' + str(ucYear) + '.csv'))

            # save results of current UC model object to GDX file to allow "Cold start" in case of error
            write2dListToCSV([[day]], os.path.join(resultsDir, 'lastdayUC' + str(ucYear) + '.csv'))
            ucModel.get_out_db().export(os.path.join(resultsDir, "last_ucmodel.gdx"))

    # write final total results
    writeHourlyResultsByPlant(genByPlant, regUpByPlant, regDownByPlant, flexByPlant, contByPlant,
                              turnonByPlant, turnoffByPlant, onOffByPlant, resultsDir, ucYear, 'UC', 'Plant')
    writeHourlyStoResults(pumphydroCharge, pumphydroSoc, resultsDir, ucYear)
    write2dListToCSV(sysResults, os.path.join(resultsDir, 'systemResultsUC' + str(ucYear) + '.csv'))
    write2dListToCSV(msAndSs, os.path.join(resultsDir, 'msAndSsUC' + str(ucYear) + '.csv'))

    return ucResultsByDay, genByPlant


def callUnitCommitment(fleetUC, hourlyCapacsUC, hourlyWindGenUC, hourlySolarGenUC, hydroPotentialUC, demandUC,
                       hoursForUC, onOffInitial, genAboveMinInitial, mdtCarriedInitial, initSoc, co2Price, regUpUC,
                       regDownUC, flexUC, contUC, genparam, reserveparam):
    """RUN UNIT COMMITMENT MODEL

    :param fleetUC:
    :param hourlyCapacsUC:
    :param hourlyWindGenUC:
    :param hourlySolarGenUC:
    :param demandUC:
    :param hoursForUC:
    :param onOffInitial:
    :param genAboveMinInitial:
    :param mdtCarriedInitial:
    :param co2Price:
    :param regUpUC:
    :param regDownUC:
    :param flexUC:
    :param contUC:
    :param genparam:
    :param reserveparam:
    :return:

    """
    currDir = os.getcwd()

    (ucFilename, dataRoot, scaleMWtoGW, scaleDollarsToThousands,
     scaleLbToShortTon) = (genparam.ucFilename, genparam.dataRoot, genparam.scaleMWtoGW,
                           genparam.scaleDollarsToThousands, genparam.scaleLbToShortTon)

    gamsFileDir = os.path.join(genparam.resultsDir, 'GAMS')
    gamsSysDir = genparam.pathSysGams.strip()

    if not gamsSysDir == '':
        wsUC = GamsWorkspace(working_directory=gamsFileDir, system_directory=gamsSysDir)
    else:
        wsUC = GamsWorkspace(working_directory=gamsFileDir)

    dbUC = wsUC.add_database()

    # Add sets and parameters to database
    cnse = 10000

    (genSet, genSymbols, hydroGenSet, hydroGenSymbols, hourSet, hourSymbols, zoneSet, zoneSymbols,
     lineSet) = addSetsToDatabaseUC(dbUC, fleetUC, hoursForUC, genparam)

    addParametersToDatabaseUC(dbUC, hourlyCapacsUC, hourlyWindGenUC, hourlySolarGenUC, hydroPotentialUC, demandUC,
                              fleetUC, genSet, genSymbols, hydroGenSet, hydroGenSymbols, hourSet, hourSymbols, zoneSet,
                              zoneSymbols, lineSet, cnse, onOffInitial, genAboveMinInitial, mdtCarriedInitial, initSoc,
                              co2Price, hoursForUC, flexUC, contUC, regUpUC, regDownUC, genparam, reserveparam)

    # Load and run GAMS model
    ucModel = wsUC.add_job_from_file(ucFilename)
    optUC = GamsOptions(wsUC)
    optUC.defines['gdxincname'] = dbUC.name
    ucModel.run(optUC, databases=dbUC, output=sys.stdout)
    ms, ss = ucModel.out_db['pModelstat'].find_record().value, ucModel.out_db['pSolvestat'].find_record().value
    if int(ms) != 8 or int(ss) != 1: print('Modelstat & solvestat:', ms, ' & ', ss, ' (should be 8 and 1)')
    return ucModel, ms, ss


def addSetsToDatabaseUC(db, fleetUC, hoursForUC, genparam):
    """ADD SETS

    :param db:
    :param fleetUC:
    :param hoursForUC:
    :return:
    """

    ipmZoneNums, lineList, ptCurtailedAll = genparam.ipmZoneNums, genparam.lineList, genparam.ptCurtailedAll

    zoneSet, zoneSymbols = addZoneSets(db, ipmZoneNums)  # create string for set IDs

    lineSet = addLineSets(db, lineList)

    (hourSet, hourSymbols) = addHourSet(db, hoursForUC)

    (genSet, genSymbols, hydroGenSet, hydroGenSymbols, pumpHydroGenSet,
     pumpHydroGenSymbols) = addGeneratorSets(db, fleetUC)

    return genSet, genSymbols, hydroGenSet, hydroGenSymbols, hourSet, hourSymbols, zoneSet, zoneSymbols, lineSet


def addParametersToDatabaseUC(db, hourlyCapacsUC, hourlyWindGenUC, hourlySolarGenUC, hydroPotentialUC, demandUC,
                              fleetUC, genSet, genSymbols, hydroGenSet, hydroGenSymbols, hourSet, hourSymbols, zoneSet,
                              zoneSymbols, lineSet, cnse, onOffInitial, genAboveMinInitial, mdtCarriedInitial, initSoc,
                              co2Price, hoursForUC, flexUC, contUC, regUpUC, regDownUC, genparam, reserveparam):
    """ADD PARAMETERS

    :param db:
    :param hourlyCapacsUC:
    :param hourlyWindGenUC:
    :param hourlySolarGenUC:
    :param hydroPotentialUC:
    :param demandUC:
    :param fleetUC:
    :param genSet:
    :param genSymbols:
    :param hydroGenSet:
    :param hydroGenSymbols:
    :param hourSet:
    :param hourSymbols:
    :param zoneSet:
    :param zoneSymbols:
    :param lineSet:
    :param cnse:
    :param onOffInitial:
    :param genAboveMinInitial:
    :param mdtCarriedInitial:
    :param initSoc: initial state of charge for pumped storage
    :param co2Price:
    :param hoursForUC:
    :param flexUC:
    :param contUC:
    :param regUpUC:
    :param regDownUC:
    :param genparam: general parameters
    :param reserveparam: reserve parameters


    
    """

    (rrToRegTime, rrToFlexTime, rrToContTime) = (reserveparam.rampRateToRegReserveScalar,
                                                 reserveparam.rampRateToFlexReserveScalar,
                                                 reserveparam.rampRateToContReserveScalar)

    # db, demandCEZonal, zoneSet, hourSet, gcmSet, hoursForCE, ipmZones, ipmZoneNums, scaleMWtoGW
    addDemandParam(db, demandUC, zoneSet, hourSet, None, hoursForUC, genparam.ipmZones, genparam.ipmZoneNums,
                   genparam.scaleMWtoGW)

    addEguParams(db, fleetUC, genSet, genSymbols, genparam.ipmZones, genparam.ipmZoneNums, genparam.scaleLbToShortTon,
                 genparam.scaleMWtoGW)

    addEguHourlyParams(db, hourlyCapacsUC, None, genSet, hourSet, hoursForUC, genparam.scaleMWtoGW)

    addEguOpCostParam(db, fleetUC, genSet, genSymbols, genparam.scaleLbToShortTon, genparam.scaleMWtoGW,
                      genparam.scaleDollarsToThousands, co2Price)

    addEguUCParams(db, fleetUC, genSet, genSymbols, genparam.scaleMWtoGW, genparam.scaleDollarsToThousands)
    addEguInitialConditions(db, genSet, genSymbols, fleetUC, onOffInitial, genAboveMinInitial, mdtCarriedInitial,
                            genparam.scaleMWtoGW)

    # TODO: Check this with Michael
    stoMarket = []
    addEguEligibleToProvideRes(db, fleetUC, genSet, genSymbols, stoMarket)

    addExistingRenewableMaxGenParams(db, None, zoneSet, genparam.ipmZones, genparam.ipmZoneNums, hourSet, hoursForUC,
                                     hourlySolarGenUC, hourlyWindGenUC, genparam.scaleMWtoGW)

    addHydroMaxGenUC(db, hydroGenSet, hydroPotentialUC, genparam.scaleMWtoGW)

    addRegReserveParameters(db, regUpUC, regDownUC, rrToRegTime, hourSet, hoursForUC, zoneSet, 'UC', genparam)

    # TODO: What was this function used for??? Check with Michael (so far i am ignoring it)
    # addEguRegCostParam(db, fleetUC, genSet, genSymbols, genparam.scaleMWtoGW, genparam.scaleDollarsToThousands)

    addFlexReserveParameters(db, flexUC, rrToFlexTime, hourSet, hoursForUC, zoneSet, 'UC', genparam)
    addContReserveParameters(db, contUC, rrToContTime, hourSet, hoursForUC, zoneSet, genparam)

    addCostNonservedEnergy(db, cnse, genparam.scaleMWtoGW, genparam.scaleDollarsToThousands)

    addCo2Price(db, co2Price, genparam.scaleDollarsToThousands)

    # Add transmission line constraints
    addLineCapacs(db, genparam.lineCapacs, lineSet, genparam.lineList, genparam.scaleMWtoGW)
    addLineSourceAndSink(db, lineSet, genparam.lineList, genparam.ipmZones, genparam.ipmZoneNums)

    # add pumped hydro parameters to database
    pumpHydroGenSet = db.get_set('pumphydroegu')
    pumpHydroGenSymbols = [symb.keys[0] for symb in pumpHydroGenSet.__iter__()]

    effDict = dict()
    for symb in pumpHydroGenSymbols: effDict[symb] = genparam.phEff
    (effname, effdesc) = ('pEfficiency', 'efficiency')
    effParam = add1dParam(db, effDict, pumpHydroGenSet, pumpHydroGenSymbols, effname, effdesc)

    socDict = getSocDict(fleetUC, genparam.phMaxSoc, pumpHydroGenSymbols, genparam.scaleMWtoGW)
    (socname, socdesc) = ('pMaxsoc', 'max state of charge (GWh)')
    socparam = add1dParam(db, socDict, pumpHydroGenSet, pumpHydroGenSymbols, socname, socdesc)

    (socname, socdesc) = ('pInitsoc', 'initial state of charge (GWh)')
    initsocparam = add1dParam(db, initSoc, pumpHydroGenSet, pumpHydroGenSymbols, socname, socdesc)


def create_description_file(genparam, curtailparam):
    """ Create a description file of the case in the output folder

    """

    outstr = '\n\n'
    outstr = outstr + '------------------ CE/UC simulation ------------------\n'
    outstr = outstr + '\n\n'
    outstr = outstr + 'Machine info         : {}\n'.format(platform.platform())
    outstr = outstr + 'Date & time          : {}\n'.format(time.strftime('%Y-%m-%d  %H:%M %Z'))
    outstr = outstr + 'CE run               : {}\n'.format(genparam.runCE)
    outstr = outstr + 'UC run               : {}\n'.format(genparam.runUC)
    outstr = outstr + 'Output Folder        : {}\n'.format(genparam.resultsDir)
    outstr = outstr + 'Input Folder         : {}\n'.format(genparam.dataRoot)
    outstr = outstr + 'Water data Folder    : {}\n'.format(curtailparam.rbmRootDir)

    outstr = outstr + '\n'

    outstr = outstr + 'Initial Year         : {}\n'.format(genparam.startYear)
    outstr = outstr + 'End Year             : {}\n'.format(genparam.endYear)
    outstr = outstr + 'Step CE              : {} years\n'.format(genparam.yearStepCE)

    outstr = outstr + '\n'

    outstr = outstr + 'Days per season      : {}\n'.format(genparam.daysPerSeason)
    outstr = outstr + 'Analysis area        : {}\n'.format(genparam.analysisArea)

    if genparam.referenceCase:
        outstr = outstr + 'GCMs considered      : {}\n'.format("Baseline case (no climate change)")
    else:
        outstr = outstr + 'GCMs considered      : {}\n'.format("{0}".format(", ".join(
            str(i) for i in curtailparam.listgcms)))

    outstr = outstr + 'CO2 Cap case         : {}\n'.format(genparam.co2CapScenario)

    outstr = outstr + '\n\n'

    outstr = outstr + 'Transmission lines capacities:\n'

    for k, value in genparam.lineCapacs.items():
        outstr = outstr + '{0:<21}: {1: .2f} MW\n'.format(k, value)

    outstr = outstr + '\n\n'

    fname = os.path.join(os.path.expanduser(genparam.resultsDir), 'DESCRIPTION.TXT')

    with open(fname, 'w') as f:
        f.write(outstr)


########### TESTING PURPOSES ONLY ####################################
# REMOVE THIS FUNCTION - JUST INCLUDED FOR TESTING PURPOSES TO EXPEDITE UC
def removeHydroForTesting(fleetUC):
    fuelTypeCol = fleetUC[0].index('Modeled Fuels')
    idxs = []
    for idx in range(len(fleetUC) - 1, 0, -1):
        if fleetUC[idx][fuelTypeCol] == 'Hydro': fleetUC.pop(idx)
