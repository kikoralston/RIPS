# Michael Craig, 31 May 2016
# Run capacity expansion and dispatch models w/ climate impacts.

import numpy as np
import shutil
import gc
import platform

try:
    from gams import *
except ImportError:
    print('gams module not found! GAMS functions will not work.')

from GAMSUtil import *
from SetupFleet import *
from PrepareData import *
from thermalderatings import *
from ProcessResults import *
from co2cap import *
from renewables import *
from ipmzones import *

from demand.ForecastDemandWithRegression import forecastZonalDemandWithReg
from reserves.ReservesWWSIS import calcWWSISReserves

sys.stdout.flush()


def masterFunction(genparam, reserveparam, curtailparam):
    """Main function of CE/UCED simulation

    This is the main function that must be called to start a CE/UCED simulation. It takes three arguments that define
    the parameters of the simulation. Each argument is an object with a set of specific parameters for the simulation.
    These objects must be created and the fields defined before being passed to this function.

    :param genparam: object of type :mod:`Generalparameters`
    :param reserveparam: object of type :mod:`Reserveparameters`
    :param curtailparam: object of type :mod:`Curtailmentparameters`
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

    (capacExpModelsEachYear, capacExpBuilds, capacExpGenByGens, capacExpRetiredUnitsByCE,
     capacExpRetiredUnitsByAge) = ([], [['TechnologyType']], [['ORIS+UnitID']], [], [])

    genFleetNoRetiredUnits, genFleetPriorCE, priorCEout_db, priorHoursCE = None, None, None, None

    t_start = time.time()
    # begin loop
    for currYear in range(genparam.startYear, genparam.endYear, genparam.yearStepCE):

        # get new curtparam only with subset of GCMs
        curtparam_year = choose_gcms_for_ce(currYear, genparam, curtailparam)

        t_year = time.time()
        print('\n----------------------------------------------------------------\n')
        print('Starting loop for year {0:4d}\n'.format(currYear))
        currCo2Cap = interpolateCO2Cap(currYear, genparam)

        print('CO2 cap in year {1:4d}: {0:,.3f} million tons'.format(currCo2Cap / 1e6, currYear))

        zonalDemandProfile, zonalTempDfs = forecastZonalDemandWithReg(currYear, genparam, curtparam_year)

        if genparam.runCE:
            # run capacity expansion model

            get_init_conditions_ce(currYear, genparam, genFleet, genFleetNoRetiredUnits, genFleetPriorCE,
                                   priorCEout_db, capacExpBuilds, capacExpGenByGens, capacExpRetiredUnitsByCE,
                                   capacExpRetiredUnitsByAge, priorHoursCE)

            # only run model in genparam.startYear if it is a cold start
            if genparam.hotStart or currYear > genparam.startYear:
                (genFleet, genFleetNoRetiredUnits, genFleetPriorCE, priorCEout_db,
                 priorHoursCE) = runCapacityExpansion(genFleet, zonalDemandProfile, currYear, currCo2Cap,
                                                      capacExpModelsEachYear, capacExpBuilds, capacExpGenByGens,
                                                      capacExpRetiredUnitsByCE, capacExpRetiredUnitsByAge,
                                                      genFleetPriorCE, priorCEout_db, priorHoursCE,
                                                      genparam, reserveparam, curtparam_year)

                # write hours for CE to pickle file
                with open(os.path.join(genparam.resultsDir, 'CE', 'hoursCE_{0}.pkl'.format(currYear)), 'wb') as f:
                    pk.dump(priorHoursCE, f, pk.HIGHEST_PROTOCOL)

            # since model was executed successfully, set hotStart to False before running next year
            genparam.hotStart = False

        if genparam.runUC:
            # Either only runs 2015, or runs in all but 2015
            if ((currYear == genparam.startYear and genparam.runFirstUCYear) or
                    (currYear > genparam.startYear and not genparam.runFirstUCYear)):
                genFleetNoRetiredUnits = loadCEFleet(currYear, genparam.resultsDir)

                uc_out = runUnitCommitment(genFleetNoRetiredUnits, zonalDemandProfile, currYear, currCo2Cap,
                                           genparam, reserveparam, curtparam_year)

        print()
        print('Elapsed Time: ' + str_elapsedtime(t_year))

        print('\n----------------------------------------------------------------')

    print()
    print('************ Total Elapsed Time: ' + str_elapsedtime(t_start) + '************')


def getInitialFleetAndDemand(genparam, reserveparam):
    """Set up initial generator fleet

    Reads folders containing data of existing power plant data and compiles 2d list with generator fleet data.

    :param genparam: object of class :mod:`Generalparameters`
    :param reserveparam: object of class :mod:`Reserveparameters`
    :return: 2d list with initial generator fleet data
    """

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
    """Run one annual iteration of the capacity expansion simulation

    This function does all the reading and preprocessing before executing capacity expansion optimization

    (capacExpModelsEachYear, capacExpBuilds, capacExpGenByGens, capacExpRetiredUnitsByCE,
     capacExpRetiredUnitsByAge) = ([], [['TechnologyType']], [['ORIS+UnitID']], [], [])

    :param genFleet: (2d list) with original existing generator fleet (in the first year of the simulation)
    :param zonalDemandProfile: (dict) nested dictionary with hourly demand by ipm zone in each gcm case being simulated.
                               Dictionary is {gcm: {zone: [hourly demand]}}
    :param currYear: (int) year of simulation
    :param currCo2Cap: float with value in MMton/year of CO2 emission limit
    :param capacExpModelsEachYear: (list) list with gams models in previous CE simulation years
    :param capacExpBuilds: (2d list) list with complete builds (by plant type) of capacity expansion model
    :param capacExpGenByGens: (2d list) list with generation of each individual plant in previous simulation year
                              (used to compute retirements)
    :param capacExpRetiredUnitsByCE: (list) list with retired units because of low CF
    :param capacExpRetiredUnitsByAge: (list) list with retired units because of age
    :param genFleetPriorCE: 2d list with existing generator fleet in current year of simulation
                            (including prior additions)
    :param priorCEout_db: gams data base object with results of capacity expansion iteration in previous simulation year
    :param priorHoursCE: (dict) dictionary with hours considered in previous CE simulation year
    :param genparam: object of class :mod:`Generalparameters`
    :param reserveparam: object of class :mod:`Reserveparameters`
    :param curtailparam: object of class :mod:`Curtailmentparameters`
    :return: resulting generator fleet (including retired units), resulting generator fleet (excluding retired units),
             generator fleet before CE expansion decisions, gams data base with results, dictionary with hours of the
             year used in CE simulation.
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

    genFleetForCE = createFleetForCurrentCELoop(genFleet, currYear, capacExpRetiredUnitsByAge, genparam)  # removes all retired units

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

        zonalGenFleet = [genFleetForCE[0]] + [row for row in genFleetForCE if row[zoneCol] == zone]

        capacCol = zonalGenFleet[0].index('Capacity (MW)')
        plantTypeCol = zonalGenFleet[0].index('PlantType')

        windCapacInZone = sum([float(row[capacCol]) for row in zonalGenFleet[1:] if row[plantTypeCol] == 'Wind'])
        solarCapacInZone = sum([float(row[capacCol]) for row in zonalGenFleet[1:] if row[plantTypeCol] == 'Solar PV'])

        print('Existing Wind: ')
        windCFs = getRenewableCFData(currZone=zone, genparam=genparam, fleetCap=50, capacInCurrFleet=windCapacInZone,
                                     type='wind', existing=True)[0]

        print('Existing Solar: ')
        solarCFs = getRenewableCFData(currZone=zone, genparam=genparam, fleetCap=50, capacInCurrFleet=solarCapacInZone,
                                      type='solar', existing=True)[0]

        # For EXISTING renewables, store GENERATION (MWh) for each zone in dictionaries
        # getRenewableCFData() returns data frames, convert to 1-d list
        zonalHourlyWindGen[zone], zonalHourlySolarGen[zone] = windCFs['gen'].tolist(), solarCFs['gen'].tolist()

        print('Got Existing RE CFs.')

        print('New Wind: ')
        newWindCFs = getRenewableCFData(currZone=zone, genparam=genparam, fleetCap=50, capacInCurrFleet=windCapacInZone,
                                        type='wind', existing=False)[0]

        print('New Solar: ')
        newSolarCFs = getRenewableCFData(currZone=zone, genparam=genparam, fleetCap=50, capacInCurrFleet=solarCapacInZone,
                                         type='solar', existing=False)[0]

        # For NEW renewables, store 1d list with CFs for each zone in dictionaries
        # getRenewableCFData() returns data frames, convert to dictionaries with segment keys and 1-d list for values
        zonalNewWindCFs[zone] = {}
        for s in newWindCFs['segment'].unique():
            zonalNewWindCFs[zone][s] = list(newWindCFs[newWindCFs['segment'] == s]['cf'].values)

        zonalNewSolarCFs[zone] = {}
        for s in newSolarCFs['segment'].unique():
            zonalNewSolarCFs[zone][s] = list(newSolarCFs[newSolarCFs['segment'] == s]['cf'].values)

        print('Got new RE CFs.')

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
    hrlyCurtailmentsAllTechsInTgtYr = determineHrlyCurtailmentsForNewTechs(eligibleCellWaterTs, newTechsCE, currYear,
                                                                           genparam, curtailparam)

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

    (capacExpModel, ms, ss) = callCapacityExpansion(genFleetForCE, hourlyCapacsCE, hourlyCurtailedTechCapacsCE,
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
    genFleetNoRetiredUnits = createFleetForCurrentCELoop(genFleet, currYear, [], genparam)

    writeCEInfoToCSVs(capacExpBuilds, capacExpGenByGens, capacExpRetiredUnitsByCE,
                      capacExpRetiredUnitsByAge, resultsDir, currYear)
    write2dListToCSV(genFleet, os.path.join(resultsDir, 'genFleetAfterCE' + str(currYear) + '.csv'))

    # Write gen fleet for UC to special folder for ease of transfer
    ceUCDir = os.path.join(genparam.resultsDir, 'CEtoUC')

    if not os.path.exists(ceUCDir):
        os.makedirs(ceUCDir)

    write2dListToCSV(genFleetNoRetiredUnits, os.path.join(ceUCDir, 'genFleetCEtoUC' + str(currYear) + '.csv'))

    return genFleet, genFleetNoRetiredUnits, genFleetForCE, capacExpModel.out_db, hoursForCE


def importHourlyThermalCurtailments(genFleet, currYear, modelName, resultsDir, genparam, curtailparam):
    """Import hourly thermal curtailments by generator

    This is a wrapper function that calls thermal curtailments simulations and cleans up the output so it
    can be used by the CE/UCED simulations

    :param genFleet: 2d list (matrix) with data of generator fleet
    :param currYear: current year of analysis
    :param modelName: string with model name ('CE' or 'UC')
    :param resultsDir: string specifying folder where output files will be saved
    :param genparam: object of type :mod:`Generalparameters`
    :param curtailparam: object of type :mod:`Curtailmentparameters`
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
    """Call capacity expansion optimization

    This funtion loads GAMS workspace and runs capacity expansion optimization

    :param genFleetForCE: 2d list (matrix) with data of generator fleet
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
    :param genparam: object of type :mod:`Generalparameters`
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

    # create lists with names of types of winds and solar (aggregated in block according to CF)
    # (these will be needed to create GAMS sets)
    #
    # pick name of first GCM (needed to get list of keys)
    g = list(newWindCFsCEZonal.keys())[0]

    # get longest list of types of wind and sort
    blocksWind = [[b for b in newWindCFsCEZonal[g][z]] for z in newWindCFsCEZonal[g]]
    blocksWind = [b for b in blocksWind if len(b) == max([len(b) for b in blocksWind])][0]
    blocksWind.sort()

    # get longest list of types of solar and sort
    blocksSolar = [[b for b in newSolarCFsCEZonal[g][z]] for z in newSolarCFsCEZonal[g]]
    blocksSolar = [b for b in blocksSolar if len(b) == max([len(b) for b in blocksSolar])][0]
    blocksSolar.sort()

    # Add sets and parameters to database
    x = addSetsToDatabaseCE(db, genFleetForCE, hoursForCE, newTechsCE, repHrsBySeason, specialHrs,
                            peakDemandHourZonal, cellsToZones, genparam, blocksWind=blocksWind,
                            blocksSolar=blocksSolar)

    addParametersToDatabaseCE(db, genparam, hoursForCE, hourlyCapacsCE, hourlyWindGenCEZonal, hourlySolarGenCEZonal,
                              demandCEZonal, newTechsCE, genFleetForCE, planningReserveZonal, newWindCFsCEZonal,
                              newSolarCFsCEZonal, currCo2Cap, seasonDemandWeights, hourlyCurtailedTechCapacsCE,
                              cellsToZones, hydroPotPerSeason)

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
    file.write("parallelmode -1\n")
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

    return capacExpModel, ms, ss


def addSetsToDatabaseCE(db, genFleetForCE, hoursForCE, newTechsCE, repHrsBySeason, specialHrs, peakDemandHourZonal,
                        cellsToZones, genparam, blocksWind=None, blocksSolar=None):
    """Add sets needed to GAMS database

    This function adds the sets to the GAMS database that are needed to run the capacity expansion optimization

    :param db: gams database
    :param genFleetForCE: 2d list (matrix) with data of generator fleet
    :param hoursForCE:
    :param newTechsCE:
    :param repHrsBySeason:
    :param specialHrs:
    :param peakDemandHourZonal:
    :param cellsToZones:
    :param genparam: object of type :mod:`Generalparameters`
    :param blocksWind:
    :param blocksSolar:
    :return:
    """
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
     techNotCurtailedSet, techNotCurtailedSymbols) = addNewTechsSets(db, newTechsCE, ptCurtailedAll,
                                                                     blocksWind=blocksWind, blocksSolar=blocksSolar)

    # create additional set with type of plant (no cooling info)
    # typeSymbols = set(map(lambda x: x.split('+')[0], techSymbols))
    # (typeSetName, typeSetDescrip, typeSetDim) = ('type', 'set with types of new plants (no cooling info)', 1)
    # typeSet = addSet(db, typeSymbols, typeSetName, typeSetDescrip, typeSetDim)

    # create additional 2d set mapping type of plant to plant+cooling
    # typeList = list(map(lambda x: x.split('+')[0], techSymbols))
    # (tech2dSetName, tech2dSetDescrip, tech2dSetDim) = ('tech2d', '2d set mapping types to techs', 2)
    # techSet2 = db.add_set(tech2dSetName, tech2dSetDim, tech2dSetDescrip)
    # for t1 in typeSymbols:
    #     idx = [i for i, x in enumerate(typeList) if x == t1]
    #     list_techs = [techSymbols[i] for i in idx]
    #    for t2 in list_techs:
    #        techSet2.add_record([t1, t2])

    peakHourSet, peakHrSymbols = addPeakHourSubset(db, peakDemandHourZonal, genparam)

    return (genSet, genSymbols, hourSet, hourSymbols, techSet, techSymbols, techCurtailedSet,
            techCurtailedSymbols, renewTechSet, renewTechSymbols, techNotCurtailedSymbols,
            hydroGenSet, hydroGenSymbols, zoneSet, zoneSymbols, lineSet, cellSet, peakHourSet,
            peakHrSymbols, pumpHydroGenSet, pumpHydroGenSymbols, gcmSet)


def addParametersToDatabaseCE(db, genparam, hoursForCE, hourlyCapacsCE, hourlyWindGenCEZonal, hourlySolarGenCEZonal,
                              demandCEZonal, newTechsCE, genFleetForCE, planningReserveZonal,  newWindCFsCEZonal,
                              newSolarCFsCEZonal, currCo2Cap, seasonDemandWeights, hourlyCurtailedTechCapacsCE,
                              cellsToZones, hydroPotPerSeason):
    """Add parameters needed to GAMS database

    This function adds the parameters to the GAMS database that are needed to run the capacity expansion optimization

    :param db: gams database object
    :param genparam: object of type :mod:`Generalparameters`
    :param hoursForCE: (dict) hours of the year considered in CE simulation
    :param hourlyCapacsCE:
    :param hourlyWindGenCEZonal:
    :param hourlySolarGenCEZonal:
    :param demandCEZonal:
    :param newTechsCE:
    :param genFleetForCE:
    :param planningReserveZonal:
    :param newWindCFsCEZonal:
    :param newSolarCFsCEZonal:
    :param currCo2Cap: (float) value of co2 cap
    :param seasonDemandWeights:
    :param hourlyCurtailedTechCapacsCE:
    :param cellsToZones:
    :param hydroPotPerSeason:
    """
    # update newTechsCE to account for blocks of Wind and Solar with different CF
    newTechsCEcopy = copy.deepcopy(newTechsCE)
    techCol = newTechsCEcopy[0].index('TechnologyType')

    renewTechSet = db.get_set('techrenew')
    renewTechSymbols = [symb.keys[0] for symb in renewTechSet.__iter__()]

    # wind
    windTechs = [t for t in renewTechSymbols if 'Wind' in t]
    windRow = [row for row in newTechsCEcopy if row[techCol] == 'Wind'][0]
    windRowNumber = newTechsCEcopy.index(windRow)
    x = newTechsCEcopy.pop(windRowNumber)

    for t in windTechs:
        row = copy.deepcopy(windRow)
        row[techCol] = t
        newTechsCEcopy.append(row)

    # solar
    solarTechs = [t for t in renewTechSymbols if 'Solar PV' in t]
    solarRow = [row for row in newTechsCEcopy if row[techCol] == 'Solar PV'][0]
    solarRowNumber = newTechsCEcopy.index(solarRow)
    x = newTechsCEcopy.pop(solarRowNumber)

    for t in solarTechs:
        row = copy.deepcopy(solarRow)
        row[techCol] = t
        newTechsCEcopy.append(row)

    addTechParams(db, newTechsCEcopy, genparam.scaleMWtoGW, genparam.scaleDollarsToThousands,
                  genparam.scaleLbToShortTon, genparam.ptCurtailedAll)

    genSet = db.get_set('egu')
    genSymbols = [symb.keys[0] for symb in genSet.__iter__()]
    addEguParams(db, genFleetForCE, genSet, genSymbols, genparam.ipmZones, genparam.ipmZoneNums,
                 genparam.scaleLbToShortTon, genparam.scaleMWtoGW)

    pumpHydroGenSet = db.get_set('pumphydroegu')
    pumpHydroGenSymbols = [symb.keys[0] for symb in pumpHydroGenSet.__iter__()]
    addPumpHydroParams(db, genFleetForCE, genparam.phEff, genparam.phMaxSoc, genparam.phInitSoc, pumpHydroGenSet,
                       pumpHydroGenSymbols, genparam.scaleMWtoGW)

    addDiscountRateParam(db, genparam.discountRate)
    addCppEmissionsCap(db, currCo2Cap)
    addSeasonDemandWeights(db, seasonDemandWeights)

#    addMaxNumNewBuilds(db, newTechsCE, zoneSet, ipmZones, ipmZoneNums, typeSet, maxAddedZonalCapacPerTech,
#                       ptCurtailedAll)

    addMaxNumNewRenew(db, newWindCFsCEZonal, newSolarCFsCEZonal, genparam.ipmZones, genparam.ipmZoneNums)

    # Add hydro max gen per time block
    hydroGenSet = db.get_set('hydroegu')
    gcmSet = db.get_set('g')
    addHydroMaxGenPerSeason(db, hydroGenSet, gcmSet, hydroPotPerSeason, genparam.scaleMWtoGW)

    # Add zonal demand, planning margin, and RE generation
    zoneSet = db.get_set('z')
    zoneSymbols = [symb.keys[0] for symb in zoneSet.__iter__()]
    hourSet = db.get_set('h')
    addDemandParam(db, demandCEZonal, zoneSet, hourSet, gcmSet, hoursForCE, genparam.ipmZones, genparam.ipmZoneNums,
                   genparam.scaleMWtoGW)

    addPlanningReserveParam(db, planningReserveZonal, genparam.ipmZones, genparam.ipmZoneNums, zoneSet, zoneSymbols,
                            genparam.scaleMWtoGW)
    #    addPeakHourToZoneParam(db, peakDemandHourZonal, peakHourSet, peakHrSymbols, ipmZones, ipmZoneNums)

    addExistingRenewableMaxGenParams(db, gcmSet, zoneSet, genparam.ipmZones, genparam.ipmZoneNums, hourSet, hoursForCE,
                                     hourlySolarGenCEZonal, hourlyWindGenCEZonal, genparam.scaleMWtoGW)

    addRenewTechCFParams(db, renewTechSet, renewTechSymbols, gcmSet, zoneSet, hourSet, hoursForCE,
                         newWindCFsCEZonal, newSolarCFsCEZonal, genparam.ipmZones, genparam.ipmZoneNums)

    # Add hourly capacity for curtailed techs
    cellSet = db.get_set('c')
    techCurtailedSet = db.get_set('techcurtailed')
    addTechCurtailedHourlyCapac(db, hourlyCurtailedTechCapacsCE, gcmSet, cellSet, techCurtailedSet, hourSet, hoursForCE,
                                genparam.scaleMWtoGW)

    # Add tech heat rate and op cost params that reflect degradation
    addEguHourlyParams(db, hourlyCapacsCE, gcmSet, genSet, hourSet, hoursForCE, genparam.scaleMWtoGW)

    addEguOpCostParam(db, genFleetForCE, genSet, genSymbols, genparam.scaleLbToShortTon, genparam.scaleMWtoGW,
                      genparam.scaleDollarsToThousands)

    # Add zone & line parameters
    lineSet = db.get_set('l')
    addLineCapacs(db, genparam.lineCapacs, lineSet, genparam.lineList, genparam.scaleMWtoGW)
    addLineSourceAndSink(db, lineSet, genparam.lineList, genparam.ipmZones, genparam.ipmZoneNums)
    addCellsToZones(db, cellSet, cellsToZones, genparam.ipmZones, genparam.ipmZoneNums)


def runUnitCommitment(genFleet, zonalDemandProfile, ucYear, currCo2Cap, genparam, reserveparam, curtailparam):
    """Run unit commitment simulation for all climate simulations in parallel for the complete simulation period

    :param genFleet: 2d list with generator fleet for UCED simulation
    :param zonalDemandProfile: (2d list) list with hourly demand in each zone
    :param ucYear: integer with year of UCED run
    :param currCo2Cap: (float) value of CO2 cap
    :param genparam: object of class :mod:`Generalparameters`
    :param reserveparam: object of type :mod:`Reserveparameters`
    :param curtailparam: object of type :mod:`Curtailmentparameters`
    :return:
    """

    (startYear, fuelPricesTimeSeries, scaleMWtoGW, scaleDollarsToThousands, scaleLbToShortTon,
     daysOpt, daysLA, projectName, resultsDir, ocAdderMin, ocAdderMax,
     windGenDataYr) = (genparam.startYear, genparam.fuelPricesTimeSeries, genparam.scaleMWtoGW,
                       genparam.scaleDollarsToThousands, genparam.scaleLbToShortTon, genparam.daysOpt, genparam.daysLA,
                       genparam.projectName, genparam.resultsDir, genparam.ocAdderMin, genparam.ocAdderMax,
                       genparam.windGenDataYr)

    resultsDir = os.path.join(resultsDir, 'UC')

    if not os.path.exists(resultsDir):
        os.makedirs(resultsDir)

    print('Entering UC loop for year {0:4d}'.format(ucYear))

    fleetUC = copy.deepcopy(genFleet)

    (startWindCapacForCFs, startSolarCapacForCFs) = (0, 0)
    zoneCol = fleetUC[0].index('Region Name')
    zonalHourlyWindGen, zonalHourlySolarGen = dict(), dict()
    zonalHourlyDfWind, zonalHourlyDfSolar = dict(), dict()
    zonalSubhourlyDfWind, zonalSubhourlyDfSolar = dict(), dict()

    print('Computing CFs for Renewables...')
    for zone in genparam.ipmZones:
        print('Zone ', zone)

        zonalGenFleet = [fleetUC[0]] + [row for row in fleetUC if row[zoneCol] == zone]

        capacCol = zonalGenFleet[0].index('Capacity (MW)')
        plantTypeCol = zonalGenFleet[0].index('PlantType')

        windCapacInZone = sum([float(row[capacCol]) for row in zonalGenFleet[1:] if row[plantTypeCol] == 'Wind'])
        solarCapacInZone = sum([float(row[capacCol]) for row in zonalGenFleet[1:] if row[plantTypeCol] == 'Solar PV'])

        print('Existing Wind. Installed Capacity = {0:.2f} MW.'.format(windCapacInZone))
        windCfsDtHr, windCfsDtSubhr = getRenewableCFData(currZone=zone, genparam=genparam, fleetCap=50,
                                                         capacInCurrFleet=windCapacInZone, type='wind', existing=True,
                                                         subHour=True)
        print('Existing Solar. Installed Capacity = {0:.2f} MW.'.format(solarCapacInZone))
        solarCfsDtHr, solarCfsDtSubhr = getRenewableCFData(currZone=zone, genparam=genparam, fleetCap=50,
                                                           capacInCurrFleet=solarCapacInZone, type='solar',
                                                           existing=True, subHour=True)

        # get 1-d list with hourly generation of existing generators
        windCFs, solarCFs = windCfsDtHr['gen'].tolist(), solarCfsDtHr['gen'].tolist()

        # convert pandas DFs to 2-d lists with necessary columns datetime and gen
        windCfsDtHr = [['datetime', 'totalGen(MWh)']] + windCfsDtHr[['datetime', 'gen']].as_matrix().tolist()
        if windCfsDtSubhr is not None:
            windCfsDtSubhr = [['datetime', 'totalGen(MWh)']] + windCfsDtSubhr[['datetime', 'gen']].as_matrix().tolist()

        solarCfsDtHr = [['datetime', 'totalGen(MWh)']] + solarCfsDtHr[['datetime', 'gen']].as_matrix().tolist()
        if solarCfsDtSubhr is not None:
            solarCfsDtSubhr = [['datetime', 'totalGen(MWh)']] + solarCfsDtSubhr[['datetime', 'gen']].as_matrix().tolist()

        start_time = time.time()

        write2dListToCSV(windCfsDtHr, os.path.join(resultsDir, 'windGenDtUC' + zone + str(ucYear) + '.csv'))
        write2dListToCSV(windCfsDtSubhr,
                         os.path.join(resultsDir, 'windGenDtSubhrUC' + zone + str(ucYear) + '.csv'))

        write2dListToCSV(solarCfsDtHr, os.path.join(resultsDir, 'solarGenDtUC' + zone + str(ucYear) + '.csv'))
        write2dListToCSV(solarCfsDtSubhr,
                         os.path.join(resultsDir, 'solarGenSubhrDtUC' + zone + str(ucYear) + '.csv'))

        # allocate wind/solar generation
        zonalHourlyWindGen[zone], zonalHourlySolarGen[zone] = windCFs, solarCFs

        zonalHourlyDfWind[zone], zonalHourlyDfSolar[zone] = windCfsDtHr, solarCfsDtHr
        zonalSubhourlyDfWind[zone], zonalSubhourlyDfSolar[zone] = windCfsDtSubhr, solarCfsDtSubhr

    #
    # run multi-core over all gcms in list of gcms
    #
    # pack arguments in list so they can be passed to pool.map
    args_list = [[gcm, fleetUC, zonalDemandProfile[gcm], ucYear, currCo2Cap, genparam, reserveparam, curtailparam,
                  zonalHourlyWindGen, zonalHourlySolarGen, zonalHourlyDfWind, zonalHourlyDfSolar,
                  zonalSubhourlyDfWind, zonalSubhourlyDfSolar] for gcm in curtailparam.listgcms]

    ncores = len(curtailparam.listgcms)

    with mp.Pool(processes=ncores) as pool:
        list_results = pool.map(runUnitCommitmentSingleGcm, args_list)

    return 0


def runUnitCommitmentSingleGcm(list_args):
    """Run unit commitment simulation for one single climate simulation for the complete simulation period

    This function was prepared to be ran using parallel simulation (multicore) so its argument is a list packing
    all arguments. To reduce simulation time, some arrays that do not depend on climate simulations (in our model)
    are executed outside of this function and passed as arguments (such as solar/wind generation)

    :param list_args: (list) a list packing all arguments for this function
                      - gcm: (string) name of gcm
                      - fleetUC: (2d list) table of generator fleet
                      - zonalDemandProfile: (2d list) hourly demand for complete simulation period for each zone
                      - ucYear: (int) year of simulation
                      - currCo2Cap: (float) value of co2 cap
                      - genparam: object of type :mod:`Generalparameters`
                      - reserveparam: object of type :mod:`Reserveparameters`
                      - curtailparam: object of type :mod:`Curtailmentparameters`
                      - zonalHourlyWindGen: (2d list) hourly wind generation for complete simulation period for each zone
                      - zonalHourlySolarGen: (2d list) hourly solar generation for complete simulation period for each zone
                      - zonalHourlyDfWind: (pd data frame) hourly wind generation for complete simulation period for each zone
                      - zonalHourlyDfSolar: (pd data frame) hourly solar generation for complete simulation period for each zone
                      - zonalSubhourlyDfWind: (pd data frame) sub-hourly wind generation for complete simulation period for each zone
                      - zonalSubhourlyDfSolar: (pd data frame) sub-hourly solar generation for complete simulation period for each zone
    :return: this function returns the integer value 0
    """
    # unpack list with arguments fot this function
    gcm, fleetUC, zonalDemandProfile, ucYear = list_args[0], list_args[1], list_args[2], list_args[3]

    currCo2Cap, genparam, reserveparam, curtailparam = list_args[4], list_args[5], list_args[6], list_args[7]
    zonalHourlyWindGen, zonalHourlySolarGen = list_args[8], list_args[9]
    zonalHourlyDfWind, zonalHourlyDfSolar = list_args[10], list_args[11]
    zonalSubhourlyDfWind, zonalSubhourlyDfSolar = list_args[12], list_args[13]

    daysForUC = list(range(genparam.ucDayInitial, (genparam.ucDayEnd + 1)))
    # daysForUC = list(range(genparam.ucDayInitial, (genparam.ucDayEnd + 1), 2))

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

    # make copy of genparam and curtailparam for this gcm run
    genparam_local = copy.deepcopy(genparam)
    curtailparam_local = copy.deepcopy(curtailparam)

    # update relevant values for local genparam and curtailparam -------------------------------------------------
    #
    # define separate folder of results for this gcm run
    resultsDir = os.path.join(genparam.resultsDir, 'UC', gcm)
    if not os.path.exists(resultsDir):
        os.makedirs(resultsDir)

    # copy folder with GAMS code to resulsdir for this gcm run
    if not os.path.exists(os.path.join(resultsDir, 'GAMS')):
        shutil.copytree(os.path.join(genparam.resultsDir, 'GAMS'),
                        os.path.join(resultsDir, 'GAMS'), symlinks=False, ignore=None)

    # copy GAMS files to GAMS dir in results folder
    try:
        shutil.copy(os.path.join(genparam.dataRoot, 'GAMS', genparam.ucFilename),
                    os.path.join(resultsDir, 'GAMS', genparam.ucFilename))
    except shutil.Error as e:
        print('Error: {}'.format(e))
        # eg. source or destination doesn't exist
    except IOError as e:
        print('Error: {}'.format(e.strerror))

    # redefine resultsDir in local genparam so UC can find GAMS code
    genparam_local.resultsDir = resultsDir

    # define list of GCMs to only respective gcm
    curtailparam_local.listgcms = [gcm]
    #
    # ------ END of updating relevant values for local genparam and curtailparam ---------------------------------

    #
    # Simulate thermal deratings/curtailments (do this before aggregating zones!!) -------------------------------
    combineWindAndSolarToSinglePlant(fleetUC, genparam_local.ipmZones, genparam_local.dataRoot)

    updateFuelPricesExistingGens(fleetUC, ucYear, fuelPricesTimeSeries)
    hourlyCapacsCurtailedGens = importHourlyThermalCurtailments(fleetUC, ucYear, 'UC', resultsDir, genparam_local,
                                                                curtailparam_local)

    # get simulated capacities for all generators and remove redundant key for GCM
    hourlyCapacsAllGens = calculateHourlyCapacsWithCurtailments(fleetUC, hourlyCapacsCurtailedGens, ucYear)[gcm]
    #
    # ------ END of updating relevant values for local genparam and curtailparam ----------------------------------

    # Aggregate to single ipm zone --------------------------------------------------------------------------------
    #
    # aggregate dict of 1d lists
    def f1(d):
        y = dict()
        y['SERC'] = np.matrix([d[z][:8760] for z in d]).sum(axis=0).tolist()[0]
        return y

    # aggregate dict of 2d lists (only 2nd "column")
    def f2(d):
        y = dict()
        df = None
        for z in d:
            if df is None:
                df = pd.DataFrame(d[z][1:], columns=d[z][0])
            else:
                df2 = pd.DataFrame(d[z][1:], columns=d[z][0])
                df = pd.concat([df, df2]).groupby(['datetime'])['totalGen(MWh)'].sum().reset_index()
        y['SERC'] = [list(df.columns)] + df.as_matrix().tolist()

        return y

    zonalDemandProfile = f1(zonalDemandProfile)
    zonalHourlyWindGen = f1(zonalHourlyWindGen)
    zonalHourlySolarGen = f1(zonalHourlySolarGen)
    zonalHourlyDfWind = f2(zonalHourlyDfWind)
    zonalSubhourlyDfWind = f2(zonalSubhourlyDfWind)
    zonalHourlyDfSolar = f2(zonalHourlyDfSolar)
    zonalSubhourlyDfSolar = f2(zonalSubhourlyDfSolar)

    # change list of ipmzones to only SERC
    genparam_local.ipmZones = ['SERC']
    genparam_local.ipmZoneNums = [1]

    # change all zones to SERC in the fleet matrix
    colPlantType = fleetUC[0].index('Region Name')
    for i in range(1, len(fleetUC)):
        fleetUC[i][colPlantType] = 'SERC'

    # remove transmission lines
    # (create dummy line with source == sink == 'SERC' and capacity == 0)
    genparam_local.lineList = ['SERC_to_SERC']
    genparam_local.lineCapacs = {'SERC_to_SERC': 0}

    #
    # ------ END aggregating to single ipm zone-------------------------------------------------------------------

    # write hourly demand for each zone for this gcm
    writeDictToCSV(zonalDemandProfile, os.path.join(resultsDir, 'demandFullYrZonalCE' + str(ucYear) + '.csv'))

    (contResHourly, regUpHourly, regDownHourly, flexResHourly, allRes, regUpWind, regDownWind, regUpSolar,
     regDownSolar, flexWind, flexSolar) = (dict(), dict(), dict(), dict(), dict(), dict(), dict(), dict(), dict(),
                                           dict(), dict())

    for zone in genparam_local.ipmZones:
        (zonalcontResHourly, zonalregUpHourly, zonalregDownHourly, zonalflexResHourly, zonalallRes, zonalregUpWind,
         zonalregDownWind, zonalregUpSolar, zonalregDownSolar, zonalflexWind,
         zonalflexSolar) = calcWWSISReserves(zonalHourlyDfWind[zone], zonalSubhourlyDfWind[zone],
                                             zonalHourlyDfSolar[zone], zonalSubhourlyDfSolar[zone],
                                             zonalDemandProfile[zone], regLoadFrac, contLoadFrac,
                                             regErrorPercentile, flexErrorPercentile)

        contResHourly[zone] = zonalcontResHourly
        regUpHourly[zone] = zonalregUpHourly
        regDownHourly[zone] = zonalregDownHourly
        flexResHourly[zone] = zonalflexResHourly
        allRes[zone] = zonalallRes
        regUpWind[zone] = zonalregUpWind
        regDownWind[zone] = zonalregDownWind
        regUpSolar[zone] = zonalregUpSolar
        regDownSolar[zone] = zonalregDownSolar
        flexWind[zone] = zonalflexWind
        flexSolar[zone] = zonalflexSolar

    # convert allRes dict to data frame for printing
    writeDictToCSV(allRes, os.path.join(resultsDir, 'reservesUC' + str(ucYear) + '.csv'))

    writeDictToCSV(contResHourly, os.path.join(resultsDir, 'reservesContUC' + str(ucYear) + '.csv'))
    writeDictToCSV(regUpHourly, os.path.join(resultsDir, 'reservesRegUpUC' + str(ucYear) + '.csv'))
    writeDictToCSV(regDownHourly, os.path.join(resultsDir, 'reservesRegDownUC' + str(ucYear) + '.csv'))
    writeDictToCSV(flexResHourly, os.path.join(resultsDir, 'reservesFlexUC' + str(ucYear) + '.csv'))

    if genparam_local.calculateCO2Price:
        co2Price = convertCo2CapToPrice(fleetUC, zonalHourlyWindGen, zonalHourlySolarGen, zonalDemandProfile,
                                        currCo2Cap, scaleMWtoGW, scaleDollarsToThousands, scaleLbToShortTon,
                                        genparam_local.dataRoot, ucYear)
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

    sysResults = initializeSystemResultsdf()

    msAndSs = [['day', 'ms', 'ss']]  # store modelstat & solvestat from GAMS

    hydroPotentials = importHydroPotentialGen(ucYear, genparam_local, gcm)

    for dayIdx in range(0, len(daysForUC), daysOpt):

        day = daysForUC[dayIdx]

        print('---------------------------------------------------------------')
        print('Starting UC simulation of day {0:4d}'.format(day))

        daysForUCAux = list(range(day, day + daysOpt + daysLA))

        (demandUC, hourlyWindGenUC, hourlySolarGenUC, hoursForUC) = getDemandAndREGenForUC(day, daysOpt, daysLA,
                                                                                           zonalDemandProfile,
                                                                                           zonalHourlyWindGen,
                                                                                           zonalHourlySolarGen)

        (regUpUC, regDownUC, flexUC, contUC) = getResForUC(day, daysOpt, daysLA, regUpHourly, regDownHourly,
                                                           flexResHourly, contResHourly)

        hourlyCapacsUC = getHourlyCapacitiesForDays(fleetUC, hourlyCapacsAllGens, hoursForUC)

        hydroPotentialUC = getDailyHydroPotentialsUC(fleetUC, hydroPotentials, daysForUCAux, ucYear)

        if daysForUC[0] == day:  # first day, therefore no initial conditions defined. MW energy values
            (onOffInitial, genAboveMinInitial, mdtCarriedInitial, initSocDict) = setInitCondsFirstUC(fleetUC, genparam_local)
        else:  # other days. MW energy values
            (onOffInitial, genAboveMinInitial, mdtCarriedInitial, initSocDict) = setInitCondsPerPriorUC(ucModel, fleetUC,
                                                                                                        hoursForUC, daysOpt, daysLA,
                                                                                                        scaleMWtoGW)
        t0 = time.time()
        ucModel, ms, ss = callUnitCommitment(fleetUC, hourlyCapacsUC, hourlyWindGenUC, hourlySolarGenUC,
                                             hydroPotentialUC, demandUC, hoursForUC, onOffInitial, genAboveMinInitial,
                                             mdtCarriedInitial, initSocDict, co2Price, regUpUC, regDownUC, flexUC,
                                             contUC, genparam_local, reserveparam)

        print('Time (secs) for UC day ' + str(day) + ': ' + str(time.time() - t0))
        print('---------------------------------------------------------------')
        print()

        msAndSs.append([day, ms, ss])

        if int(ms) not in [1, 8] or int(ss) != 1:
            print('******* ERROR IN OPTIMIZATION! Day {0}. GCM {1} *******'.format(day, gcm))
            # if GAMS/CPLEX returns invalid results, replace results in ucModel with the one from previous day
            ucModel = copy.copy(lastUcModel)
        else:
            lastUcModel = copy.copy(ucModel)

        # just saves GAMS model of current day to list
        # ucResultsByDay.append((day, ucModel))

        saveHourlyResultsByPlant(genByPlant, regUpByPlant, regDownByPlant, flexByPlant, contByPlant, turnonByPlant,
                                 turnoffByPlant, onOffByPlant, genToRow, hourToColPlant, ucModel, day, daysOpt)

        saveHourlyPumpedHydroResults(pumphydroSoc, pumphydroCharge, ucModel, day, daysOpt)

        sysResults = saveHourlySystemResults(sysResults, ucModel, day, daysOpt)

        # check if file with flag to save GDX files exists in GAMS folder
        gamsFileDir = os.path.join(genparam_local.resultsDir, 'GAMS')
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

            sysResults.to_csv(os.path.join(resultsDir, 'systemResultsUC' + str(ucYear) + '.csv'))

            write2dListToCSV(msAndSs, os.path.join(resultsDir, 'msAndSsUC' + str(ucYear) + '.csv'))

            # save results of current UC model object to GDX file to allow "Cold start" in case of error
            write2dListToCSV([[day]], os.path.join(resultsDir, 'lastdayUC' + str(ucYear) + '.csv'))
            ucModel.get_out_db().export(os.path.join(resultsDir, "last_ucmodel.gdx"))

    # write final total results
    writeHourlyResultsByPlant(genByPlant, regUpByPlant, regDownByPlant, flexByPlant, contByPlant, turnonByPlant,
                              turnoffByPlant, onOffByPlant, resultsDir, ucYear, 'UC', 'Plant')

    writeHourlyStoResults(pumphydroCharge, pumphydroSoc, resultsDir, ucYear)

    sysResults.to_csv(os.path.join(resultsDir, 'systemResultsUC' + str(ucYear) + '.csv'))

    write2dListToCSV(msAndSs, os.path.join(resultsDir, 'msAndSsUC' + str(ucYear) + '.csv'))

    return 0


def callUnitCommitment(fleetUC, hourlyCapacsUC, hourlyWindGenUC, hourlySolarGenUC, hydroPotentialUC, demandUC,
                       hoursForUC, onOffInitial, genAboveMinInitial, mdtCarriedInitial, initSoc, co2Price, regUpUC,
                       regDownUC, flexUC, contUC, genparam, reserveparam):
    """Run the unit commintment and economic dispatch optimization for a single day

    This function is called sequentially in order to simulate dispatch over all days in the complete simulation
    horizon

    :param fleetUC: 2d list with generator fleet for UCED simulation
    :param hourlyCapacsUC:
    :param hourlyWindGenUC:
    :param hourlySolarGenUC:
    :param hydroPotentialUC:
    :param demandUC:
    :param hoursForUC: 1d list with hours (1-8760) that will be used in this UCED run
    :param onOffInitial:
    :param genAboveMinInitial:
    :param mdtCarriedInitial:
    :param initSoc:
    :param co2Price:
    :param regUpUC:
    :param regDownUC:
    :param flexUC:
    :param contUC:
    :param genparam: object of class :mod:`Generalparameters`
    :param reserveparam: object of type :mod:`Reserveparameters`
    :return:

    """
    gamsFileDir = os.path.join(genparam.resultsDir, 'GAMS')
    gamsSysDir = genparam.pathSysGams.strip()

    if not gamsSysDir == '':
        wsUC = GamsWorkspace(working_directory=gamsFileDir, system_directory=gamsSysDir)
    else:
        wsUC = GamsWorkspace(working_directory=gamsFileDir)

    dbUC = wsUC.add_database()

    # Add sets and parameters to GAMS database

    x = addSetsToDatabaseUC(dbUC, fleetUC, hoursForUC, genparam)

    addParametersToDatabaseUC(dbUC, hourlyCapacsUC, hourlyWindGenUC, hourlySolarGenUC, hydroPotentialUC, demandUC,
                              fleetUC, onOffInitial, genAboveMinInitial, mdtCarriedInitial, initSoc, co2Price,
                              hoursForUC, flexUC, contUC, regUpUC, regDownUC, genparam, reserveparam)

    # Load and run GAMS model
    ucModel = wsUC.add_job_from_file(genparam.ucFilename)
    optUC = GamsOptions(wsUC)
    optUC.defines['gdxincname'] = dbUC.name
    ucModel.run(optUC, databases=dbUC, output=sys.stdout)
    ms, ss = ucModel.out_db['pModelstat'].find_record().value, ucModel.out_db['pSolvestat'].find_record().value

    if int(ms) not in [1, 8] or int(ss) != 1:
        print('************************************************************************************')
        print('Modelstat & solvestat:', ms, ' & ', ss, ' (should be 8 and 1)')
        print('************************************************************************************')

    return ucModel, ms, ss


def addSetsToDatabaseUC(db, fleetUC, hoursForUC, genparam):
    """Add sets to GAMS database

    This function adds the sets to the GAMS database needed to run the UCED optimization

    :param db: gams database object
    :param fleetUC: 2d list with generator fleet for UCED simulation
    :param hoursForUC: 1d list with hours (1-8760) that will be used in this UCED run
    :param genparam: object of class :mod:`Generalparameters`
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
                              fleetUC, onOffInitial, genAboveMinInitial, mdtCarriedInitial, initSoc, co2Price,
                              hoursForUC, flexUC, contUC, regUpUC, regDownUC, genparam, reserveparam):
    """Add parameters to GAMS database

    This function adds the parameters to the GAMS database needed to run the UCED optimization

    :param db: gams database object
    :param hourlyCapacsUC:
    :param hourlyWindGenUC:
    :param hourlySolarGenUC:
    :param hydroPotentialUC:
    :param demandUC:
    :param fleetUC: 2d list with generator fleet for UCED simulation
    :param onOffInitial:
    :param genAboveMinInitial:
    :param mdtCarriedInitial:
    :param initSoc: initial state of charge for pumped storage
    :param co2Price: (float) price of co2
    :param hoursForUC: 1d list with hours (1-8760) that will be used in this UCED run
    :param flexUC:
    :param contUC:
    :param regUpUC:
    :param regDownUC:
    :param genparam: object of class :mod:`Generalparameters`
    :param reserveparam: object of type :mod:`Reserveparameters`
    """

    # get needed sets and set symbols from data base
    genSet = db.get_set('egu')
    genSymbols = [symb.keys[0] for symb in genSet.__iter__()]

    hydroGenSet = db.get_set('hydroegu')

    hourSet = db.get_set('h')

    zoneSet = db.get_set('z')

    lineSet = db.get_set('l')

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

    addCostNonservedEnergy(db, genparam.cnse, genparam.scaleMWtoGW, genparam.scaleDollarsToThousands)

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

    :param genparam: object of class :mod:`Generalparameters`
    :param reserveparam: object of type :mod:`Reserveparameters`
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

    outstr = outstr + 'CO2 Cap case         : {:.0%}\n'.format(genparam.co2CapPercentage)

    outstr = outstr + '\n\n'

    outstr = outstr + 'Transmission lines capacities:\n'

    for k, value in genparam.lineCapacs.items():
        outstr = outstr + '{0:<21}: {1: .2f} MW\n'.format(k, value)

    outstr = outstr + '\n\n'

    fname = os.path.join(os.path.expanduser(genparam.resultsDir), 'DESCRIPTION.TXT')

    with open(fname, 'w') as f:
        f.write(outstr)
