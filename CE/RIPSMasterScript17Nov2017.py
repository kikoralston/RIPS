# Michael Craig, 31 May 2016
# Run capacity expansion and dispatch models w/ climate impacts.

import os, csv, operator, copy, time, random
import numpy as np
import datetime as dt
from gams import *
from AuxFuncs import *
from GAMSAuxFuncs import *
from SetupGeneratorFleet import setupGeneratorFleet, aggregatePlantTypeToORIS
# from ImportDemandProfile import importZonalHourlyDemand
from ForecastDemandWithRegression import forecastZonalDemandWithReg
from TransmissionLineFuncs import *
from UpdateFuelPriceFuncs import *
from DemandFuncs import *
from DemandFuncsCE import *
from GetHydroMaxGenPotential import getHydroEPotential
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
from SetupResultLists import setupHourlyResultsByPlant, setupHourlySystemResults
from SaveHourlyResults import saveHourlyResultsByPlant, saveHourlySystemResults
from WriteUCResults import writeHourlyResultsByPlant, writeHourlyStoResults
from ReservesWWSIS import calcWWSISReserves
from SaveCEOperationalResults import saveCapacExpOperationalData
from LoadCEFleet import loadCEFleet
from ModifyGeneratorCapacityWithWaterTData import (determineHrlyCurtailmentsForExistingGens,
                                                   getGenToCellAndCellToGenDictionaries,
                                                   createBaseFilenameToReadOrWrite, getCellLatAndLongFromFolderName,
                                                   processRBMDataIntoIndividualCellFiles)
from ModifyNewTechCapacityWithWaterTData import determineHrlyCurtailmentsForNewTechs


################################################################################
###### UNIVERSAL PARAMETERS ####################################################
################################################################################
def setKeyParameters():
    # KEY PARAMETERS

    # path to folder where data is located. Script will assume a specific folder tree structure
    # to find different files.
    dataRoot = '/Volumes/KIKO64/Data'

    # whether running on XSEDE supercomputer (affects sys_dir in GAMS call)
    xsedeRun = False

    runCE, runUC = True, False
    runFirstUCYear = False
    # NOTE: if run every 5 years, leap year is in curtailment data.
    (startYear, endYear, yearStepCE) = (2015, 2030, 10)
    firstUCYear = 2015
    daysPerSeason = 10
    daysForUC = [val for val in range(1, 2)]  # start @ 1, go to 1 past end day (e.g., 366 for full year)
    analysisArea = 'TVA'  # coreSERC, allSERC, onlyTN, test, TVA
    # KEY CURTAILMENT PARAMETERS
    incCurtailments, incRegs = True, True  # whether to model curtailments,whether to model env regs on water T
    coolDesignT = 100  # design temperature of cooling techs
    # PTs curtailed via regression (ptCurtailed) and via enviro regulations (ptCurtailedRegs)
    ptCurtailed = {'Coal Steam', 'Combined Cycle'}
    ptCurtailedRegs = {'Coal Steam', 'Coal Steam CCS', 'Combined Cycle', 'Combined Cycle CCS', 'Nuclear'}
    # Determines which cells new plants can be sited in; can be all all cells or those w/ gens already
    cellsEligibleForNewPlants = 'withGens'  # 'all' (all cells) or 'withGens' (only cells already with gen inside)
    # Of those cells given above parameter, whether input all of them or just the one w/ max wtr T per zone.
    cellNewTechCriteria = 'all'  # 'all' or 'maxWaterT' (maxWaterT: only include cell w/ max water T in CE.)
    # GENERAL PARAMETERS
    compressFleet = True
    co2CapScenario = 'none'
    scenario = 'normal'
    co2CapEndYr, co2CapEnd = getCo2Cap(co2CapScenario)
    fuelPricesTimeSeries = importFuelPrices(dataRoot, scenario)
    # RESULTS DIRECTORY
    resultsDir = '/Users/kiko/Documents/CE'
    folderName = ('Area' + analysisArea + 'Cells' + cellNewTechCriteria +
                  ('Curtail' if incCurtailments == True else 'NoCurtail') +
                  ('EnvRegs' if incRegs == True else 'NoRegs') +
                  'C' + co2CapScenario + 'S' + scenario[:3])
    resultsDir = os.path.join(resultsDir, folderName)
    if not os.path.exists(resultsDir): os.makedirs(resultsDir)
    # THERMAL CURTAILMENT PARAMETERS
    processRBMData = False  # set to True if first time improting set of RBM data
    # RENEWABLE CAPACITY FACTOR PARAMETERS
    tzAnalysis = 'EST'
    projectName = 'rips'
    windGenDataYr = 2009  # picked somewhat arbitrarily - average CF of all years in WIND data
    # CAPACITY EXPANSION PARAMETERS
    # capacExpFilename = 'CERIPS24July2017.gms'
    capacExpFilename = 'CERIPS24July2017PS.gms'
    maxAddedZonalCapacPerTech = 10000
    incITC = True
    retirementCFCutoff = .3  # retire units w/ CF lower than given value
    ptEligRetCF = ['Coal Steam']
    selectCurtailDays = True
    planningReserve = 0.15  # fraction of peak demand
    discountRate = 0.07  # fraction
    allowCoalWithoutCCS, onlyNSPSUnits, permitOncethru = True, False, True
    # UNIT COMMITMENT PARAMETERS
    ucFilename = 'UCRIPS18April2017.gms'
    calculateCO2Price = True
    (daysOpt, daysLA) = (1, 1)
    ocAdderMin, ocAdderMax = 0, 0.05  # $/MWh
    # PUMPED HYDRO PARAMETERS
    phEff, phMaxSoc, phInitSoc = .81, 5, .5  # max soc as multiple of capacity; init SOC as fraction of max soc
    # CONVERSION PARAMETERS
    scaleMWtoGW = 1000
    scaleDollarsToThousands = 1000
    scaleLbToShortTon = 2000
    # TEMPORARY PARAMETERS
    annualDemandGrowth = 0.02  # fraction per year
    # AREA OF ANALYSIS
    if analysisArea == 'coreSERC':
        states = ['North Carolina', 'South Carolina', 'Georgia', 'Mississippi', 'Alabama',
                  'Kentucky', 'Tennessee']
        statesAbbrev = ['NC', 'SC', 'GA', 'MS', 'AL', 'KY', 'TN']
        ipmZones = ['S_SOU', 'S_VACA', 'S_C_KY', 'S_C_TVA']
    elif analysisArea == 'allSERC':
        states = ['North Carolina', 'South Carolina', 'Georgia', 'Mississippi', 'Alabama',
                  'Kentucky', 'Tennessee', 'Missouri', 'Virginia', 'Illinois', 'Louisiana']
        statesAbbrev = ['NC', 'SC', 'GA', 'MS', 'AL', 'KY', 'TN', 'MO', 'VA', 'IL', 'LA']
        ipmZones = ['S_SOU', 'S_VACA', 'S_C_KY', 'S_C_TVA', 'S_D_WOTA',
                    'S_D_N_AR', 'S_D_AMSO', 'S_D_REST']
    elif analysisArea == 'onlyTN':
        (states, statesAbbrev) = (['Tennessee'], ['TN'])
        ipmZones = ['S_C_TVA']
    elif analysisArea == 'test':
        states = ['North Carolina', 'South Carolina', 'Georgia', 'Mississippi', 'Alabama']
        statesAbbrev = ['NC', 'SC', 'GA', 'MS', 'AL']
        ipmZones = ['S_SOU', 'S_VACA']
    elif analysisArea == 'TVA':
        states = ['North Carolina', 'Georgia', 'Mississippi', 'Alabama', 'Kentucky', 'Tennessee']
        statesAbbrev = ['NC', 'GA', 'MS', 'AL', 'KY', 'TN']
        ipmZones = ['S_C_TVA']
    fipsToZones, fipsToPolys = getIPMPolys(dataRoot, ipmZones)  # get IPM polygons
    statePolys = getStatePolys(dataRoot, states)  # get state polygons
    # CREATE LIST W/ IDXS THAT PARALLEL ZONES - NEED TO PRESERVE ORDER SO DON'T USE A DICT
    ipmZoneNums = [idx for idx in range(1, len(ipmZones) + 1)]
    write2dListToCSV([['Zone', 'ZoneNum']] + rotate([ipmZones, ipmZoneNums]),
                     os.path.join(resultsDir, 'zoneNamesToNumbers.csv'))
    # IMPORT TRANSMISSION SYSTEM DATA
    lineList, lineCapacs = setLines(dataRoot)
    # OLD VARIABLES
    testModel = False  # use dummy test system; not currently working
    return (dataRoot, xsedeRun, runCE, runUC, resultsDir, startYear, endYear, yearStepCE, daysPerSeason,
            testModel, cellsEligibleForNewPlants, compressFleet, co2CapEndYr, co2CapEnd,
            fuelPricesTimeSeries, processRBMData, tzAnalysis, projectName, windGenDataYr, capacExpFilename,
            retirementCFCutoff, ptEligRetCF, daysPerSeason, selectCurtailDays, planningReserve,
            discountRate, allowCoalWithoutCCS, onlyNSPSUnits, ucFilename, calculateCO2Price,
            daysOpt, daysLA, scaleMWtoGW, scaleDollarsToThousands, scaleLbToShortTon, annualDemandGrowth,
            states, statesAbbrev, ipmZones, ocAdderMin, ocAdderMax, maxAddedZonalCapacPerTech,
            lineList, lineCapacs, ptCurtailed, ptCurtailedRegs, cellNewTechCriteria, firstUCYear,
            scenario, incITC, ipmZoneNums, fipsToZones, fipsToPolys, statePolys, permitOncethru, incCurtailments,
            incRegs, phEff, phMaxSoc, phInitSoc, runFirstUCYear, coolDesignT)


def setThermalCurtailmentParameters():
    # Env reg parameters: dict of state : max water T (deg C)
    envRegMaxT = {'Alabama': 32, 'Georgia': 32, 'Mississippi': 32, 'South Carolina': 32.2,
                  'North Carolina': 32, 'Tennessee': 30.5, 'Kentucky': 31.7}
    # RBM processing parameters
    outputHeaders = ['Year', 'Month', 'Day', 'Streamflow(cfs)', 'StreamT(degC)', 'HeadwaterT(degC)', 'AirT(degC)']
    locPrecision = 4
    numCellsToProcess = 4
    tempAndSpatFilename = 'bcc-csm1-1-m_rcp45_r1i1p1'
    nsegFilename = 'Tennessee.nseg_nday'
    # Water T directories
    rbmRootDir = '/Volumes/KIKO64/Data/DatafromUW'
    rbmDataDir = os.path.join(rbmRootDir, 'RBMRawWaterTData10Aug2016')
    rbmOutputDir = os.path.join(rbmRootDir, 'RBMProcessedWaterT25Aug2016', tempAndSpatFilename)
    return (outputHeaders, locPrecision, tempAndSpatFilename, nsegFilename,
            rbmRootDir, rbmDataDir, rbmOutputDir, numCellsToProcess, envRegMaxT)


def defineReserveParameters():
    # Requirement parameters - based on WWSIS Phase 2
    regLoadFrac = .01  # frac of hourly load in reg up & down
    contLoadFrac = .03  # frac of hourly load in contingency
    regErrorPercentile = 95  # percentile of 10-m wind & 5-m solar forecast errors used in reg reserves
    flexErrorPercentile = 70  # percentile of hourly wind & solar forecast errors used in reg reserves
    # Cost coeff - from Denholm et al. 2013, val of E sto in grid apps
    regUpCostCoeffs = {'Combined Cycle': 6, 'Combined Cycle CCS': 6, 'O/G Steam': 4,
                       'Coal Steam': 10, 'Coal Steam CCS': 10}  # $/MWh
    # Timeframes
    regReserveMinutes = 5  # reg res must be provided w/in 5 mins
    flexReserveMinutes = 10  # spin reserves must be provided w/in 10 minutes
    contingencyReserveMinutes = 30  # contingency res must be provided w/in 30 minutes
    minutesPerHour = 60
    rampRateToRegReserveScalar = regReserveMinutes / minutesPerHour  # ramp rate in MW/hr
    rampRateToFlexReserveScalar = flexReserveMinutes / minutesPerHour  # ramp rate in MW/hr
    rampRateToContReserveScalar = contingencyReserveMinutes / minutesPerHour
    return (regLoadFrac, contLoadFrac, regErrorPercentile, flexErrorPercentile,
            regUpCostCoeffs, rampRateToRegReserveScalar, rampRateToFlexReserveScalar,
            rampRateToContReserveScalar)


def importFuelPrices(dataRoot, scenario):
    fuelPriceDir = os.path.join(dataRoot, 'FuelPricesCapacityExpansion')

    if scenario == 'ng':
        fuelFileName = 'FuelPriceTimeSeries2Aug2016LowNG.csv'
    else:
        fuelFileName = 'FuelPriceTimeSeries2Aug2016.csv'

    return readCSVto2dList(os.path.join(fuelPriceDir, fuelFileName))


################################################################################
###### MASTER FUNCTION #########################################################
################################################################################
def masterFunction():
    # Load parameters
    (dataRoot, xsedeRun, runCE, runUC, resultsDir, startYear, endYear, yearStepCE, daysPerSeason,
     testModel, cellsEligibleForNewPlants, compressFleet, co2CapEndYr, co2CapEnd,
     fuelPricesTimeSeries, processRBMData, tzAnalysis, projectName, windGenDataYr, capacExpFilename,
     retirementCFCutoff, ptEligRetCF, daysPerSeason, selectCurtailDays, planningReserve,
     discountRate, allowCoalWithoutCCS, onlyNSPSUnits, ucFilename, calculateCO2Price,
     daysOpt, daysLA, scaleMWtoGW, scaleDollarsToThousands, scaleLbToShortTon, annualDemandGrowth,
     states, statesAbbrev, ipmZones, ocAdderMin, ocAdderMax, maxAddedZonalCapacPerTech,
     lineList, lineCapacs, ptCurtailed, ptCurtailedRegs, cellNewTechCriteria, firstUCYear,
     scenario, incITC, ipmZoneNums, fipsToZones, fipsToPolys, statePolys, permitOncethru,
     incCurtailments, incRegs, phEff, phMaxSoc, phInitSoc, runFirstUCYear, coolDesignT) = setKeyParameters()

    (regLoadFrac, contLoadFrac, regErrorPercentile, flexErrorPercentile, regUpCostCoeffs,
     rrToRegTime, rrToFlexTime, rrToContTime) = defineReserveParameters()

    genFleet = getInitialFleetAndDemand(testModel, states, ipmZones, endYear, startYear,
                                        fuelPricesTimeSeries, compressFleet, dataRoot, resultsDir, ocAdderMin,
                                        ocAdderMax, regUpCostCoeffs,
                                        ptCurtailed | ptCurtailedRegs)  # combine curtailed sets
    print('Set up initial data')

    if processRBMData: processRBMDataToUsableFormat()

    (capacExpModelsEachYear, capacExpBuilds, capacExpGenByGens, capacExpRetiredUnitsByCE,
     capacExpRetiredUnitsByAge) = ([], [['TechnologyType']], [['ORIS+UnitID']], [], [])

    # begin loop
    for currYear in range(startYear, endYear, yearStepCE):
        currCo2Cap = interpolateCO2Cap(currYear, co2CapEndYr, co2CapEnd) * 1E3
        print('Curr CO2 cap:', currCo2Cap)
        zonalDemandProfile, zonalTempDfs = forecastZonalDemandWithReg(currYear, dataRoot, ipmZones, resultsDir)

        if currYear > firstUCYear and runCE == True:
            if currYear == startYear + yearStepCE:
                priorCapacExpModel, priorHoursCE, genFleetPriorCE = None, None, None  # first CE run

            (genFleet, genFleetNoRetiredUnits, genFleetPriorCE, priorCapacExpModel,
             priorHoursCE) = runCapacityExpansion(genFleet, zonalDemandProfile, startYear, currYear, endYear,
                                                  planningReserve, discountRate, fuelPricesTimeSeries, states,
                                                  statesAbbrev, scaleMWtoGW, scaleDollarsToThousands, currCo2Cap,
                                                  allowCoalWithoutCCS, capacExpFilename, onlyNSPSUnits, daysPerSeason,
                                                  retirementCFCutoff, scaleLbToShortTon, tzAnalysis, projectName,
                                                  dataRoot, resultsDir, capacExpModelsEachYear, capacExpBuilds,
                                                  capacExpGenByGens, capacExpRetiredUnitsByCE,
                                                  capacExpRetiredUnitsByAge, ocAdderMin, ocAdderMax,
                                                  maxAddedZonalCapacPerTech, windGenDataYr, regUpCostCoeffs, xsedeRun,
                                                  scenario, ptEligRetCF, genFleetPriorCE, priorCapacExpModel,
                                                  priorHoursCE, incITC, selectCurtailDays, cellsEligibleForNewPlants,
                                                  ipmZones, lineList, lineCapacs, cellNewTechCriteria, ptCurtailed,
                                                  ptCurtailedRegs, ipmZoneNums, fipsToZones, fipsToPolys, statePolys,
                                                  permitOncethru, incCurtailments, incRegs, phEff, phMaxSoc,
                                                  phInitSoc, coolDesignT)
        if runUC:
            # Either only runs 2015, or runs in all but 2015
            if ((currYear == firstUCYear and runFirstUCYear is True) or
                    (currYear > firstUCYear and runFirstUCYear is False)):
                if currYear == firstUCYear:
                    genFleetNoRetiredUnits = genFleet
                else:
                    genFleetNoRetiredUnits = loadCEFleet(currYear, resultsDir)
                (ucResultsByDay, hourlyGenerationByPlants) = runUnitCommitment(genFleetNoRetiredUnits,
                                                                               demandWithGrowth, startYear, currYear,
                                                                               fuelPricesTimeSeries, states,
                                                                               statesAbbrev, scaleMWtoGW,
                                                                               scaleDollarsToThousands, currCo2Cap,
                                                                               calculateCO2Price, scaleLbToShortTon,
                                                                               daysForUC, daysOpt, daysLA, tzAnalysis,
                                                                               projectName, dataRoot, resultsDir,
                                                                               ocAdderMin, ocAdderMax, windGenDataYr,
                                                                               regLoadFrac, contLoadFrac,
                                                                               regErrorPercentile, flexErrorPercentile,
                                                                               rrToRegTime, rrToFlexTime, rrToContTime,
                                                                               copy.deepcopy(regUpCostCoeffs), xsedeRun,
                                                                               runCE, scenario, ucFilename, ipmZones,
                                                                               incCurtailments)


################################################################################
####### SET UP INITIAL FLEET AND DEMAND ########################################
################################################################################
def getInitialFleetAndDemand(testModel, states, ipmZones, endYear,
                             startYear, fuelPricesTimeSeries, compressFleet, dataRoot, resultsDir, ocAdderMin,
                             ocAdderMax, regUpCostCoeffs, ptCurtailedAll):
    genFleet = setupGeneratorFleet(testModel, states, ipmZones, endYear, startYear,
                                   fuelPricesTimeSeries, compressFleet, dataRoot, ocAdderMin, ocAdderMax,
                                   regUpCostCoeffs,
                                   ptCurtailedAll)
    aggregatePlantTypeToORIS(genFleet, 'Hydro')
    aggregatePlantTypeToORIS(genFleet, 'Pumped Storage')
    # Add placeholder columns for additions & retirements by CE
    ceHeaders = ['YearAddedCE', 'YearRetiredByCE', 'YearRetiredByAge']
    genFleet[0].extend(ceHeaders)
    for row in genFleet[1:]: row.extend([''] * len(ceHeaders))
    write2dListToCSV(genFleet, os.path.join(resultsDir, 'genFleetInitial.csv'))
    return genFleet


################################################################################
####### PROCESS RBM WATER AND TEMPERATURE DATA #################################
################################################################################
def processRBMDataToUsableFormat():
    print('Processing raw RBM data')
    (outputHeaders, locPrecision, tempAndSpatFilename, nsegFilename,
     rbmRootDir, rbmDataDir, rbmOutputDir, numCellsToProcess, envRegMaxT) = setThermalCurtailmentParameters()
    processRBMDataIntoIndividualCellFiles(rbmDataDir, tempAndSpatFilename,
                                          rbmOutputDir, nsegFilename, locPrecision, outputHeaders, numCellsToProcess)
    print('Processed RBM data')


################################################################################
####### RUN CAPACITY EXPANSION #################################################
################################################################################
def runCapacityExpansion(genFleet, zonalDemandProfile, startYear, currYear, endYear,
                         planningReserve, discountRate, fuelPricesTimeSeries, states, statesAbbrev,
                         scaleMWtoGW, scaleDollarsToThousands, currCo2Cap, allowCoalWithoutCCS,
                         capacExpFilename, onlyNSPSUnits, daysPerSeason, retirementCFCutoff, scaleLbToShortTon,
                         tzAnalysis, projectName, dataRoot, resultsDirOrig, capacExpModelsEachYear, capacExpBuilds,
                         capacExpGenByGens, capacExpRetiredUnitsByCE, capacExpRetiredUnitsByAge,
                         ocAdderMin, ocAdderMax, maxAddedZonalCapacPerTech, windGenDataYr, regUpCostCoeffs,
                         xsedeRun, scenario, ptEligRetCF, genFleetPriorCE, priorCapacExpModel,
                         priorHoursCE, incITC, selectCurtailDays, cellsEligibleForNewPlants, ipmZones,
                         lineList, lineCapacs, cellNewTechCriteria, ptCurtailed, ptCurtailedRegs, ipmZoneNums,
                         fipsToZones, fipsToPolys, statePolys, permitOncethru, incCurtailments, incRegs,
                         phEff, phMaxSoc, phInitSoc, coolDesignT):
    # get curtailment parameters
    (outputHeaders, locPrecision, tempAndSpatFilename, nsegFilename, rbmRootDir,
     rbmDataDir, rbmOutputDir, numCellsToProcess, envRegMaxT) = setThermalCurtailmentParameters()

    resultsDir = os.path.join(resultsDirOrig, 'CE')

    if not os.path.exists(resultsDir): os.makedirs(resultsDir)
    print('Entering CE loop for year {0:4d}'.format(currYear))

    write2dListToCSV([[currCo2Cap]], os.path.join(resultsDir, 'co2CapCE' + str(currYear) + '.csv'))
    newTechsCE = getNewTechs(allowCoalWithoutCCS, onlyNSPSUnits, regUpCostCoeffs, currYear,
                             dataRoot, resultsDir, scenario, incITC, permitOncethru)
    updateFuelPrices(genFleet, newTechsCE, currYear, fuelPricesTimeSeries)
    write2dListToCSV(newTechsCE, os.path.join(resultsDir, 'newTechsCE' + str(currYear) + '.csv'))
    if priorCapacExpModel is not None:  # if not in first CE loop
        unitsRetireCFPriorCE = retireUnitsCFPriorCE(genFleet, genFleetPriorCE, retirementCFCutoff,
                                                    priorCapacExpModel, priorHoursCE, scaleMWtoGW, ptEligRetCF,
                                                    currYear)
        print('Units that retire due to econ from prior CE ' + str(currYear) + ':', unitsRetireCFPriorCE)
        write2dListToCSV([unitsRetireCFPriorCE],
                         os.path.join(resultsDir, 'genRetirementsEconCEPrior' + str(currYear) + '.csv'))
    genFleetForCE = createFleetForCurrentCELoop(genFleet, currYear, capacExpRetiredUnitsByAge,
                                                dataRoot, scenario)  # removes all retired units
    print('Num units that retire due to age in ' + str(currYear) + ':' + str(len(capacExpRetiredUnitsByAge[-1]) - 1))
    write2dListToCSV(genFleetForCE, os.path.join(resultsDir, 'genFleetForCEPreRECombine' + str(currYear) + '.csv'))
    combineWindAndSolarToSinglePlant(genFleetForCE, ipmZones, dataRoot)  # combines by zone
    write2dListToCSV(genFleetForCE, os.path.join(resultsDir, 'genFleetForCE' + str(currYear) + '.csv'))
    (startWindCapacForCFs, startSolarCapacForCFs) = (0, 0)
    writeDictToCSV(zonalDemandProfile, os.path.join(resultsDir, 'demandFullYrZonalCE' + str(currYear) + '.csv'))
    zoneCol = genFleetForCE[0].index('Region Name')
    zonalHourlyWindGen, zonalHourlySolarGen = dict(), dict()
    zonalNewWindCFs, zonalNewSolarCFs, zonalNetDemand = dict(), dict(), dict()
    for zone in ipmZones:
        print('Zone ', zone)
        # Check 5/18/17: next 2 lines are working properly!
        zonalGenFleet = [genFleetForCE[0]] + [row for row in genFleetForCE if row[zoneCol] == zone]
        zoneDemand = zonalDemandProfile[zone]
        (windCFs, windCfsDtHr, windCfsDtSubhr, ewdIdAndCapac, solarCFs, solarCfsDtHr, solarCfsDtSubhr,
         solarFilenameAndCapac) = getRenewableCFs(zonalGenFleet, startWindCapacForCFs, startSolarCapacForCFs,
                                                  tzAnalysis, dataRoot, windGenDataYr, zone, fipsToZones, fipsToPolys)
        print('Got RE CFs')
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
         addedWindCapac, addedSolarCapac) = getNewWindAndSolarCFs(zonalGenFleet,
                                                                  currYear, 'CE', tzAnalysis, dataRoot, resultsDir,
                                                                  windGenDataYr, zone, fipsToZones, fipsToPolys)
        print('Got new RE CFs')
        write2dListToCSV(newWindIdAndCapac,
                         os.path.join(resultsDir, 'windNewIdAndCapacCE' + zone + str(currYear) + '.csv'))
        write2dListToCSV(newSolarFilenameAndCapac,
                         os.path.join(resultsDir, 'solarNewIdAndCapacCE' + zone + str(currYear) + '.csv'))
        (netDemand, hourlyWindGen, hourlySolarGen) = getNetDemand(zoneDemand, windCFs,
                                                                  ewdIdAndCapac, solarCFs, solarFilenameAndCapac,
                                                                  currYear, 'CE', resultsDir)
        zonalHourlyWindGen[zone], zonalHourlySolarGen[zone] = hourlyWindGen, hourlySolarGen
        zonalNewWindCFs[zone], zonalNewSolarCFs[zone] = newWindCFs, newSolarCFs
        zonalNetDemand[zone] = netDemand

    writeDictToCSV(zonalHourlyWindGen, os.path.join(resultsDir, 'windGenFullYrCE' + str(currYear) + '.csv'))
    writeDictToCSV(zonalHourlySolarGen, os.path.join(resultsDir, 'solarGenFullYrCE' + str(currYear) + '.csv'))
    writeDictToCSV(zonalNewWindCFs, os.path.join(resultsDir, 'windNewCFsFullYrCE' + str(currYear) + '.csv'))
    writeDictToCSV(zonalNewSolarCFs, os.path.join(resultsDir, 'solarNewCFsFullYrCE' + str(currYear) + '.csv'))
    writeDictToCSV(zonalNetDemand, os.path.join(resultsDir, 'demandNetFullYrCE' + str(currYear) + '.csv'))
    hrlyCurtailmentsAllGensInTgtYr = importHourlyThermalCurtailments(genFleetForCE,  # existing gens only!
                                                                     currYear, 'CE', ptCurtailed, ptCurtailedRegs,
                                                                     resultsDir, incCurtailments, incRegs,
                                                                     envRegMaxT, dataRoot, coolDesignT)
    (demandCE, hoursForCE, repHrsBySeason, specialHrs, regHrsBySeason, demandCEZonal,
     hourlyWindGenCEZonal, hourlySolarGenCEZonal, peakDemandHourZonal,
     planningReserveZonal, repAndSpeHoursDict) = selectWeeksForExpansion(zonalDemandProfile,
                                                                         zonalNetDemand, zonalHourlyWindGen,
                                                                         zonalHourlySolarGen, daysPerSeason,
                                                                         selectCurtailDays,
                                                                         hrlyCurtailmentsAllGensInTgtYr, currYear,
                                                                         resultsDir, planningReserve)
    write2dListToCSV([hoursForCE], os.path.join(resultsDir, 'hoursCE' + str(currYear) + '.csv'))
    write2dListToCSV([demandCE], os.path.join(resultsDir, 'demandCE' + str(currYear) + '.csv'))
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
    # Import hydropower generation potential (dict of season:ORIS ID: potential gen)
    hydroPotPerSeason = getHydroEPotential(genFleetForCE, zonalDemandProfile, repAndSpeHoursDict, currYear, dataRoot)
    for season in hydroPotPerSeason:
        writeDictToCSV(hydroPotPerSeason[season],
                       os.path.join(resultsDir, 'hydroPot' + season + 'CE' + str(currYear) + '.csv'))
    # Handle curtailments
    hourlyCapacsCE, hourlyCapacsCEList = getHourlyNonRECapacsForCE(genFleetForCE, hrlyCurtailmentsAllGensInTgtYr,
                                                                   hoursForCE, currYear)
    write2dListToCSV(hourlyCapacsCEList, os.path.join(resultsDir, 'curtailedHourlyCapacsCE' + str(currYear) + '.csv'))
    eligibleCellWaterTs = loadEligibleCellWaterTs(cellsEligibleForNewPlants, genFleet, rbmOutputDir,
                                                  locPrecision, ipmZones, fipsToZones, fipsToPolys)

    # hrlyCurtailmentsAllTechsInTgtYr is a dict of (tech,cell):[hourlycapac]
    hrlyCurtailmentsAllTechsInTgtYr, hrlyTechCurtailmentsList = \
        determineHrlyCurtailmentsForNewTechs(eligibleCellWaterTs, newTechsCE, currYear, cellNewTechCriteria,
                                             ptCurtailed, ptCurtailedRegs, fipsToZones, fipsToPolys, incCurtailments,
                                             incRegs, rbmOutputDir, locPrecision, statePolys, dataRoot, coolDesignT,
                                             envRegMaxT, resultsDir)
    write2dListToCSV(hrlyTechCurtailmentsList,
                     os.path.join(resultsDir, 'curtailmentsTechHourly' + str(currYear) + '.csv'))
    (hourlyCurtailedTechCapacsCE, curtailedTechHourlyCapacsCEList,
     curtailedTechHourlyCapacsList) = getHourlyCurtailedTechCapacsForCE(newTechsCE,
                                                                        hrlyCurtailmentsAllTechsInTgtYr, hoursForCE,
                                                                        currYear, ptCurtailed | ptCurtailedRegs)
    write2dListToCSV(curtailedTechHourlyCapacsList,
                     os.path.join(resultsDir, 'curtailedTechHourlyCapacs' + str(currYear) + '.csv'))
    write2dListToCSV(curtailedTechHourlyCapacsCEList,
                     os.path.join(resultsDir, 'curtailedTechHourlyCapacsCE' + str(currYear) + '.csv'))
    print('Got hourly generation profiles with curtailments for existing and new techs')
    # Extract set of cells for techs & assign to zones
    cellsForNewTechs = set()
    for (tech, cell) in hrlyCurtailmentsAllTechsInTgtYr: cellsForNewTechs.add(cell)
    cellsToZones = assignCellsToIPMZones(cellsForNewTechs, fipsToZones, fipsToPolys)
    writeDictToCSV(cellsToZones, os.path.join(resultsDir, 'cellsToZonesCE' + str(currYear) + '.csv'))

    # Run CE
    print('Set inputs, running CE model...')
    t0 = time.time()
    (capacExpModel, ms, ss, curtailedTechs, renewTechs, notCurtailedTechs,
     pumpedHydroSymbs) = callCapacityExpansion(genFleetForCE, hourlyCapacsCE, hourlyCurtailedTechCapacsCE,
                                               hourlyWindGenCEZonal, hourlySolarGenCEZonal, demandCEZonal, newTechsCE,
                                               planningReserveZonal, discountRate, hoursForCE, newWindCFsCEZonal,
                                               newSolarCFsCEZonal, scaleMWtoGW, scaleDollarsToThousands, currCo2Cap,
                                               capacExpFilename, seasonDemandWeights, repHrsBySeason, specialHrs,
                                               peakDemandHourZonal, scaleLbToShortTon, dataRoot,
                                               maxAddedZonalCapacPerTech, xsedeRun, ipmZones, lineList, lineCapacs,
                                               cellsToZones, ipmZoneNums, ptCurtailed | ptCurtailedRegs,
                                               hydroPotPerSeason, phEff, phMaxSoc, phInitSoc)

    write2dListToCSV([['ms', 'ss'], [ms, ss]], os.path.join(resultsDir, 'msAndSsCE' + str(currYear) + '.csv'))
    print('Time (secs) for CE year ' + str(currYear) + ': ' + str(time.time() - t0))

    # Write and Save resulting data
    (genByPlant, genByCTechAndCell, genByRETechAndZone, genByNCTechAndZone, sysResults,
     co2EmAndCostResults, flowByLine, chargeByPH, socByPH) = saveCapacExpOperationalData(capacExpModel,
                                                                                         genFleetForCE, curtailedTechs,
                                                                                         renewTechs, notCurtailedTechs,
                                                                                         cellsForNewTechs, ipmZoneNums,
                                                                                         hoursForCE, lineList,
                                                                                         pumpedHydroSymbs)
    write2dListToCSV(genByPlant, os.path.join(resultsDir, 'genByPlantCE' + str(currYear) + '.csv'))
    write2dListToCSV(genByCTechAndCell, os.path.join(resultsDir, 'genByCurtTechCE' + str(currYear) + '.csv'))
    write2dListToCSV(genByRETechAndZone, os.path.join(resultsDir, 'genByRETechCE' + str(currYear) + '.csv'))
    write2dListToCSV(genByNCTechAndZone, os.path.join(resultsDir, 'genByNotCurtTechCE' + str(currYear) + '.csv'))
    write2dListToCSV(sysResults, os.path.join(resultsDir, 'systemResultsCE' + str(currYear) + '.csv'))
    write2dListToCSV(co2EmAndCostResults, os.path.join(resultsDir, 'co2EmsAndCostsCE' + str(currYear) + '.csv'))
    write2dListToCSV(flowByLine, os.path.join(resultsDir, 'lineFlowsCE' + str(currYear) + '.csv'))
    write2dListToCSV(chargeByPH, os.path.join(resultsDir, 'chargeByPumpHydroCE' + str(currYear) + '.csv'))
    write2dListToCSV(socByPH, os.path.join(resultsDir, 'socByPumpHydroCE' + str(currYear) + '.csv'))
    capacExpModelsEachYear.append((currYear, capacExpModel))
    newCurtTech, newRETech, newNotCurtTech = saveCapacExpBuilds(capacExpBuilds, capacExpModel, currYear)

    genFleet = addNewGensToFleet(genFleet, newCurtTech, newRETech, newNotCurtTech, newTechsCE,
                                 currYear, ipmZones, ipmZoneNums, ocAdderMin, ocAdderMax, cellsToZones,
                                 ptCurtailed | ptCurtailedRegs,
                                 statePolys)
    genFleet = selectAndMarkUnitsRetiredByCE(genFleet, genFleetForCE, retirementCFCutoff, capacExpModel, currYear,
                                             capacExpGenByGens, capacExpRetiredUnitsByCE, scaleMWtoGW, hoursForCE,
                                             planningReserveZonal, endYear, capacExpRetiredUnitsByAge, demandCEZonal,
                                             hourlyWindGenCEZonal, hourlySolarGenCEZonal, newWindCFsCEZonal,
                                             newSolarCFsCEZonal, ptEligRetCF)

    # removes all retired units; [] is dummy list b/c not adding ret age units to list
    genFleetNoRetiredUnits = createFleetForCurrentCELoop(genFleet, currYear, [], dataRoot, scenario)

    writeCEInfoToCSVs(capacExpBuilds, capacExpGenByGens, capacExpRetiredUnitsByCE,
                      capacExpRetiredUnitsByAge, resultsDir, currYear)
    write2dListToCSV(genFleet, os.path.join(resultsDir, 'genFleetAfterCE' + str(currYear) + '.csv'))

    # Write gen fleet for UC to special folder for ease of transfer
    ceUCDir = os.path.join(resultsDirOrig, 'CEtoUC')

    if not os.path.exists(ceUCDir): os.makedirs(ceUCDir)
    write2dListToCSV(genFleetNoRetiredUnits, os.path.join(ceUCDir, 'genFleetCEtoUC' + str(currYear) + '.csv'))

    return genFleet, genFleetNoRetiredUnits, genFleetForCE, capacExpModel, hoursForCE


########### IMPORT HOURLY THERMAL CURTAILMENTS BY GENERATOR ####################
# Output: dict mapping each gen to 2d list of datetime for year of run to hourly
# net capacity curtailments (MW)
def importHourlyThermalCurtailments(genFleet, currYear, modelName, ptCurtailed, ptCurtailedRegs,
                                    resultsDir, incCurtailments, incRegs, envRegMaxT, dataRoot, coolDesignT):
    (outputHeaders, locPrecision, tempAndSpatFilename, nsegFilename,
     rbmRootDir, rbmDataDir, rbmOutputDir, numCellsToProcess, envRegMaxT) = setThermalCurtailmentParameters()
    (hrlyCurtailmentsAllGensInTgtYr, genToCellLatLongsList,
     hrlyCurtailmentsList) = determineHrlyCurtailmentsForExistingGens(locPrecision,
                                                                      genFleet, rbmOutputDir, currYear, modelName,
                                                                      ptCurtailed, ptCurtailedRegs,
                                                                      incCurtailments, incRegs, envRegMaxT, dataRoot,
                                                                      coolDesignT, resultsDir)
    write2dListToCSV(genToCellLatLongsList, os.path.join(resultsDir,
                                                         'mapGensToCells' + modelName + str(currYear) + '.csv'))
    write2dListToCSV(hrlyCurtailmentsList, os.path.join(resultsDir,
                                                        'curtailmentsHourly' + str(currYear) + '.csv'))
    return hrlyCurtailmentsAllGensInTgtYr


########### CALL CAPACITY EXPANSION ############################################
def callCapacityExpansion(genFleetForCE, hourlyCapacsCE, hourlyCurtailedTechCapacsCE,
                          hourlyWindGenCEZonal, hourlySolarGenCEZonal, demandCEZonal, newTechsCE, planningReserveZonal,
                          discountRate, hoursForCE, newWindCFsCEZonal, newSolarCFsCEZonal, scaleMWtoGW,
                          scaleDollarsToThousands, currCo2Cap, capacExpFilename, seasonDemandWeights, repHrsBySeason,
                          specialHrs, peakDemandHourZonal, scaleLbToShortTon, dataRoot, maxAddedZonalCapacPerTech,
                          xsedeRun, ipmZones, lineList, lineCapacs, cellsToZones, ipmZoneNums, ptCurtailedAll,
                          hydroPotPerSeason, phEff, phMaxSoc, phInitSoc):
    currDir = os.getcwd()

    gamsFileDir = os.path.join(dataRoot, 'GAMS')
    gamsSysDir = 'C:\\GAMS\\win64\\24.7'

    if not xsedeRun:
        ws = GamsWorkspace(working_directory=gamsFileDir, system_directory=gamsSysDir)
    elif xsedeRun:
        ws = GamsWorkspace(working_directory=gamsFileDir)

    db = ws.add_database()
    # Add sets and parameters to database
    (genSet, genSymbols, hourSet, hourSymbols, techSet, techSymbols, techCurtailedSet, techCurtailedSymbols,
     renewTechSet, renewTechSymbols, techNotCurtailedSymbols, hydroGenSet, hydroGenSymbols, zoneSet, zoneSymbols,
     lineSet, cellSet, peakHourSet, peakHrSymbols, pumpHydroGenSet,
     pumpHydroGenSymbols) = addSetsToDatabaseCE(db, genFleetForCE,
                                                hoursForCE, newTechsCE, repHrsBySeason, specialHrs, peakDemandHourZonal,
                                                ipmZoneNums,
                                                lineList, cellsToZones, ptCurtailedAll)
    addParametersToDatabaseCE(db, hourlyCapacsCE, hourlyWindGenCEZonal, hourlySolarGenCEZonal, demandCEZonal,
                              newTechsCE, genFleetForCE, genSet, genSymbols, hydroGenSet, hydroGenSymbols, hourSet,
                              hourSymbols, techSet,
                              techSymbols, techCurtailedSet, techCurtailedSymbols, renewTechSet, renewTechSymbols,
                              planningReserveZonal,
                              discountRate, newWindCFsCEZonal, newSolarCFsCEZonal, scaleMWtoGW, scaleDollarsToThousands,
                              currCo2Cap, seasonDemandWeights, hourlyCurtailedTechCapacsCE, scaleLbToShortTon,
                              maxAddedZonalCapacPerTech,
                              zoneSet, zoneSymbols, lineSet, ipmZones, ipmZoneNums, lineList, lineCapacs, cellsToZones,
                              cellSet, peakDemandHourZonal, peakHourSet, peakHrSymbols, hydroPotPerSeason,
                              ptCurtailedAll,
                              pumpHydroGenSet, pumpHydroGenSymbols, phEff, phMaxSoc, phInitSoc)
    # Load GAMS model
    capacExpFile = capacExpFilename
    capacExpModel = ws.add_job_from_file(capacExpFile)
    # Run GAMS model
    opt = GamsOptions(ws)
    opt.defines['gdxincname'] = db.name
    capacExpModel.run(opt, databases=db)
    ms = capacExpModel.out_db['pModelstat'].find_record().value
    ss = capacExpModel.out_db['pSolvestat'].find_record().value
    if int(ms) != 8 or int(ms) != 1 or int(ss) != 1:
        print('Modelstat & solvestat:', ms, ' & ', ss,
              '(should be 8 (integer solution model found) or 1 (optimal solution found) & 1)')
    return capacExpModel, ms, ss, techCurtailedSymbols, renewTechSymbols, techNotCurtailedSymbols, pumpHydroGenSymbols


################################### ADD SETS
def addSetsToDatabaseCE(db, genFleetForCE, hoursForCE, newTechsCE, repHrsBySeason, specialHrs,
                        peakDemandHourZonal, ipmZoneNums, lineList, cellsToZones, ptCurtailedAll):
    (genSet, genSymbols, hydroGenSet, hydroGenSymbols, pumpHydroGenSet,
     pumpHydroGenSymbols) = addGeneratorSets(db, genFleetForCE)
    (hourSet, hourSymbols) = addHourSet(db, hoursForCE)
    addHourSeasonSubsets(db, repHrsBySeason)
    addHourSpecialSubset(db, specialHrs)
    peakHourSet, peakHrSymbols = addPeakHourSubset(db, peakDemandHourZonal)
    (techSet, techSymbols, techCurtailedSet, techCurtailedSymbols, renewTechSet, renewTechSymbols,
     techNotCurtailedSet, techNotCurtailedSymbols) = addNewTechsSets(db, newTechsCE, ptCurtailedAll)
    zoneSet, zoneSymbols = addZoneSets(db, ipmZoneNums)  # create string for set IDs
    lineSet = addLineSets(db, lineList)
    cellSet = addCellSet(db, cellsToZones)
    return (genSet, genSymbols, hourSet, hourSymbols, techSet, techSymbols, techCurtailedSet,
            techCurtailedSymbols, renewTechSet, renewTechSymbols, techNotCurtailedSymbols,
            hydroGenSet, hydroGenSymbols, zoneSet, zoneSymbols, lineSet, cellSet, peakHourSet,
            peakHrSymbols, pumpHydroGenSet, pumpHydroGenSymbols)


################################### ADD PARAMETERS
def addParametersToDatabaseCE(db, hourlyCapacsCE, hourlyWindGenCEZonal, hourlySolarGenCEZonal, demandCEZonal,
                              newTechsCE, genFleetForCE, genSet, genSymbols, hydroGenSet, hydroGenSymbols, hourSet,
                              hourSymbols, techSet,
                              techSymbols, techCurtailedSet, techCurtailedSymbols, renewTechSet, renewTechSymbols,
                              planningReserveZonal,
                              discountRate, newWindCFsCEZonal, newSolarCFsCEZonal, scaleMWtoGW, scaleDollarsToThousands,
                              currCo2Cap, seasonDemandWeights, hourlyCurtailedTechCapacsCE, scaleLbToShortTon,
                              maxAddedZonalCapacPerTech,
                              zoneSet, zoneSymbols, lineSet, ipmZones, ipmZoneNums, lineList, lineCapacs, cellsToZones,
                              cellSet, peakDemandHourZonal, peakHourSet, peakHrSymbols, hydroPotPerSeason,
                              ptCurtailedAll,
                              pumpHydroGenSet, pumpHydroGenSymbols, phEff, phMaxSoc, phInitSoc):
    addTechParams(db, newTechsCE, techSet, techSymbols, hourSet, hourSymbols,
                  scaleMWtoGW, scaleDollarsToThousands, scaleLbToShortTon, ptCurtailedAll)
    addEguParams(db, genFleetForCE, genSet, genSymbols, ipmZones, ipmZoneNums, scaleLbToShortTon, scaleMWtoGW)
    addPumpHydroParams(db, genFleetForCE, phEff, phMaxSoc, phInitSoc, pumpHydroGenSet, pumpHydroGenSymbols, scaleMWtoGW)
    addDiscountRateParam(db, discountRate)
    addCppEmissionsCap(db, currCo2Cap)
    addSeasonDemandWeights(db, seasonDemandWeights)
    addMaxNumNewBuilds(db, newTechsCE, zoneSet, ipmZones, ipmZoneNums, techSet, maxAddedZonalCapacPerTech,
                       ptCurtailedAll)
    # Add hydro max gen per time block
    addHydroMaxGenPerSeason(db, hydroGenSet, hydroGenSymbols, hydroPotPerSeason, scaleMWtoGW)
    # Add zonal demand, planning margin, and RE generation
    addDemandParam(db, demandCEZonal, zoneSet, hourSet, hourSymbols, ipmZones, ipmZoneNums, scaleMWtoGW)
    addPlanningReserveParam(db, planningReserveZonal, ipmZones, ipmZoneNums, zoneSet, zoneSymbols, scaleMWtoGW)
    addPeakHourToZoneParam(db, peakDemandHourZonal, peakHourSet, peakHrSymbols, ipmZones, ipmZoneNums)
    addExistingRenewableMaxGenParams(db, zoneSet, ipmZones, ipmZoneNums, hourSet, hourSymbols,
                                     hourlySolarGenCEZonal, hourlyWindGenCEZonal, scaleMWtoGW)
    addRenewTechCFParams(db, renewTechSet, renewTechSymbols, zoneSet, hourSet, hourSymbols,
                         newWindCFsCEZonal, newSolarCFsCEZonal, ipmZones, ipmZoneNums)
    # Add hourly capacity for curtailed techs
    addTechCurtailedHourlyCapac(db, hourlyCurtailedTechCapacsCE, cellSet, techCurtailedSet, hourSet, hourSymbols,
                                scaleMWtoGW)
    # Add tech heat rate and op cost params that reflect degradation
    addEguHourlyParams(db, hourlyCapacsCE, genSet, hourSet, hourSymbols, scaleMWtoGW)
    addEguOpCostParam(db, genFleetForCE, genSet, genSymbols, scaleLbToShortTon, scaleMWtoGW, scaleDollarsToThousands)
    # Add zone & line parameters
    addLineCapacs(db, lineCapacs, lineSet, lineList, scaleMWtoGW)
    addLineSourceAndSink(db, lineSet, lineList, ipmZones, ipmZoneNums)
    addCellsToZones(db, cellSet, cellsToZones, ipmZones, ipmZoneNums)


################################################################################
####### RUN UNIT COMMITMENT ####################################################
################################################################################
def runUnitCommitment(genFleet, demandScaled, startYear, ucYear, fuelPricesTimeSeries,
                      states, statesAbbrev, scaleMWtoGW, scaleDollarsToThousands, currCo2Cap, calculateCO2Price,
                      scaleLbToShortTon, daysForUC, daysOpt, daysLA, tzAnalysis, projectName, dataRoot,
                      resultsDir, ocAdderMin, ocAdderMax, windGenDataYr, regLoadFrac, contLoadFrac,
                      regErrorPercentile, flexErrorPercentile, rrToRegTime, rrToFlexTime, rrToContTime,
                      regUpCostCoeffsUC, xsedeRun, runCE, scenario, ucFilename, ipmZones, incCurtailments):
    resultsDir = os.path.join(resultsDir, 'UC')
    if not os.path.exists(resultsDir): os.makedirs(resultsDir)
    print('Entering UC loop for year ' + str(ucYear))
    fleetUC = copy.deepcopy(genFleet)

    (startWindCapacForCFs, startSolarCapacForCFs) = (0, 0)
    (windCFs, windCfsDtHr, windCfsDtSubhr, windIdAndCapac, solarCFs, solarCfsDtHr, solarCfsDtSubhr,
     solarFilenameAndCapac) = getRenewableCFs(fleetUC, startWindCapacForCFs, startSolarCapacForCFs,
                                              states, statesAbbrev, tzAnalysis, projectName, dataRoot, windGenDataYr)
    write2dListToCSV([demandScaled], os.path.join(resultsDir, 'demandUC' + str(ucYear) + '.csv'))
    write2dListToCSV(windCFs, os.path.join(resultsDir, 'windCFsUC' + str(ucYear) + '.csv'))
    write2dListToCSV(windCfsDtHr, os.path.join(resultsDir, 'windCFsDtUC' + str(ucYear) + '.csv'))
    write2dListToCSV(windCfsDtSubhr, os.path.join(resultsDir, 'windCFsDtSubhrUC' + str(ucYear) + '.csv'))
    write2dListToCSV(windIdAndCapac, os.path.join(resultsDir, 'windIdAndCapacUC' + str(ucYear) + '.csv'))
    write2dListToCSV(solarCFs, os.path.join(resultsDir, 'solarCFsUC' + str(ucYear) + '.csv'))
    write2dListToCSV(solarCfsDtHr, os.path.join(resultsDir, 'solarCFsDtUC' + str(ucYear) + '.csv'))
    write2dListToCSV(solarCfsDtSubhr, os.path.join(resultsDir, 'solarCFsDtSubhrUC' + str(ucYear) + '.csv'))
    write2dListToCSV(solarFilenameAndCapac, os.path.join(resultsDir, 'solarIdAndCapacUC' + str(ucYear) + '.csv'))
    (contResHourly, regUpHourly, regDownHourly, flexResHourly, allRes, regUpWind,
     regDownWind, regUpSolar, regDownSolar, flexWind, flexSolar) = calcWWSISReserves(windCfsDtHr,
                                                                                     windCfsDtSubhr, windIdAndCapac,
                                                                                     solarCfsDtHr, solarCfsDtSubhr,
                                                                                     solarFilenameAndCapac,
                                                                                     demandScaled,
                                                                                     regLoadFrac, contLoadFrac,
                                                                                     regErrorPercentile,
                                                                                     flexErrorPercentile)
    write2dListToCSV(allRes, os.path.join(resultsDir, 'reservesUC' + str(ucYear) + '.csv'))
    write2dListToCSV([contResHourly], os.path.join(resultsDir, 'reservesContUC' + str(ucYear) + '.csv'))
    write2dListToCSV([regUpHourly], os.path.join(resultsDir, 'reservesRegUpUC' + str(ucYear) + '.csv'))
    write2dListToCSV([regDownHourly], os.path.join(resultsDir, 'reservesRegDownUC' + str(ucYear) + '.csv'))
    write2dListToCSV([flexResHourly], os.path.join(resultsDir, 'reservesFlexUC' + str(ucYear) + '.csv'))
    (netDemand, hourlyWindGen, hourlySolarGen) = getNetDemand(demandScaled, windCFs, windIdAndCapac, solarCFs,
                                                              solarFilenameAndCapac, ucYear, 'UC', resultsDir)
    write2dListToCSV([[val] for val in hourlyWindGen], os.path.join(resultsDir, 'windGenUC' + str(ucYear) + '.csv'))
    write2dListToCSV([[val] for val in hourlySolarGen], os.path.join(resultsDir, 'solarGenUC' + str(ucYear) + '.csv'))
    updateFuelPricesExistingGens(fleetUC, ucYear, fuelPricesTimeSeries)
    combineWindAndSolarToSinglePlant(fleetUC, ipmZones, dataRoot)
    dailyCurtailmentsAllGensInTgtYr = importDailyThermalCurtailments(fleetUC, ucYear,
                                                                     'UC', ptCurtailed, resultsDir,
                                                                     incCurtailments)  # dict w/ gen IDs, then 2d list
    (hourlyCapacsAllGens, hourlyCapacsAllGensList) = calculateHourlyCapacsWithCurtailments(fleetUC,
                                                                                           dailyCurtailmentsAllGensInTgtYr,
                                                                                           ucYear)
    write2dListToCSV(hourlyCapacsAllGensList, 'curtailedHourlyCapacsUC' + str(ucYear) + '.csv')
    if calculateCO2Price:
        co2Price = convertCo2CapToPrice(fleetUC, hourlyWindGen, hourlySolarGen, demandScaled, currCo2Cap,
                                        scaleMWtoGW, scaleDollarsToThousands, scaleLbToShortTon, dataRoot)
    else:
        co2Price = 0
    print('CO2 price:', co2Price, '$/ton')
    write2dListToCSV([[co2Price]], os.path.join(resultsDir, 'co2PriceUC' + str(ucYear) + '.csv'))
    write2dListToCSV(fleetUC, os.path.join(resultsDir, 'genFleetUC' + str(ucYear) + '.csv'))
    # Set up results lists
    ucResultsByDay = []  # list of UC GAMS models
    (genByPlant, regUpByPlant, flexByPlant, contByPlant, turnonByPlant, turnoffByPlant, regDownByPlant,
     onOffByPlant, genToRow, hourToColPlant) = setupHourlyResultsByPlant(daysForUC, fleetUC)
    (sysResults, resultToRow, hourToColSys) = setupHourlySystemResults(daysForUC)
    msAndSs = [['day', 'ms', 'ss']]  # store modelstat & solvestat from GAMS
    for dayIdx in range(0, len(daysForUC), daysOpt):
        day = daysForUC[dayIdx]
        (demandUC, hourlyWindGenUC, hourlySolarGenUC, hoursForUC) = getDemandAndREGenForUC(day, daysOpt, daysLA,
                                                                                           demandScaled, hourlyWindGen,
                                                                                           hourlySolarGen)
        (regUpUC, regDownUC, flexUC, contUC) = getResForUC(day, daysOpt, daysLA, regUpHourly, regDownHourly,
                                                           flexResHourly, contResHourly)
        hourlyCapacsUC = getHourlyCapacitiesForDays(fleetUC, hourlyCapacsAllGens, hoursForUC)
        if daysForUC[0] == day:  # first day, therefore no initial conditions defined. MW energy values
            (onOffInitial, genAboveMinInitial, mdtCarriedInitial) = setInitCondsFirstUC(fleetUC)
        else:  # other days. MW energy values
            (onOffInitial, genAboveMinInitial, mdtCarriedInitial) = setInitCondsPerPriorUC(ucModel, fleetUC,
                                                                                           hoursForUC, daysOpt, daysLA,
                                                                                           scaleMWtoGW)
        t0 = time.time()
        ucModel, ms, ss = callUnitCommitment(fleetUC, ucFilename, hourlyCapacsUC, hourlyWindGenUC,
                                             hourlySolarGenUC, demandUC, hoursForUC, onOffInitial, genAboveMinInitial,
                                             mdtCarriedInitial,
                                             scaleMWtoGW, scaleDollarsToThousands, co2Price, scaleLbToShortTon, regUpUC,
                                             regDownUC, dataRoot, rrToRegTime, rrToFlexTime, rrToContTime, flexUC,
                                             contUC,
                                             xsedeRun)
        print('Time (secs) for UC day ' + str(day) + ': ' + str(time.time() - t0))
        ucResultsByDay.append((day, ucModel))  # just saves GAMS model
        saveHourlyResultsByPlant(genByPlant, regUpByPlant, regDownByPlant, flexByPlant, contByPlant,
                                 turnonByPlant, turnoffByPlant, onOffByPlant, genToRow, hourToColPlant, ucModel, day,
                                 daysOpt)
        saveHourlySystemResults(sysResults, resultToRow, hourToColSys, ucModel, day, daysOpt)
        msAndSs.append([day, ms, ss])
    writeHourlyResultsByPlant(genByPlant, regUpByPlant, regDownByPlant, flexByPlant, contByPlant,
                              turnonByPlant, turnoffByPlant, onOffByPlant, resultsDir, ucYear, 'UC', 'Plant')
    write2dListToCSV(sysResults, os.path.join(resultsDir, 'systemResultsUC' + str(ucYear) + '.csv'))
    write2dListToCSV(msAndSs, os.path.join(resultsDir, 'msAndSsUC' + str(ucYear) + '.csv'))
    return (ucResultsByDay, genByPlant)


########### TESTING PURPOSES ONLY ####################################
# REMOVE THIS FUNCTION - JUST INCLUDED FOR TESTING PURPOSES TO EXPEDITE UC
def removeHydroForTesting(fleetUC):
    fuelTypeCol = fleetUC[0].index('Modeled Fuels')
    idxs = []
    for idx in range(len(fleetUC) - 1, 0, -1):
        if fleetUC[idx][fuelTypeCol] == 'Hydro': fleetUC.pop(idx)


########### RUN UNIT COMMITMENT MODEL ##########################################
def callUnitCommitment(fleetUC, ucFilename, hourlyCapacsUC, hourlyWindGenUC,
                       hourlySolarGenUC, demandUC, hoursForUC, onOffInitial, genAboveMinInitial, mdtCarriedInitial,
                       scaleMWtoGW, scaleDollarsToThousands, co2Price, scaleLbToShortTon, regUpUC,
                       regDownUC, dataRoot, rrToRegTime, rrToFlexTime, rrToContTime, flexUC, contUC, xsedeRun):
    currDir = os.getcwd()

    gamsFileDir = os.path.join(dataRoot, 'GAMS')
    gamsSysDir = 'C:\\GAMS\\win64\\24.7'

    if xsedeRun == False:
        wsUC = GamsWorkspace(working_directory=gamsFileDir, system_directory=gamsSysDir)
    elif xsedeRun == True:
        wsUC = GamsWorkspace(working_directory=gamsFileDir)

    dbUC = wsUC.add_database()
    # Add sets and parameters to database
    cnse = 10000
    (genSet, genSymbols, hourSet, hourSymbols, hydroGenSet,
     hydroGenSymbols) = addSetsToDatabaseUC(dbUC, fleetUC, hoursForUC)
    addParametersToDatabaseUC(dbUC, hourlyCapacsUC, hourlyWindGenUC, hourlySolarGenUC, demandUC, fleetUC,
                              genSet, genSymbols, hydroGenSet, hydroGenSymbols, hourSet, hourSymbols, cnse,
                              rrToFlexTime,
                              onOffInitial, genAboveMinInitial, mdtCarriedInitial, scaleMWtoGW, scaleDollarsToThousands,
                              co2Price, hoursForUC, scaleLbToShortTon, rrToRegTime, regUpUC, regDownUC, rrToContTime,
                              flexUC, contUC, lineSet, lineList, lineCapacs)
    # Load and run GAMS model
    ucModel = wsUC.add_job_from_file(ucFilename)
    optUC = GamsOptions(wsUC)
    optUC.defines['gdxincname'] = dbUC.name
    ucModel.run(optUC, databases=dbUC)
    ms, ss = ucModel.out_db['pModelstat'].find_record().value, ucModel.out_db['pSolvestat'].find_record().value
    if int(ms) != 8 or int(ss) != 1: print('Modelstat & solvestat:', ms, ' & ', ss, ' (should be 8 and 1)')
    return ucModel, ms, ss


################################### ADD SETS
def addSetsToDatabaseUC(db, fleetUC, hoursForUC):
    (genSet, genSymbols, windGenSet, windGenSymbols, solarGenSet, solarGenSymbols,
     hydroGenSet, hydroGenSymbols) = addGeneratorSets(db, fleetUC)
    (hourSet, hourSymbols) = addHourSet(db, hoursForUC)
    return (genSet, genSymbols, hydroGenSet, hydroGenSymbols, hourSet, hourSymbols)


################################### ADD PARAMETERS
def addParametersToDatabaseUC(db, hourlyCapacsUC, hourlyWindGenUC, hourlySolarGenUC, demandUC, fleetUC,
                              genSet, genSymbols, hydroGenSet, hydroGenSymbols, hourSet, hourSymbols, cnse,
                              rrToFlexTime,
                              onOffInitial, genAboveMinInitial, mdtCarriedInitial, scaleMWtoGW, scaleDollarsToThousands,
                              co2Price, hoursForUC, scaleLbToShortTon, rrToRegTime, regUpUC, regDownUC, rrToContTime,
                              flexUC, contUC, lineSet, lineList, lineCapacs):
    addDemandParam(db, demandUC, hourSet, hourSymbols, scaleMWtoGW)
    addEguParams(db, fleetUC, genSet, genSymbols, scaleLbToShortTon)
    addEguHourlyParams(db, hourlyCapacsCE, hourlyHrsCE, genSet, hourSet, hourSymbols, scaleMWtoGW)
    addEguOpCostParam(db, fleetUC, genSet, genSymbols, scaleLbToShortTon, scaleMWtoGW, scaleDollarsToThousands,
                      co2Price)
    addEguUCParams(db, fleetUC, genSet, genSymbols, scaleMWtoGW, scaleDollarsToThousands)
    addEguInitialConditions(db, genSet, genSymbols, fleetUC, onOffInitial, genAboveMinInitial,
                            mdtCarriedInitial, scaleMWtoGW)
    addEguEligibleToProvideRes(db, fleetUC, genSet, genSymbols, stoMarket)
    addExistingRenewableMaxGenParams(db, hourSet, hourSymbols, hourlySolarGenUC,
                                     hourlyWindGenUC, scaleMWtoGW)
    addRegReserveParameters(db, regUpUC, regDownUC, rrToRegTime, hourSet, hourSymbols, scaleMWtoGW, 'UC')
    addEguRegCostParam(db, fleetUC, genSet, genSymbols, scaleMWtoGW, scaleDollarsToThousands)
    addFlexReserveParameters(db, flexUC, rrToFlexTime, hourSet, hourSymbols, scaleMWtoGW, 'UC')
    addContReserveParameters(db, contUC, rrToContTime, hourSet, hourSymbols, scaleMWtoGW)
    addCostNonservedEnergy(db, cnse, scaleMWtoGW, scaleDollarsToThousands)
    addCo2Price(db, co2Price, scaleDollarsToThousands)
    # Add transmission line constraints
    addLineCapacs(db, lineCapacs, lineSet, lineList, scaleMWtoGW)
    addLineSourceAndSink(db, lineSet, lineList)


################################################################################
################################################################################
################################################################################

masterFunction()
