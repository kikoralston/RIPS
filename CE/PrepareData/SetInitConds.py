# Michael Craig
# October 4, 2016
# Sets values for some initial condition parameters for CE and UC run,
# either for first run of entire year using assumed values (UC run) or based on last values in prior run.

import os
from GAMSUtil.GAMSAuxFuncs import *
from GAMSUtil.GAMSAddSetToDatabaseFuncs import isolateGenSymbols
from GAMSUtil.GAMSAddParamToDatabaseFuncs import getSocDict
import pickle as pk
from AuxFuncs import readCSVto2dList

def setInitCondsFirstUC(fleetUC, genparam):
    """SET INITIAL CONDITION PARAMETERS FOR FIRST UC RUN

    For first UC run of year. Assume all plants initially off w/ no carried MDT


    :param fleetUC: 2d list with gen fleet
    :param genparam: general parameters
    :return: 1d list of initial on/off, gen above min (MWh), & carried MDT values,
             Dict with initial state of charge for pumped hydro
    """
    onOffInitial = [0 for i in range(1, len(fleetUC))]
    genAboveMinInitial = [0 for i in range(1, len(fleetUC))] #MW
    mdtCarriedInitial = [0 for i in range(1, len(fleetUC))]

    # initial state of charge in day 1 is defined as a fraction of max SOC
    # max SOC is defined as phMaxSoc * Capacity
    pumpHydroGenSymbols = isolateGenSymbols(fleetUC, 'Pumped Storage')
    socDict = getSocDict(fleetUC, genparam.phMaxSoc, pumpHydroGenSymbols, genparam.scaleMWtoGW)
    initSocDict = {gen: socDict[gen] * genparam.phInitSoc for gen in socDict}

    return onOffInitial, genAboveMinInitial, mdtCarriedInitial, initSocDict
            

def setInitCondsPerPriorUC(ucModel, fleetUC, hoursForUC, daysOpt, daysLA, scaleMWtoGW):
    """SET INITIAL CONDITION PARAMETERS PER PRIOR UC RUN

    Set values for init cond params based on prior UC run

    :param ucModel: prior UC run results as GAMS object
    :param fleetUC: 2d list with gen fleet
    :param hoursForUC: hours for curr UC run
    :param daysOpt: num days for optim horiz (i.e., keep those results)
    :param daysLA: num days look ahead (need to skip these)
    :param scaleMWtoGW: scale factor MW to GW
    :return: 1d lists of initial on/off & gen above min (MWh) & carried MDT vals.
             Dict with initial state of charge for pumped hydro
    """

    # For genAboveMin & onOff, just need variable value in prior hour
    # lastHourSymbolPriorUCRun = 'h' + str((min(hoursForUC) - 1))

    # get list with hours from previous UC run and get last hour of optimization period in previous run
    hoursPreviousUC = [h.keys[0] for h in ucModel.out_db.get_set('h')]
    lastHourSymbolPriorUCRun = hoursPreviousUC[daysOpt*24 - 1]

    onOffDict = extract2dVarResultsIntoDict(ucModel, 'vOnoroff')
    onOffInitial = getInitCondValues(onOffDict, fleetUC, lastHourSymbolPriorUCRun)
    genAboveMinDict = extract2dVarResultsIntoDict(ucModel, 'vGenabovemin')
    genAboveMinInitial = getInitCondValues(genAboveMinDict, fleetUC, lastHourSymbolPriorUCRun, 1/scaleMWtoGW) #MW

    # For mdtCarriedInitial, get last turnoff decision, subtract # hours from then to end of last time period, and
    # subtract that from MDT.
    mdtCarriedInitial = getMdtCarriedInitial(onOffInitial, ucModel, fleetUC, hoursForUC, daysOpt, daysLA)

    # initial state of charge for pumped storage
    socDict = extract2dVarResultsIntoDict(ucModel, 'vSoc')
    initSocDict = {key[0]: socDict[key] for key in socDict if key[1] == lastHourSymbolPriorUCRun}

    return onOffInitial, genAboveMinInitial, mdtCarriedInitial, initSocDict


def getMdtCarriedInitial(onOffInitial, ucModel, fleetUC, hoursForUC, daysOpt, daysLA):
    """Determines carried MDT hours based on when unit turned off (if at all) in prior UC run

    :param onOffInitial: whether initially on/off in curr UC run
    :param ucModel: GAMS object with prior UC run results
    :param fleetUC: 2d list with gen fleet
    :param hoursForUC: hours included in curr UC run
    :param daysOpt: num days in optimization horizon (i.e., to keep)
    :param daysLA: num days included as LA
    :return: 1d list of carried MDT hours
    """
    mdtCarriedInitial = []
    fleetMDTCol = fleetUC[0].index('MinDownTime(hrs)')
    turnOffDict = extract2dVarResultsIntoDict(ucModel, 'vTurnoff')

    # lastHourPriorUCRun = min(hoursForUC) - 1
    # get list with hours from previous UC run and get last hour of optimization period in previous run (as an integer)
    hoursPreviousUC = [h.keys[0] for h in ucModel.out_db.get_set('h')]
    lastHourPriorUCRun = int(hoursPreviousUC[daysOpt*24 - 1].replace('h', ''))

    for rowNum in range(1, len(fleetUC)):
        if onOffInitial[rowNum-1] == 1:  # on @ start, so no MDT carried
            mdtCarriedInitial.append(0)
        else:
            genSymbol = createGenSymbol(fleetUC[rowNum], fleetUC[0])
            turnOff = 0
            for hr in range(lastHourPriorUCRun, lastHourPriorUCRun - (24*daysOpt) + 1, -1):
                if turnOffDict[(genSymbol, 'h' + str(hr))] == 1 and turnOff == 0:
                    turnOff = 1
                    genMDT = float(fleetUC[rowNum][fleetMDTCol])
                    # Hr that turn off counts toward MDT; therefore +1 to hr.
                    mdtCarriedInitial.append(max(0, genMDT - (lastHourPriorUCRun - hr + 1)))

            if turnOff == 0:
                mdtCarriedInitial.append(0)  # never turned off in last UC

    return mdtCarriedInitial


def getInitCondValues(initCondDict, fleetUC, lastHourSymbolPriorUCRun, *args):
    """Convert dict of output UC values into 1d list

    :param initCondDict: dictionary of UC output (genID:val)
    :param fleetUC: gen fleet
    :param lastHourSymbolPriorUCRun: last hour symbol of prior UC run
    :param args:
    :return: 1d list of init cond values for curr UC run. For energy values, outputs in MW.
    """
    initCondValues = []

    if len(args) > 0:
        scalar = args[0]
    else:
        scalar = 1

    for row in fleetUC[1:]:
        initCondValues.append(initCondDict[(createGenSymbol(row, fleetUC[0]), lastHourSymbolPriorUCRun)]*scalar)

    return initCondValues


def choose_gcms_for_ce(currYear, genparam, curtailparam):
    """Get names of GCMs for CE annual iteration

    In each iteration of the CE model, this function ranks GCM simulations according to annual peak demand values and
    chooses GCMs in the rank indexes predefined in `genparam`. Returns a new instance of `curtailparam` that includes
    only the chosen GCMs for the CE run

    :param currYear: (int) current year of CE simulation
    :param genparam: object of type :mod:`Generalparameters`
    :param curtailparam: object of type :mod:`Curtailmentparameters`
    :return: object of type :mod:`Curtailmentparameters` with only the desired gcms in the field `listgcms`
    """

    if not genparam.referenceCase:

        if genparam.gcmranking is not None:

            file_demand = 'df_demand_{1:}_{0:4d}.pk'.format(currYear, genparam.rcp)
            file_curtailment = 'curtailments_{1}_{0:4d}.pk'.format(currYear, genparam.rcp)

            # read demand file
            with open(os.path.join(curtailparam.rbmRootDir, file_demand), 'rb') as f:
                df_demand = pk.load(f)

            # compute total hourly system demand
            df_demand = df_demand.groupby(['gcm', 'hour']).agg({'demand': sum}).reset_index()

            df_demand = df_demand.groupby(['gcm']).agg({'demand': max}).reset_index().sort_values(
                by=['demand']).reset_index(drop=True)

            # get names of GCms in according to given ranking in genparam
            gcms_chosen = list(df_demand.iloc[genparam.gcmranking, ]['gcm'].astype('str'))
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
        curtparam_year.listgcms = ['ref{0:02d}'.format(i) for i, p in enumerate(genparam.gcmranking)]

    return curtparam_year


def get_init_conditions_ce(currYear, genparam, genFleet, genFleetNoRetiredUnits, genFleetPriorCE, priorCEout_db,
                           capacExpBuilds, capacExpGenByGens, capacExpRetiredUnitsByCE, capacExpRetiredUnitsByAge,
                           priorHoursCE):
    """Defines some initial conditions of CE run

    This function defines some of the initial conditions necessary the first iteration of the CE simulation. Because
    the CE simulation can crash in the middle of the multi-year sequential simulation, this function also reads data from
    a previous iteration in order to restart a CE simulation in the middle of the simulation horizon.

    For example, if the CE simulation has an original simulation horizon between 2015 and 2050 but it crashes during the
    2035 iteration, the user can restart it by setting the `genparam.coldStart = True` and `genparam.startYear = 2035`.
    This function will read the relevant results from the previous completed iteration (e.g. 2030) and will use them as
    starting condition for the 2035 iteration.

    :param currYear: (int) current year of CE simulation
    :param genparam: object of type :mod:`Generalparameters`
    :param genFleet: (2d list)
    :param genFleetNoRetiredUnits: (2d list)
    :param genFleetPriorCE:
    :param priorCEout_db: (gams database) gams database with
    :param capacExpBuilds:
    :param capacExpGenByGens:
    :param capacExpRetiredUnitsByCE:
    :param capacExpRetiredUnitsByAge:
    :param priorHoursCE:
    """
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
