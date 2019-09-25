# Michael Craig
# October 4, 2016
# Sets values for initial condition parameters for UC run, either for first run of entire year using assumed values
# or based on last values in prior run.

from GAMSAuxFuncs import *
from GAMSAddSetToDatabaseFuncs import isolateGenSymbols
from GAMSAddParamToDatabaseFuncs import getSocDict


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
