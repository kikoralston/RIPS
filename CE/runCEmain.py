import os, sys
import pandas as pd
from memory_profiler import profile
from RIPSMasterScript import *
import Parameters
sys.stdout.flush()


@profile
def runCEmain(cwd=os.getcwd(), case=None, runUC=False):

    print('Loading parameters and setting up initial data')
    # Load parameters
    genparam = Parameters.Generalparameters()
    genparam.load(fname=os.path.join(cwd, 'generalparameters.txt'))

    reserveparam = Parameters.Reserveparameters()
    reserveparam.load(fname=os.path.join(cwd, 'reserveparameters.txt'))

    curtailparam = Parameters.Curtailmentparameters()
    curtailparam.load(fname=os.path.join(cwd, 'curtailmentparameters.txt'))

    if case is not None:

        case = int(case) - 1

        df = pd.read_csv(os.path.join(cwd, 'list_cases.csv'), comment='#', sep=',', skipinitialspace=True,
                         dtype={'case': int, 'resultsDir': object, 'listgcms': object, 'referenceCase': bool,
                                'co2CapPercentage': float, 'rcp': object, 'useLineLimits': bool,
                                'gcmranking': object, 'ncores_py': int, 'ncores_gams': int,
                                'incDemand': bool, 'incHydroCC': bool, 'incThermoCC': bool},
                         keep_default_na=False)

        genparam.resultsDir = df['resultsDir'].iloc[case].strip()
        genparam.referenceCase = df['referenceCase'].iloc[case]
        curtailparam.listgcms = list(map(str.strip, df['listgcms'].iloc[case].split(';')))
        genparam.co2CapPercentage = df['co2CapPercentage'].iloc[case]
        genparam.rcp = df['rcp'].iloc[case].lower().strip()
        genparam.useLineLimits = df['useLineLimits'].iloc[case]

        genparam.gcmranking = list(map(str.strip, df['gcmranking'].iloc[case].split(';')))
        # convert ranking of gcms to expected format
        if genparam.gcmranking == ['']:
            # change empty list to None
            genparam.gcmranking = None
        else:
            # change string values to integer or float
            if genparam.referenceCase:
                genparam.gcmranking = list(map(float, genparam.gcmranking))
            else:
                genparam.gcmranking = list(map(int, genparam.gcmranking))

        genparam.ncores_py = df['ncores_py'].iloc[case]
        genparam.ncores_gams = df['ncores_gams'].iloc[case]

        genparam.incDemandCC = df['incDemandCC'].iloc[case]
        genparam.incHydroCC = df['incHydroCC'].iloc[case]

        genparam.incCurtailments = df['incThermoCC'].iloc[case]
        genparam.incRegs = df['incThermoCC'].iloc[case]
        genparam.selectCurtailDays = df['incThermoCC'].iloc[case]

    # read line constraints parameters *AFTER* updating `useLineLimits`
    genparam.setLines()

    # BASE LINE CASE
    if genparam.referenceCase:
        # reference case. No impact of climate change! Set all flags to false
        genparam.incCurtailments = False
        genparam.incRegs = False
        genparam.selectCurtailDays = False
        genparam.incDemandCC = False
        genparam.incHydroCC = False
        # if baseline case we just need a valid gcm for indexing the dictionaries (GCM data will not be used)
        curtailparam.listgcms = [curtailparam.listgcms[0]]

    if not os.path.exists(genparam.resultsDir):
        os.makedirs(genparam.resultsDir)

    print()
    print('------------------------------------------')

    print('folder out: {}'.format(genparam.resultsDir))
    print('list gcms: {}'.format(curtailparam.listgcms))
    print('referenceCase: {}'.format(genparam.referenceCase))
    print('selectCurtailDays: {}'.format(genparam.selectCurtailDays))
    print('incRegs: {}'.format(genparam.incRegs))
    print('incDemandCC: {}'.format(genparam.incDemandCC))
    print('incHydroCC: {}'.format(genparam.incHydroCC))
    print('CO2 Cap Scenario: {}'.format(genparam.co2CapPercentage))
    print('RCP: {}'.format(genparam.rcp))

    print('------------------------------------------')

    print()
    print('Parameters loaded! Calling CE masterFunction...')
    print()
    print()

    masterFunction(genparam, reserveparam, curtailparam)


if __name__ == "__main__":

    if len(sys.argv) == 2:
        runUC = (sys.argv[1] == 'True')
        runCEmain(runUC=runUC)
    elif len(sys.argv) == 3:
        runUC = (sys.argv[1] == 'True')
        runCEmain(case=sys.argv[2], runUC=runUC)
    else:
        runCEmain()
