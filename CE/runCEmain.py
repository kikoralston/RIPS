import os, sys
import pandas as pd
from memory_profiler import profile
from RIPSMasterScript17Nov2017 import *
from Parameters import *
sys.stdout.flush()


@profile
def runCEmain(cwd=os.getcwd(), case=None, runUC=False):

    print('Loading parameters and setting up initial data')
    # Load parameters
    genparam = Generalparameters.Generalparameters()
    genparam.load(fname=os.path.join(cwd, 'generalparameters.txt'))

    reserveparam = Reserveparameters.Reserveparameters()
    reserveparam.load(fname=os.path.join(cwd, 'reserveparameters.txt'))

    curtailparam = Curtailmentparameters.Curtailmentparameters()
    curtailparam.load(fname=os.path.join(cwd, 'curtailmentparameters.txt'))

    if runUC:
        # temporary: run a single reference case UC simulation for 2015
        # I will add this to the list_cases file eventually

        genparam.referenceCase = True
        curtailparam.listgcms = [curtailparam.listgcms[0]]
        genparam.rcp = 'rcp45'

        genparam.incCurtailments = False
        genparam.incRegs = False
        genparam.selectCurtailDays = False

        genparam.runCE = False
        genparam.runFirstUCYear = True
        genparam.runUC = True

        genparam.startYear = 2015
        genparam.endYear = 2020

        if not os.path.exists(genparam.resultsDir):
            os.mkdir(genparam.resultsDir)

    else:
        if case is not None:

            case = int(case) - 1

            df = pd.read_csv(os.path.join(cwd, 'list_cases.csv'), comment='#', sep=',', skipinitialspace=True,
                             dtype={'case': int, 'resultsDir': object, 'listgcms': object, 'referenceCase': bool,
                                    'co2CapScenario': object, 'rcp': object, 'useLineLimits': bool,
                                    'gcmranking': object, 'ncores_py': int, 'ncores_gams': int},
                             keep_default_na=False)

            genparam.resultsDir = df['resultsDir'].iloc[case].strip()
            genparam.referenceCase = df['referenceCase'].iloc[case]
            curtailparam.listgcms = list(map(str.strip, df['listgcms'].iloc[case].split(';')))
            genparam.co2CapScenario = df['co2CapScenario'].iloc[case].lower().strip()
            genparam.rcp = df['rcp'].iloc[case].lower().strip()
            genparam.useLineLimits = df['useLineLimits'].iloc[case]

            genparam.gcmranking = list(map(str.strip, df['gcmranking'].iloc[case].split(';')))
            # convert ranking of gcms to expected format
            if genparam.gcmranking == ['']:
                # change empty list to None
                genparam.gcmranking = None
            else:
                # change string values to integer
                genparam.gcmranking = list(map(int, genparam.gcmranking))

            genparam.ncores_py = df['ncores_py'].iloc[case]
            genparam.ncores_gams = df['ncores_gams'].iloc[case]

        # BASE LINE CASE
        if genparam.referenceCase:
            genparam.incCurtailments = False
            genparam.incRegs = False
            genparam.selectCurtailDays = False
            # if baseline case we just need a valid gcm for indexing the dictionaries (GCM data will not be used)
            curtailparam.listgcms = [curtailparam.listgcms[0]]

        if not os.path.exists(genparam.resultsDir):
            os.mkdir(genparam.resultsDir)

        print()
        print('------------------------------------------')

        print('folder out: {}'.format(genparam.resultsDir))
        print('list gcms: {}'.format(curtailparam.listgcms))
        print('referenceCase: {}'.format(genparam.referenceCase))
        print('incCurtailments: {}'.format(genparam.incCurtailments))
        print('selectCurtailDays: {}'.format(genparam.selectCurtailDays))
        print('incRegs: {}'.format(genparam.incRegs))
        print('CO2 Cap Scenario: {}'.format(genparam.co2CapScenario))
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
