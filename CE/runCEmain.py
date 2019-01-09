import os
import pandas as pd
from RIPSMasterScript17Nov2017 import *
from Parameters import *


def runCEmain(cwd=os.getcwd(), case=None):

    print('Loading parameters and setting up initial data')
    # Load parameters
    genparam = Generalparameters.Generalparameters()
    genparam.load(fname=os.path.join(cwd, 'generalparameters.txt'))

    reserveparam = Reserveparameters.Reserveparameters()
    reserveparam.load(fname=os.path.join(cwd, 'reserveparameters.txt'))

    curtailparam = Curtailmentparameters.Curtailmentparameters()
    curtailparam.load(fname=os.path.join(cwd, 'curtailmentparameters.txt'))

    if case is not None:

        case = int(case) - 1

        df = pd.read_csv(os.path.join(cwd, 'list_cases.csv'), comment='#', sep=',')

        genparam.resultsDir = df['resultsDir'].iloc[case].strip()
        genparam.referenceCase = df['referenceCase'].iloc[case].strip()
        curtailparam.listgcms = list(map(str.strip, df['listgcms'].iloc[0].split(';')))

        genparam.co2CapScenario = df['co2CapScenario'].iloc[case].lower().strip()

    # BASE LINE CASE
    if genparam.referenceCase:
        genparam.incCurtailments = False
        genparam.incRegs = False
        curtailparam.listgcms = [curtailparam.listgcms[0]]

    if not os.path.exists(genparam.resultsDir):
        os.mkdir(genparam.resultsDir)

    print()
    print('------------------------------------------')
    print('folder out: {}'.format(genparam.resultsDir))
    print('list gcms: {}'.format(curtailparam.listgcms))
    print('referenceCase: {}'.format(genparam.referenceCase))
    print('incCurtailments: {}'.format(genparam.incCurtailments))
    print('incRegs: {}'.format(genparam.incRegs))
    print('CO2 Cap Scenario: {}'.format(genparam.co2CapScenario))

    print('------------------------------------------')

    print()
    print('Parameters loaded! Calling CE masterFunction...')
    print()
    print()
    masterFunction(genparam, reserveparam, curtailparam)


if __name__ == "__main__":
    runCEmain()
