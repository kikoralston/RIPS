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
        curtailparam.listgcms = list(map(str.strip, df['listgcms'].iloc[0].split(';')))

        # BASE LINE CASE
        if curtailparam.listgcms[0] == '' and len(curtailparam.listgcms) == 1:
            # baseline case. change to 'na' or 'NA'
            curtailparam.listgcms = ['na']

        if curtailparam.listgcms[0] == 'na' and len(curtailparam.listgcms) == 1:
            genparam.incCurtailments = False
            genparam.incRegs = False

        print()
        print('------------------------------------------')
        print('folder out: {}'.format(genparam.resultsDir))
        print('list gcms: {}'.format(curtailparam.listgcms))
        print('incCurtailments: {}'.format(curtailparam.incCurtailments))
        print('incRegs: {}'.format(curtailparam.incRegs))
        print('------------------------------------------')

        if not os.path.exists(genparam.resultsDir):
            os.mkdir(genparam.resultsDir)

    print()
    print('Parameters loaded! Calling CE masterFunction...')
    print()
    print()
    masterFunction(genparam, reserveparam, curtailparam)


if __name__ == "__main__":
    runCEmain()
