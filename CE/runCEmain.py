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

        df = pd.read_csv(fname=os.path.join(cwd, 'list_cases.csv'))

        genparam.resultsDir = df['resultsDir'].iloc[case]
        curtailparam.listgcms = [df['listgcms'].iloc[case]]

        print()
        print('------------------------------------------')
        print('folder out: {}'.format(genparam.resultsDir))
        print('list gcms: {}'.format(curtailparam.listgcms))
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
