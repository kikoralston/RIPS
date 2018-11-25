import os
import pandas as pd
from RIPSMasterScript17Nov2017 import *
from Parameters import *


def runCEmain(case=None):

    print('Loading parameters and setting up initial data')
    # Load parameters
    genparam = Generalparameters.Generalparameters()
    genparam.load(fname='./generalparameters.txt')

    reserveparam = Reserveparameters.Reserveparameters()
    reserveparam.load(fname='./reserveparameters.txt')

    curtailparam = Curtailmentparameters.Curtailmentparameters()
    curtailparam.load(fname='./curtailmentparameters.txt')

    if case is not None:

        case = int(case) - 1

        df = pd.read_csv('../../scripts/list_cases.csv')

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
