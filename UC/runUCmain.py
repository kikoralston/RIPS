import os, sys
from memory_profiler import profile
from RIPSMasterScript17Nov2017 import *
from Parameters import *
import multiprocessing as mp

sys.stdout.flush()


def list_all_gcms(rcp=None):

    listgcms_total = ['bcc-csm1-1_rcp45', 'bcc-csm1-1_rcp85', 'bcc-csm1-1-m_rcp45', 'bcc-csm1-1-m_rcp85',
                      'BNU-ESM_rcp45', 'BNU-ESM_rcp85', 'CanESM2_rcp45', 'CanESM2_rcp85', 'CCSM4_rcp45',
                      'CCSM4_rcp85', 'CNRM-CM5_rcp45', 'CNRM-CM5_rcp85', 'CSIRO-Mk3-6-0_rcp45',
                      'CSIRO-Mk3-6-0_rcp85', 'GFDL-ESM2G_rcp45', 'GFDL-ESM2G_rcp85', 'GFDL-ESM2M_rcp45',
                      'GFDL-ESM2M_rcp85', 'HadGEM2-CC365_rcp45', 'HadGEM2-CC365_rcp85', 'HadGEM2-ES365_rcp45',
                      'HadGEM2-ES365_rcp85', 'inmcm4_rcp45', 'inmcm4_rcp85', 'IPSL-CM5A-LR_rcp45',
                      'IPSL-CM5A-LR_rcp85', 'IPSL-CM5A-MR_rcp45', 'IPSL-CM5A-MR_rcp85', 'IPSL-CM5B-LR_rcp45',
                      'IPSL-CM5B-LR_rcp85', 'MIROC-ESM_rcp45', 'MIROC-ESM_rcp85', 'MIROC-ESM-CHEM_rcp45',
                      'MIROC-ESM-CHEM_rcp85', 'MIROC5_rcp45', 'MIROC5_rcp85', 'MRI-CGCM3_rcp45',
                      'MRI-CGCM3_rcp85', 'NorESM1-M_rcp45', 'NorESM1-M_rcp85']

    if rcp is not None:
        listgcms_total = [g for g in listgcms_total if rcp in g]

    return listgcms_total

@profile
def runUCgcm(gcm, rcp, cwd=os.getcwd(), yearUC=2050):

    print('Loading parameters and setting up initial data')
    # Load parameters
    genparam = Generalparameters.Generalparameters()
    genparam.load(fname=os.path.join(cwd, 'generalparameters.txt'))

    reserveparam = Reserveparameters.Reserveparameters()
    reserveparam.load(fname=os.path.join(cwd, 'reserveparameters.txt'))

    curtailparam = Curtailmentparameters.Curtailmentparameters()
    curtailparam.load(fname=os.path.join(cwd, 'curtailmentparameters.txt'))

    curtailparam.listgcms = [gcm]
    genparam.rcp = rcp

    genparam.runCE = False
    genparam.runFirstUCYear = True
    genparam.runUC = True

    # BASE LINE CASE
    if genparam.referenceCase:
        genparam.incCurtailments = False
        genparam.incRegs = False
        genparam.selectCurtailDays = False

    genparam.startYear = yearUC
    genparam.endYear = yearUC+1

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
    runUCgcm('bcc-csm1-1_rcp45', 'rcp45')
