import os, sys
import multiprocessing as mp
from memory_profiler import profile
from RIPSMasterScript17Nov2017 import *
from Parameters import *

#sys.stdout.flush()


def wrapper(arg_list):
    x = runUCsinglegcm(arg_list[0], arg_list[1], arg_list[2], arg_list[3])
    
    return x


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
def runUC(gcmIdxlist, pathresultroot, rcp, pathCEtoUC, cwd=os.getcwd(), yearUC=2050):

    resultsDir = os.path.join(pathresultroot)

    # redirect stdout to file
    #original_stdout = sys.stdout
    #sys.stdout = open(os.path.join(resultsDir, 'stdout.txt'), 'a')

    # copy CE results to results folder for this gcm (if not already copied)
    if not os.path.exists(os.path.join(resultsDir, 'CEtoUC')):
        shutil.copytree(pathCEtoUC, os.path.join(resultsDir, 'CEtoUC'), symlinks=False, ignore=None)

    print('Loading parameters and setting up initial data')
    # Load parameters
    genparam = Generalparameters.Generalparameters()
    genparam.load(fname=os.path.join(cwd, 'generalparameters.txt'))

    reserveparam = Reserveparameters.Reserveparameters()
    reserveparam.load(fname=os.path.join(cwd, 'reserveparameters.txt'))

    curtailparam = Curtailmentparameters.Curtailmentparameters()
    curtailparam.load(fname=os.path.join(cwd, 'curtailmentparameters.txt'))

    genparam.rcp = rcp

    if rcp == 'ref':
        genparam.referenceCase = True
        genparam.incCurtailments = False
        genparam.incRegs = False
        genparam.selectCurtailDays = False
        genparam.incDemandCC = False
        genparam.incHydroCC = False
    else:
        genparam.referenceCase = False
        genparam.incCurtailments = True
        genparam.incRegs = True
        genparam.selectCurtailDays = True
        genparam.incDemandCC = True
        genparam.incHydroCC = True

    genparam.runCE = False
    genparam.runFirstUCYear = True
    genparam.runUC = True
    
    genparam.ncores_py = 4
    genparam.gams = 1

    genparam.gcmranking = gcmIdxlist
    genparam.resultsDir = resultsDir

    genparam.startYear = yearUC
    genparam.endYear = yearUC + 1

    print()
    print('Parameters loaded! Calling CE masterFunction...')
    print()
    print()

    masterFunction(genparam, reserveparam, curtailparam)
    
    #sys.stdout = original_stdout

    return 0


if __name__ == "__main__":

    # runUCsinglegcm(11, '', 'rcp45')
        
    # first argument is range/list of gcms rankings. Two possible formats:
    #   - 'rankstart-rankend'. 'rankend' is included
    #   - 'rank1, rank2, rank3'
    if '-' in sys.argv[1]:
        start, end = sys.argv[1].split('-')
        gcmIdxlist = list(range(int(start), int(end)+1))
    else:
        gcmIdxlist = [float(x) for x in sys.argv[1].split(',')]

    rcp = sys.argv[2]
    pathresultroot = sys.argv[3]
    pathCEtoUC = sys.argv[4]

    if rcp != 'ref':
        gcmIdxlist = [int(x) for x in gcmIdxlist]

    print("----------------------------------------------")
    print(gcmIdxlist)
    print(rcp)
    print(pathresultroot)
    print(pathCEtoUC)

    # create pathresult for this UC run
    if not os.path.exists(pathresultroot):
        os.makedirs(pathresultroot)

    runUC(gcmIdxlist, pathresultroot, rcp, pathCEtoUC)
