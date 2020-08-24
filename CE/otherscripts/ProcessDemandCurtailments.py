"""
Francisco February 2019

This script runs a pre-processing of the demand simulations and curtailment simulations for all
GCM scenarios available and all years.

It saves the outputs to pandas data frames

"""
#
# ----- import custom functions and classes from CE/UC -----
from RIPSMasterScript17Nov2017 import getInitialFleetAndDemand
from ForecastDemandWithRegression import *
from ModifyGeneratorCapacityWithWaterTData import *
from ModifyNewTechCapacityWithWaterTData import *
from LoadEligibleCellWaterTs import *
from ImportNewTechs import *
from AuxFuncs import *
from UpdateFuelPriceFuncs import *
from Parameters import *
#
# ---- import python packages -----
import multiprocessing as mp
import pandas as pd
import numpy as np
import pickle as pk
import netCDF4 as nc
import re
import gc


def wrapper_forecast_demand(list_args):
    curr_year, gen_param, curt_param = list_args

    zonalDemandProfile, zonalTempDfs = forecastZonalDemandWithReg(curr_year, gen_param, curt_param)

    gc.collect()

    return zonalDemandProfile


def wrapper_new_tech_curtailments(list_args):
    """ Worker function for using multiprocessing with function calculateGeneratorCurtailments

    :param list_args: list with arguments for function
    :return:
    """

    # create data frame of 2d list with new techs
    df_newtechs = pd.DataFrame(list_args[1][1:], columns=list_args[1][0])

    curtailments_out = calculateTechsCurtailments(cellWaterTsForNewTechs=list_args[0], newTechsCE=list_args[1],
                                                  currYear=list_args[2], genparam=list_args[3],
                                                  curtailparam=list_args[4], gcm=list_args[5])

    # convert to a data frame and summarize values to Mean, Min and Summer mean
    df = pd.DataFrame(curtailments_out)
    df['hour'] = df.index
    df = pd.melt(df, id_vars=['hour'], var_name=['type', 'cell'])
    df['season'] = df['hour'].apply(lambda x: 'Summer' if 3624 <= x <= 5831 else 'NotSummer')

    del curtailments_out
    gc.collect()

    # downcast strings to categorical to save space
    df['type'] = df['type'].astype('category')
    df['cell'] = df['cell'].astype('category')
    df['season'] = df['season'].astype('category')

    gc.collect()

    df2 = df.groupby(['type', 'cell'])['value'].agg([np.mean, np.min]).reset_index()

    df3 = df.groupby(['type', 'cell', 'season'])['value'].agg([np.mean]).\
        reset_index().rename(columns={'mean': 'summer_mean'})
    df3 = df3[df3['season'] == 'Summer']
    del df3['season']

    result = pd.merge(df2, df3, how='left', on=['type', 'cell'])

    # add column with capacity info of each type
    result['type2'] = result['type'].apply(lambda x: x.split(sep='+')[0])
    df_newtechs_short = df_newtechs[['TechnologyType', 'Capacity(MW)']]
    result = pd.merge(result, df_newtechs_short, how='left', left_on=['type2'], right_on=['TechnologyType'])
    result = pd.DataFrame(result[['type', 'Capacity(MW)', 'cell', 'mean', 'amin', 'summer_mean']])

    # downcast strings to categorical to save space
    result['Capacity(MW)'] = result['Capacity(MW)'].astype('float')

    return result


def processCurtailmentsForNewTechs(eligibleCellWaterTs, newTechsCE, currYear, genparam, curtailparam, resultsDir,
                                   pbar=True):
    """ Returns dict of (plant+cooltype,cell):hrly tech curtailments

    :param eligibleCellWaterTs:
    :param newTechsCE:
    :param currYear:
    :param genparam:
    :param curtailparam:
    :param resultsDir:
    :param pbar:
    :return:
    """
    # Isolate water Ts to year of analysis
    eligibleCellWaterTsCurrYear = getWaterTsInCurrYear(currYear, eligibleCellWaterTs)
    cellWaterTsForNewTechs = selectCells(eligibleCellWaterTsCurrYear, genparam.cellNewTechCriteria,
                                         genparam.fipsToZones, genparam.fipsToPolys)

    # compile list of arguments to call multiprocessing
    args_list = [[cellWaterTsForNewTechs[gcm], newTechsCE, currYear, genparam, curtailparam, gcm]
                 for gcm in curtailparam.listgcms]

    ncores = min(len(curtailparam.listgcms), genparam.ncores_py)

    with mp.Pool(processes=ncores) as pool:
        list_curtailments = pool.map(wrapper_new_tech_curtailments, args_list)

    # combine all dfs into one single data frame
    df_out = pd.DataFrame()

    for i, df in enumerate(list_curtailments):
        # get gcm name
        gcm = curtailparam.listgcms[i]

        df_aux = pd.DataFrame(df)
        df_aux['gcm'] = gcm
        df_aux['gcm'] = pd.Categorical(df_aux['gcm'], categories=curtailparam.listgcms)

        df_aux = df_aux[['gcm', 'type', 'Capacity(MW)', 'cell', 'mean', 'amin', 'summer_mean']]

        df_out = df_out.append(df_aux, ignore_index=True)

    return df_out


def preProcessAll():

    cwd = os.getcwd()

    print('Loading parameters and setting up initial data')
    # Load parameters
    genparam = Generalparameters.Generalparameters()
    genparam.load(fname=os.path.join(cwd, 'generalparameters.txt'))

    reserveparam = Reserveparameters.Reserveparameters()
    reserveparam.load(fname=os.path.join(cwd, 'reserveparameters.txt'))

    curtailparam = Curtailmentparameters.Curtailmentparameters()
    curtailparam.load(fname=os.path.join(cwd, 'curtailmentparameters.txt'))

    # change cells eligible for new plants to 'all' (will compute for ALL grid cells that are not MASKED)
    genparam.cellsEligibleForNewPlants = 'maxflow'
    genparam.cellNewTechCriteria = 'all'
    genparam.resultsDir = '/sharedstorage/user/fralston/UWData/'

    genparam.startYear = 2045
    genparam.endYear = 2051
    genparam.ncores_py = 4

    genparam.referenceCase = False
    genparam.incCurtailments = True
    genparam.incRegs = True

    n_cells = 400

    list_rcps = ['rcp45', 'rcp85']

    listgcms_total=['bcc-csm1-1_rcp45', 'bcc-csm1-1_rcp85', 'bcc-csm1-1-m_rcp45', 'bcc-csm1-1-m_rcp85',
                    'BNU-ESM_rcp45', 'BNU-ESM_rcp85', 'CanESM2_rcp45', 'CanESM2_rcp85', 'CCSM4_rcp45',
                    'CCSM4_rcp85', 'CNRM-CM5_rcp45', 'CNRM-CM5_rcp85', 'CSIRO-Mk3-6-0_rcp45',
                    'CSIRO-Mk3-6-0_rcp85', 'GFDL-ESM2G_rcp45', 'GFDL-ESM2G_rcp85', 'GFDL-ESM2M_rcp45',
                    'GFDL-ESM2M_rcp85', 'HadGEM2-CC365_rcp45', 'HadGEM2-CC365_rcp85', 'HadGEM2-ES365_rcp45',
                    'HadGEM2-ES365_rcp85', 'inmcm4_rcp45', 'inmcm4_rcp85', 'IPSL-CM5A-LR_rcp45',
                    'IPSL-CM5A-LR_rcp85', 'IPSL-CM5A-MR_rcp45', 'IPSL-CM5A-MR_rcp85', 'IPSL-CM5B-LR_rcp45',
                    'IPSL-CM5B-LR_rcp85', 'MIROC-ESM_rcp45', 'MIROC-ESM_rcp85', 'MIROC-ESM-CHEM_rcp45',
                    'MIROC-ESM-CHEM_rcp85', 'MIROC5_rcp45', 'MIROC5_rcp85', 'MRI-CGCM3_rcp45',
                    'MRI-CGCM3_rcp85', 'NorESM1-M_rcp45', 'NorESM1-M_rcp85']

    # curtailparam.listgcms = ['NorESM1-M_rcp85', 'CSIRO-Mk3-6-0_rcp85']

    curtailparam.rbmRootDir ='/sharedstorage/user/fralston/UWData/'

    genFleet = getInitialFleetAndDemand(genparam, reserveparam)

    scp_command1 = ('rsync -avhe ssh --progress ' +
                    'pi@128.2.70.151:/home/pi/seagate/UWdata/forcing_maca_hourly_*_{0}_{1:4d}0101-{1:4d}1231.nc {2}')

    for currYear in range(genparam.startYear, genparam.endYear, genparam.yearStepCE):

        for rcp in list_rcps:

            curtailparam.listgcms = [g for g in listgcms_total if rcp in g]

            # first copy all GCM files from my external disk to resultsDir
            start_time = time.time()
            print('Downloading weather files...')
            os.system(scp_command1.format(rcp, currYear, genparam.resultsDir))
            print('Done' + str_elapsedtime(start_time))

            newTechsCE = getNewTechs(currYear, genparam, reserveparam)

            updateFuelPrices(genFleet, newTechsCE, currYear, genparam.fuelPricesTimeSeries)

            start_time = time.time()
            print('Computing hourly demand for each zone in year {0:4d}...'.format(currYear), flush=True)

            list_args = []
            for g in curtailparam.listgcms:
                cp = copy.deepcopy(curtailparam)
                cp.listgcms = [g]

                list_args.append([currYear, genparam, cp])

            with mp.Pool(processes=genparam.ncores_py) as pool:
                list_demand = pool.map(wrapper_forecast_demand, list_args)

            # list_demand is a list where each element is a dictionary with the projections of zonal hourly demand
            # for one GCM. Convert this list to a pandas DF with columns [GCM, zone, hour, demand]
            print('Converting results to panda data frame')
            dem_final = pd.DataFrame()
            for d in list_demand:
                # get gcm name
                gcm_name = list(d.keys())[0]

                # get dictionary with zonal demand and convert to pandas data frame with zones as column names
                dem_aux = pd.DataFrame(d[gcm_name])
                dem_aux['hour'] = dem_aux.index
                dem_aux = dem_aux.melt(id_vars=['hour'], var_name='zone', value_name='demand')
                dem_aux['gcm'] = gcm_name
                dem_aux = dem_aux[['gcm', 'zone', 'hour', 'demand']]

                # downcast to categorical
                dem_aux['gcm'] = pd.Categorical(dem_aux['gcm'], categories=curtailparam.listgcms)
                dem_aux['zone'] = dem_aux['zone'].astype('category')

                dem_final = dem_final.append(dem_aux, ignore_index=True)

            pd.to_pickle(dem_final, path=os.path.join(genparam.resultsDir, 'df_demand_{0}_{1}.pk'.format(rcp, currYear)))

            print('Done! Elapsed time: ' + str_elapsedtime(start_time))
            print('-----------------------------------------------------')

            start_time = time.time()
            print('Computing hourly curtailments for new generators...', flush=True)
            eligibleCellWaterTs = loadEligibleCellWaterTs(genFleet=None, currYear=currYear, genparam=genparam,
                                                          curtailparam=curtailparam, n_cells=n_cells)

            # hrlyCurtailmentsAllTechsInTgtYr is a dict of (tech,cell):[hourlycapac]
            df_curtail = processCurtailmentsForNewTechs(eligibleCellWaterTs, newTechsCE, currYear, genparam,
                                                        curtailparam, '')

            pd.to_pickle(df_curtail, path=os.path.join(genparam.resultsDir, 'curtailments_{0}_{1}.pk'.format(rcp, currYear)))

            print('Done! Elapsed time: ' + str_elapsedtime(start_time))
            print('-----------------------------------------------------')
            print()

            # delete GCM files to free up space
            print('Removing weather files...')
            os.system('rm {0}/forcing_maca_hourly*'.format(genparam.resultsDir))
            print('Done')
            print()
