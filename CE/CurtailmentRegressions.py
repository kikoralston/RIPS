# Michael Craig
# August 7, 2017
# Script loads regression parameters from Aviva for curtailments,
# and

import copy
import numpy as np
import pandas as pd
import json
import os
from GAMSAuxFuncs import *
import matplotlib as plt


def setEnvRegCurtailments(coolType, capac, pt, ppDeltaT, metAndWaterData, maxT, streamAvailFrac, genparam):
    """
    Calculate environmental regulatory curtailments

    Calculate environmental regulatory curtailments for given time series of water temperatures and flows for a
    single generator. Curtailment values are multiplied by capacity to get available capacity, so value output
    by this function of 0 means full curtailment.

    Inputs: met and water data (PD Dataframe, hourly timescale), whether to model curtailments and curtailments from
    regulations, plant types that are curtailed, plant cooling type, plant capcity, plant plant type, plant efficiency,
    plant fuel type, plant cooling system designed change in T, , regulatory limit on stream T, and fraction of stream water available
    for cooling.

    :param coolType: (string) cooling type
    :param capac: capacity (MW)
    :param pt: plant type
    :param ppDeltaT:
    :param metAndWaterData: (pd data frame) water T data and water flow data (pd series of hourly data)
    :param maxT: regulatory stream temperature threshold
    :param streamAvailFrac: Maximum share of stream flow that can be used by power plant
    :param genparam: object of class General Parameters
    :return:
    """

    waterC = metAndWaterData['waterT'].values

    waterFlow = metAndWaterData['flow'].values

    # water flow is in cfs. Convert to gal/h
    waterFlow = waterFlow * (7.48052 * 3600)

    hrlyCurtsRegs = np.ones(metAndWaterData.shape[0])

    capacs = [np.concatenate((np.arange(0, float(capac), 50), np.array([float(capac)]))) for hr in range(len(waterC))]

    if genparam.incRegs and (pt in genparam.ptCurtailedRegs) and (coolType == 'once through'):

        # read set of regression coefficients on water intensity
        regcoeffs = loadRegCoeffs(dataRoot=genparam.dataRoot, fname='water.json')

        # get coefficients for this class of plant
        coeffs = getCoeffsForGenOrTech(pt, coolType, genparam.ptCurtailedRegs, regcoeffs, ppDeltaT)

        # plantFlows is 2d list of np arrays with hourly flows corresponding to each capac level
        plantFlowsAllHrs = calculateGeneratorDischargeFlow(coeffs, metAndWaterData, capacs, streamAvailFrac)

        mixTInputs = [(waterC[hr], waterFlow[hr], int(ppDeltaT), plantFlowsAllHrs[hr]) for hr in range(len(waterC))]

        # mixedTs is 2d list of np arrays w/ mixed Ts corresponding to each capac level
        mixedTsAllHrs = np.array(list(map(calculateMixedTemperature, mixTInputs)))

        idxWaterTExceedMaxT = [(np.where(mixedTsInHr > maxT)[0][0] if np.any(mixedTsInHr >= maxT)
        else -1) for mixedTsInHr in mixedTsAllHrs]

        hrlyCurtsRegs = np.array([capacs[ir][idx] for ir, idx in enumerate(idxWaterTExceedMaxT)]).flatten()

        # divide hourly curtailment by nominal capacity to get available capacity in %
        hrlyCurtsRegs = hrlyCurtsRegs / float(capac)

    return hrlyCurtsRegs


def calculateMixedTemperature(inputs):
    """
    Computes final stream temperature after water discharge by power plant

    :param inputs: list with stream T, stream flow, temperature rise through power plant
    cooling system, power plant discharge flow.
    :return: mixed stream T
    """

    streamT, streamFlow, ppDeltaT, ppFlows = inputs[0], inputs[1], inputs[2], inputs[3]
    # print([type(x) for x in [streamT,streamFlow,ppDeltaT,ppFlows]])
    return streamT + ((ppFlows / streamFlow) / (ppFlows / streamFlow + 1)) * ppDeltaT


def calculateGeneratorDischargeFlow(coeffs, metAndWaterData, capacs, streamAvailFrac):
    """
    Computes Generator discharge flow using Aviva's regression

    :param coeffs: dictionary with regression coefficients for this power plant
    :param metAndWaterData: pandas data frame with meteo and water data for each period
    :param capacs: list of numpy arrays with capacity levels in MW for each period
    :param streamAvailFrac: Max share of stream flow that can be used by power plant (1d numpy array length n_hours)
    :return: 2-d numpy array with water discharge in (gal) for each capacity level (num_per vs num_capacs)
    """

    n_hours = metAndWaterData.shape[0]

    waterFlow = metAndWaterData['flow'].values

    # water flow is in cfs. Convert to gal/h
    waterFlow = waterFlow * (7.48052 * 3600)

    # get water intensity values from Aviva's regression equations
    waterIntensity = setAvivaCurtailments(coeffs, metAndWaterData)

    # multiply each water intensity value by the numpy arrays of capacity (MW) for each hour
    # result is hourly flows corresponding to each capac level (n_hours vs n_capac_levels)
    hrlydischarge = np.array(list(map(lambda x, y: x*y, waterIntensity, capacs)))

    # reshape max discharge for each hour so it is (n_hours vs 1)
    maxdischarge = (streamAvailFrac*waterFlow).reshape(n_hours, 1)

    # compute hourly discharge values
    hrlydischarge = np.minimum(maxdischarge, hrlydischarge)

    return hrlydischarge


def calcCurtailmentForGenOrTech(plantType, fuelAndCoalType, coolType, state, capac, coolDesignT, metAndWaterData,
                                coeffs, genparam, curtailparam):
    """
    CALCULATE CURTAILMENT FOR A SINGLE GENERATOR

    Output 1d list of hourly final curtailments for year for given generator (or generator type) with or without
    curtailments and/or environmental regulations.

    :param plantType: (string) plat type
    :param fuelAndCoalType: (string) fuel type
    :param coolType: (string)  cooling type
    :param state: (string) state where plant is located
    :param capac: (number) capacity of plant
    :param coolDesignT:
    :param metAndWaterData: (pd data frame) data frame with meteo and water data
    :param coeffs: (dict) dictionary with curtailment regression coefficients (plantType,coolType,coolDesignT)
    :param genparam: object of class Generalparameters
    :param curtailparam: object of class Curtailmentparameters
    :return: 1d numpy array
    """

    # add interaction term used by Aviva's regressions
    metAndWaterData['airF:rh'] = metAndWaterData['airF'] * metAndWaterData['rh']

    hourlyCurtailments = runCurtailRegression(metAndWaterData, coeffs, genparam.incCurtailments, plantType,
                                              coolType, genparam.ptCurtailed)

    if state in curtailparam.envRegMaxT.keys():
        # create numpy array with max share of stream flow available in each hour
        streamAvailFrac = np.array([curtailparam.maxFracFlow[str(d.month)] for d in metAndWaterData['date']])

        # get regulatory threshold for state where plant is located
        envRegMaxT = curtailparam.envRegMaxT[state]
        hourlyCurtailmentsRegs = setEnvRegCurtailments(coolType, capac, plantType, coolDesignT, metAndWaterData,
                                                       envRegMaxT, streamAvailFrac, genparam)
    else:
        # if state is not on list. assume no regulatory curtailment
        hourlyCurtailmentsRegs = len(hourlyCurtailments) * [1]

    # get the more restrictive one
    hourlyCurtailments = np.minimum(hourlyCurtailments, hourlyCurtailmentsRegs)

    return hourlyCurtailments


def runCurtailRegression(metAndWaterData, coeffs, incCurtailments, pt, coolType, ptCurtailed):
    """
    Output 1d list w/ hourly curtailment (% of capacity) for generator

    :param metAndWaterData:
    :param coeffs:
    :param incCurtailments:
    :param pt:
    :param coolType:
    :param ptCurtailed:
    :return:
    """

    hrlyCurts = np.ones(metAndWaterData.shape[0])

    if incCurtailments and pt in ptCurtailed:
        hrlyCurts = setAvivaCurtailments(coeffs, metAndWaterData)

    return list(hrlyCurts)


def setAvivaCurtailments(coeffs, metAndWaterData):
    """
    Use Aviva's regressions to curtail generation. Note that curtailment values are fraction of total capacity
    (val of 0.3 means max capac is 30% of total capac)

    :param coeffs: dictionary with regression coefficients for this power plant
    :param metAndWaterData: pandas data frame with meteo and water data
    :return: array with hourly curtailment values (fraction of total capacity)
    """

    hrlyCurtsAviva = np.zeros(metAndWaterData.shape[0])

    # remove values of upper bound and lower bound from dictionary of coefficients
    ubound = coeffs.pop('ubound')
    lbound = coeffs.pop('lbound')

    for param in coeffs:
        if param != 'intercept':
            hrlyCurtsAviva = np.add(hrlyCurtsAviva, coeffs[param] * np.array(metAndWaterData[param]))
        else:
            hrlyCurtsAviva = np.add(hrlyCurtsAviva, coeffs[param])

    # apply lower and upper bound cut-off values
    return np.maximum(np.minimum(hrlyCurtsAviva, ubound), lbound)


def loadRegCoeffs(dataRoot, fname):
    """
    Load regression coefficients from JSON file.

    This function works for both available capacity (%) and water intensity (gal/MWh) files. These files are assumed
    to be inside the folder /dataRoot/Curtailment/.

    :param dataRoot: root of data folder (see pdf 'structure of data folder' for more details on folder structure)
    :param fname: name of JSON file (with extension)
    :return: nested dictionary of cooling type: reg coeffs
    """

    regCoeffs = json.load(fp=open(os.path.join(dataRoot, 'Curtailment', fname)))

    return regCoeffs


def getCoeffsForGenOrTech(plantType, coolType, plantTypesCurtailed, regCoeffs, coolDesignT):
    """
    GET COEFFICIENTS FOR PARTICULAR GEN OR TECH. If plant not eligible for curtailment or cool type not included in
    reg coeffs (e.g., wind), assume no curtailment. See loadRegCoeffs for structure of coeffs.

    Inputs: parameters for getting reg coeffs, and regCoeffs (dict of cool type: cool design T: parameter:coefficient)

    :param plantType: (string) with type of plant ('Coal Steam', 'Hydro', 'Combined Cycle', ...)
    :param coolType: (string) cooling technology type
    :param plantTypesCurtailed: (list) list with plant types that are curtailed
    :param regCoeffs: (dict) nested dictionary with regression coefficient for all (plantType,coolType,coolDesignT)
    :param coolDesignT: (string) cooling technology design temperature
    :return: (dict) dictionary with regression coeffs for this specific combination of (plantType,coolType,coolDesignT)
    """

    # check if this plant has curtailment data
    if (plantType not in plantTypesCurtailed or plantType not in list(regCoeffs.keys()) or
            coolType not in list(regCoeffs[plantType].keys()) or
            coolDesignT not in list(regCoeffs[plantType][coolType].keys())):
        coeffs = None
    else:
        coeffs = copy.deepcopy(regCoeffs[plantType][coolType][coolDesignT])

    # returns deep copy of original dictionary
    return coeffs


def getKeyCurtailParams(gen, genFleet):
    """GET KEY PARAMETERS OF INDIVIDUAL GENERATOR FOR CURTAILMENT OF EXISTING GENERATORS.

    Reads data base of generators and returns Parameters that affect curtailment:
    coal steam vs NGCC, bit vs. subbit vs. lignite, HR, once through vs. recirc vs. dry cooling,
    wet FGD vs. lime spray dryer

    :param gen: (String) id of generator
    :param genFleet: 2d list with generator fleet
    :return: tuple with parameters of generator: (plantType, heat rate, fuel Type, cooling Type,
                                                 fgdType, state, capacity)
    """

    genRow = getGenRowInFleet(gen, genFleet)
    heads = genFleet[0]
    (plantTypeCol, hrCol, coalTypeCol, coolTypeCol, so2ControlCol, stateCol,
     capacCol) = (heads.index('PlantType'), heads.index('Heat Rate (Btu/kWh)'), heads.index('Modeled Fuels'),
                  heads.index('Cooling Tech'), heads.index('Wet/DryScrubber'), heads.index('State Name'),
                  heads.index('Capacity (MW)'))

    coolDesignCol = heads.index('coolingDesignT')


    (plantType, hr, state, capac) = (genRow[plantTypeCol], genRow[hrCol], genRow[stateCol], genRow[capacCol])
    fuelAndCoalType = isolateFirstFuelType(genRow[coalTypeCol])
    coolType = getCoolType(genRow[coolTypeCol])
    fgdType = getSO2Control(genRow[so2ControlCol])
    coolDesignT = genRow[coolDesignCol]

    return plantType, hr, fuelAndCoalType, coolType, fgdType, state, capac, coolDesignT


def getGenRowInFleet(gen, genFleet):
    (orisID, unitID) = separateGenSymbol(gen)
    (orisCol, unitIdCol) = (genFleet[0].index('ORIS Plant Code'), genFleet[0].index('Unit ID'))
    return [row for row in genFleet if str(row[orisCol]) == orisID and row[unitIdCol] == unitID][0]  # return 1d list


def isolateFirstFuelType(fuel):
    multiFuelDivider = '&'  # some plants have multiple modeled fuels divided by &
    if multiFuelDivider in fuel: fuel = fuel[:fuel.index(multiFuelDivider)]
    return fuel


def getCoolType(coolingType):
    possibleCoolingTypes = ['once through', 'dry cooling', 'recirculating']
    finalCoolType = 'NA'
    for possibleCoolingType in possibleCoolingTypes:
        if possibleCoolingType in coolingType: finalCoolType = possibleCoolingType
    return finalCoolType


def getSO2Control(so2Control):
    possibleFGDTypes = ['wet', 'dry']
    fgdType = 'NA'
    for possibleFGDType in possibleFGDTypes:
        if possibleFGDType in so2Control.lower(): fgdType = possibleFGDType
    return fgdType


def getKeyCurtailParamsNewTechs(newTechsCE, techRow):
    """GET KEY PARAMS FOR REGRESSION

    Function returns the key parameters of individual NEW generators needed to do the curtailment simulations

    :param newTechsCE:
    :param techRow:
    :return:
    """

    heads = newTechsCE[0]
    (plantTypeCol, hrCol, fuelTypeCol, coolTypeCol, fgdCol, capCol, coolDesignTCol) = (heads.index('TechnologyType'),
                                                                                       heads.index('HR(Btu/kWh)'),
                                                                                       heads.index('FuelAndCoalType'),
                                                                                       heads.index('Cooling Tech'),
                                                                                       heads.index('SO2 Scrubber'),
                                                                                       heads.index('Capacity(MW)'),
                                                                                       heads.index('coolingDesignT'))

    (plantType, hr, coolType, fgdType, cap, coolingDesignT) = (techRow[plantTypeCol], techRow[hrCol],
                                                               techRow[coolTypeCol], techRow[fgdCol], techRow[capCol],
                                                               techRow[coolDesignTCol])
    fuelAndCoalType = techRow[fuelTypeCol]

    return plantType, hr, fuelAndCoalType, coolType, fgdType, cap, coolingDesignT


def test_results():
    """testing and plotting curtailment results to compare to aviva's plots

    :return:
    """

    (plantType, hr, fuelAndCoalType, coolType, fgdType, state, capac) = getKeyCurtailParams(gen, genFleet)
    coeffs = getCoeffsForGenOrTech(plantType, coolType, genparam.ptCurtailed, regCoeffs,
                                   genparam.coolDesignT)

    x = np.arange(start=70, stop=110.1, step=0.1)
    y = np.arange(start=20, stop=110.1, step=0.1)

    xy = np.array(np.meshgrid(x, y)).reshape(2, x.shape[0]*y.shape[0]).T

    metAndWaterData = pd.DataFrame({'airF': xy[:, 0], 'rh': xy[:, 1], 'airF:rh': xy[:, 0]*xy[:, 1]})

    hourlyCurtailments = runCurtailRegression(metAndWaterData, coeffs, genparam.incCurtailments, plantType,
                                              coolType, genparam.ptCurtailed)

    xx, yy = np.meshgrid(x, y)
    zz = np.array(hourlyCurtailments).reshape(xx.shape)

    fig, ax = plt.subplots(1)

    N = 10
    base = plt.cm.get_cmap(plt.get_cmap('YlGn'))
    color_list = base(np.linspace(0, 1, N))
    cmap_name = base.name + str(N)
    my_map = base.from_list(cmap_name, color_list, N)

    bins = np.concatenate(([-np.inf], np.arange(start=0.1, stop=1, step=0.1), [np.inf]))
    values = np.digitize(zz, bins)

    im = ax.pcolormesh(xx, yy, values, cmap=my_map)

    plt.xlim(70, 110)
    plt.ylim(20, 100)

    cbar = fig.colorbar(im, orientation='vertical')
    cbar.set_ticks([])
    for j in np.arange(N+1, step=2):
        cbar.ax.text(1, j/N, j/N, ha='left', va='center')

    plt.savefig('./example.png')
    plt.close(fig)


