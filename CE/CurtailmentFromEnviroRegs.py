# Michael Craig
# November 11, 2017
# Equations to calculate generator curtailment due to enviro regulations
# on stream temperatures. Past papers have enforced enviro regulations
# at the outlet of the power plant, but regulatiosn actually apply to the
# stream temperature. Thus, my method uses a mass balance equation to
# estimate the stream temperature after mixing with thermal discharge from
# the power plant. Since this temperature depends on output, I determine the
# max power output of a power plant before it trips a regulatory limit.

import pandas as pd
import numpy as np


def setEnvRegCurtailments(coolType, capac, pt, fuelType, ppDeltaT, metAndWaterData, maxT, streamAvailFrac):
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
    :param fuelType: fuel type
    :param ppDeltaT:
    :param metAndWaterData: (pd data frame) water T data and water flow data (pd series of hourly data)
    :param maxT: regulatory stream temperature threshold
    :param streamAvailFrac: Maximum share of stream flow that can be used by power plant
    :return:
    """

    waterC = metAndWaterData['waterC'].values
    waterFlow = metAndWaterData['waterC'].values

    hrlyCurtsRegs = np.ones(metAndWaterData.shape[0])

    capacs = [np.arange(0, float(capac) + .1, 100) for hr in range(len(waterC))]

    if incCurtailments and incRegs and (pt in ptCurtailedRegs) and (coolType == 'once through'):

        # plantFlows is 2d list of np arrays with hourly flows corresponding to each capac level
        plantFlowsAllHrs = calculateGeneratorDischargeFlow(coeffs, metAndWaterData, capacs, streamAvailFrac)

        mixTInputs = [(waterC[hr], waterFlow[hr], ppDeltaT, plantFlowsAllHrs[hr]) for hr in range(len(waterC))]

        # mixedTs is 2d list of np arrays w/ mixed Ts corresponding to each capac level
        mixedTsAllHrs = np.array(list(map(calculateMixedTemperature, mixTInputs)))

        idxWaterTExceedMaxT = [(np.where(mixedTsInHr >= maxT)[0][0] if np.any(mixedTsInHr >= maxT)
        else -1) for mixedTsInHr in mixedTsAllHrs]

        hrlyCurtsRegs = np.array([capacs[0][idx] for idx in idxWaterTExceedMaxT]).flatten()

    return hrlyCurtsRegs


def calculateMixedTemperature(inputs):
    """
    Computes final stream temperature after water discharge by power plant

    Inputs: stream T, stream flow, temperature rise through power plant
    cooling system, power plant discharge flow (array).

    Outputs: mixed stream T

    :param inputs:
    :return:
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
    :param streamAvailFrac: Maximum share of stream flow that can be used by power plant
    :return: 2-d numpy array with water discharge in (gal) for each capacity level (num_per vs num_capacs)
    """

    num_periods = metAndWaterData.shape[0]
    waterIntensity = np.zeros(num_periods)
    waterFlow = metAndWaterData['waterC'].values

    for param in coeffs:
        if param != 'intercept':
            waterIntensity = np.add(hrlyCurtsAviva, coeffs[param] * np.array(metAndWaterData[param]))
        else:
            waterIntensity = np.add(hrlyCurtsAviva, coeffs[param])

    # multiply each water intensity value by the numpy arrays of capacity (MW) for each hour
    # result is hourly flows corresponding to each capac level
    hrlydischarge = np.array(list(map(lambda x, y: x*y, waterIntensity, capacs)))

    maxdischarge = (streamAvailFrac*waterFlow).reshape(num_periods, 1)

    hrlydischarge = np.minimum(maxdischarge, hrlydischarge)

    return hrlydischarge


# def calculateGeneratorDischargeFlow(inputs):
#     """
#     Computes Generator discharge flow
#
#     Inputs: stream T, fraction of stream flow available for cooling, stream flow,
#     gen capacity [MW] (array), gen efficiency, k_os (fraction of waste heat lost to other sinks),
#     temperature rise through generator cooling system.
#     Units: min(m3/s,MW*1/(g/cm^3*J/(g*C)*C)) = min(m^3/s,MW*cm^3/J) = min(m^3/s,1E6J/s*cm^3/J)
#     = min(m^3/s,1E6cm^3/s) = min(m^3/s,m^3/s)
#
#     Outputs: flow of generator discharge
#
#     :param inputs:
#     :return:
#     """
#     streamT, streamFlow, streamAvailFrac, capacs, eff, kos, ppDeltaT = (inputs[0], inputs[1],
#                                                                         inputs[2], inputs[3], inputs[4], inputs[5],
#                                                                         inputs[6])
#     rhow = 1  # g/cm^3 (water density)
#     cp = 4.184  # J/(g*C) (specific heat of water)
#     return [min(streamAvailFrac * streamFlow, capac * ((1 - eff - kos) / eff) * (1 / (rhow * cp * ppDeltaT)))
#             for capac in capacs]
# TEST FUNCTIONS
# def testFuncs():
#     incCurtailments,incRegs,ptCurtailedRegs,coolType,maxT = True,True,['Coal'],'once through',32
#     capac,eff,fuelType,ppDeltaT,pt = 200,.3,'Coal',20,'Coal'
#     metAndWaterData = pd.DataFrame({'test':[1,2,3],'waterC':[10,20,31.9999]})
#     streamAvailFrac = 1
#     waterFlow = [10000,10000,10000]
#     x = setEnvRegCurtailments(incCurtailments,incRegs,ptCurtailedRegs,coolType,capac,pt,
#         eff,fuelType,ppDeltaT,metAndWaterData,waterFlow,maxT,streamAvailFrac)
#     print(x)
#     waterFlow = [1,2,10000]
#     x = setEnvRegCurtailments(incCurtailments,incRegs,ptCurtailedRegs,coolType,capac,pt,
#         eff,fuelType,ppDeltaT,metAndWaterData,waterFlow,maxT,streamAvailFrac)
#     print(x)
#     waterFlow = [10000,10000,10000]
#     metAndWaterData = pd.DataFrame({'test':[1,2,3],'waterC':[50,50,50]})
#     x = setEnvRegCurtailments(incCurtailments,incRegs,ptCurtailedRegs,coolType,capac,pt,
#         eff,fuelType,ppDeltaT,metAndWaterData,waterFlow,maxT,streamAvailFrac)
#     print(x)
#     metAndWaterData = pd.DataFrame({'test':[1,2,3],'waterC':[10,20,30]})
#     x = setEnvRegCurtailments(incCurtailments,incRegs,ptCurtailedRegs,coolType,capac,pt,
#         eff,fuelType,ppDeltaT,metAndWaterData,waterFlow,maxT,streamAvailFrac)
#     print(x)

# testFuncs()
