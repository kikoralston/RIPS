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


# Calculate environmental regulatory curtailments for given time series
# of water temperatures and flows for a single generator. Curtailment
# values are multiplied by capacity to get available capacity, so value output
# by this function of 0 means full curtailment.
# Inputs: met and water data (PD Dataframe, hourly timescale), whether to
# model curtailments and curtailments from regulations, plant types that are
# curtailed, plant cooling type, plant capcity, plant plant type, plant efficiency,
# plant fuel type, plant cooling system designed change in T, water T data (
# pd series of hourly data), water flow data (pd series of hourly data), regulatory
# limit on stream T, and fraction of stream water available for cooling.
def setEnvRegCurtailments(incCurtailments, incRegs, ptCurtailedRegs, coolType, capac,
                          pt, eff, fuelType, ppDeltaT, metAndWaterData, waterFlow, maxT, streamAvailFrac):
    # Set k_os based on Bartos & Chester values of .12 for coal & .2 for gas;
    # assume non-gas plants all have .12 value
    kos = .2 if fuelType == 'Natural Gas' else .12
    waterC = metAndWaterData['waterC'].values
    hrlyCurtsRegs = np.array([1 for i in range(metAndWaterData.shape[0])])
    capacs = [np.arange(0, capac + .1, 100) for hr in range(len(waterC))]
    if (incCurtailments == True and incRegs == True and pt in ptCurtailedRegs
            and coolType == 'once through'):
        dischargeInputs = [(waterC[hr], waterFlow[hr], streamAvailFrac, capacs[hr], eff,
                            kos, ppDeltaT) for hr in range(len(waterC))]
        # plantFlows is 2d list of np arrays with hourly flows corresponding to each capac level
        plantFlowsAllHrs = np.array(list(map(calculateGeneratorDischargeFlow, dischargeInputs)))
        mixTInputs = [(waterC[hr], waterFlow[hr], ppDeltaT, plantFlowsAllHrs[hr]) for hr in range(len(waterC))]
        # mixedTs is 2d list of np arrays w/ mixed Ts corresponding to each capac level
        mixedTsAllHrs = np.array(list(map(calculateMixedTemperature, mixTInputs)))
        idxWaterTExceedMaxT = [(np.where(mixedTsInHr >= maxT)[0][0] if np.any(mixedTsInHr >= maxT)
        else -1) for mixedTsInHr in mixedTsAllHrs]
        hrlyCurtsRegs = np.array([capacs[0][idx] for idx in idxWaterTExceedMaxT]).flatten()
    return hrlyCurtsRegs


# Inputs: stream T, stream flow, temperature rise through power plant
# cooling system, power plant discharge flow (array).
# Outputs: mixed stream T
def calculateMixedTemperature(inputs):
    streamT, streamFlow, ppDeltaT, ppFlows = inputs[0], inputs[1], inputs[2], inputs[3]
    # print([type(x) for x in [streamT,streamFlow,ppDeltaT,ppFlows]])
    return streamT + ((ppFlows / streamFlow) / (ppFlows / streamFlow + 1)) * ppDeltaT


# Inputs: stream T, fraction of stream flow available for cooling, stream flow,
# gen capacity [MW] (array), gen efficiency, k_os (fraction of waste heat lost to other sinks),
# temperature rise through generator cooling system.
# Units: min(m3/s,MW*1/(g/cm^3*J/(g*C)*C)) = min(m^3/s,MW*cm^3/J) = min(m^3/s,1E6J/s*cm^3/J)
# = min(m^3/s,1E6cm^3/s) = min(m^3/s,m^3/s)
# Outputs: flow of generator discharge
def calculateGeneratorDischargeFlow(inputs):
    streamT, streamFlow, streamAvailFrac, capacs, eff, kos, ppDeltaT = (inputs[0], inputs[1],
                                                                        inputs[2], inputs[3], inputs[4], inputs[5],
                                                                        inputs[6])
    rhow = 1  # g/cm^3 (water density)
    cp = 4.184  # J/(g*C) (specific heat of water)
    return [min(streamAvailFrac * streamFlow, capac * ((1 - eff - kos) / eff) * (1 / (rhow * cp * ppDeltaT)))
            for capac in capacs]

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
