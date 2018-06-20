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
from CurtailmentRegressions import loadRegCoeffs, getCoeffsForGenOrTech

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
