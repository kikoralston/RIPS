# Michael Craig
# August 7, 2017
# Script loads regression parameters from Aviva for curtailments,
# and

import copy
import numpy as np
import pandas as pd
import json
import os
from CurtailmentFromEnviroRegs import *
from GAMSAuxFuncs import *
from AuxFuncs import *


def calcCurtailmentForGenOrTech(plantType, hr, fuelAndCoalType, coolType, fgdType, state, capac,
                                ptCurtailed, ptCurtailedRegs, metAndWaterData, incCurtailments, incRegs,
                                envRegMaxT, coolDesignT, coeffs):
    """
    CALCULATE CURTAILMENT FOR NEW AND EXISTING GENS

    Output 1d list of hourly final curtailments for year for given generator (or generator type) with or without
    curtailments and/or environmental regulations.

    :param plantType:
    :param hr:
    :param fuelAndCoalType:
    :param coolType:
    :param fgdType:
    :param state:
    :param capac:
    :param ptCurtailed:
    :param ptCurtailedRegs:
    :param metAndWaterData:
    :param incCurtailments:
    :param incRegs:
    :param envRegMaxT:
    :param coolDesignT:
    :param coeffs:
    :return:
    """

    # hourlyCurtailments = runCurtailRegression(metAndWaterData, coeffs, incCurtailments, incRegs, envRegMaxT[state],
    #  plantType, coolType, ptCurtailed, ptCurtailedRegs)
    hourlyCurtailments = runCurtailRegression(metAndWaterData, coeffs, incCurtailments, plantType, coolType,
                                              ptCurtailed)

    # missing variables: eff, fuelType, ppDeltaT, waterFlow, streamAvailFrac

    # according to the documentation on this method (see pdf)
    eff = 3.412/hr

    # NEED TO CHECK THIS!!!!
    streamAvailFrac = 0.3

    # fuelType = fuelAndCoalType
    # ppDeltaT is the $Delta_g$ in the docuemtn. Where does it come from?
    # waterFlow
    # streamAvailFrac: max fraction of the river flow that can be extracted ($gamma$ in doc).

    hourlyCurtailmentsRegs = setEnvRegCurtailments(incCurtailments, incRegs, ptCurtailedRegs, coolType, capac,
                                                   plantType, eff, fuelAndCoalType, ppDeltaT, metAndWaterData,
                                                   waterFlow, envRegMaxT, streamAvailFrac)
    hourlyCurtailments = np.minimum(hourlyCurtailments, hourlyCurtailmentsRegs)

    return hourlyCurtailments


def runCurtailRegression(metAndWaterData, coeffs, incCurtailments, pt, coolType,
                         ptCurtailed):
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

    :param coeffs:
    :param metAndWaterData:
    :return: array with hourly curtailment values (fraction of total capacity)
    """

    hrlyCurtsAviva = np.zeros(metAndWaterData.shape[0])

    for param in coeffs:
        if param != 'intercept':
            hrlyCurtsAviva = np.add(hrlyCurtsAviva, coeffs[param] * np.array(metAndWaterData[param]))
        else:
            hrlyCurtsAviva = np.add(hrlyCurtsAviva, coeffs[param])

    # regressions can result in val > 1; constrain to 1.
    return np.minimum(hrlyCurtsAviva, 1)


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


def getCoeffsForGenOrTech(plantType, hr, fuelAndCoalType, coolType, fgdType, plantTypesCurtailed,
                          regCoeffs, coolDesignT):
    """
    GET COEFFICIENTS FOR PARTICULAR GEN OR TECH. If plant not eligible for curtailment or cool type not included in
    reg coeffs (e.g., wind), assume no curtailment. See loadRegCoeffs for structure of coeffs.

    Inputs: parameters for getting reg coeffs, and regCoeffs (dict of cool type: cool design T: parameter:coefficient)

    :param plantType: (string) with type of plant ('Coal Steam', 'Hydro', 'Combined Cycle', ...)
    :param hr: (float) heat rate
    :param fuelAndCoalType: (string) fuel type
    :param coolType: (string) cooling technology type
    :param fgdType:
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
        coeffs = regCoeffs[plantType][coolType][coolDesignT]
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
     capacCol) = (heads.index('PlantType'), heads.index('Heat Rate (Btu/kWh)'),
                  heads.index('Modeled Fuels'), heads.index('Cooling Tech'),
                  heads.index('Wet/DryScrubber'), heads.index('State Name'), heads.index('Capacity (MW)'))
    (plantType, hr, state, capac) = (genRow[plantTypeCol], genRow[hrCol], genRow[stateCol], genRow[capacCol])
    fuelAndCoalType = isolateFirstFuelType(genRow[coalTypeCol])
    coolType = getCoolType(genRow[coolTypeCol])
    fgdType = getSO2Control(genRow[so2ControlCol])
    return (plantType, hr, fuelAndCoalType, coolType, fgdType, state, capac)


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

    Function returns the key parameters of individual NEW generators needed to define which regression function to use

    :param newTechsCE:
    :param techRow:
    :return:
    """

    heads = newTechsCE[0]
    (plantTypeCol, hrCol, fuelTypeCol, coolTypeCol, fgdCol) = (heads.index('TechnologyType'),
                                                               heads.index('HR(Btu/kWh)'),
                                                               heads.index('FuelAndCoalType'),
                                                               heads.index('Cooling Tech'), heads.index('SO2 Scrubber'))
    (plantType, hr, coolType, fgdType) = (techRow[plantTypeCol], techRow[hrCol],
                                          techRow[coolTypeCol], techRow[fgdCol])
    fuelAndCoalType = techRow[fuelTypeCol]
    return plantType, hr, fuelAndCoalType, coolType, fgdType