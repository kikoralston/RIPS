#Michael Craig
#August 7, 2017
#Script loads regression parameters from Aviva for curtailments, 
#and 

import copy
import numpy as np
import pandas as pd
from CurtailmentFromEnviroRegs import runEnvRegCurtail
from GAMSAuxFuncs import *
from AuxFuncs import *

###### CALCULATE CURTAILMENT FOR NEW AND EXISTING GENS #########################
#Output 1d list of hourly curtailments for year for given generator with or without
#curtialments and/or environmental regulations.
def calcCurtailmentForGenOrTech(plantType,hr,fuelAndCoalType,coolType,fgdType,state,capac,
            ptCurtailed,ptCurtailedRegs,metAndWaterData,incCurtailments,incRegs,
            envRegMaxT,coolDesignT,coeffs):    
    hourlyCurtailments = runCurtailRegression(metAndWaterData,coeffs,incCurtailments,
                                    incRegs,envRegMaxT[state],plantType,coolType,ptCurtailed,ptCurtailedRegs)
    hourlyCurtailmentsRegs = setEnvRegCurtailments(incCurtailments,incRegs,ptCurtailedRegs,coolType,capac,
                pt,eff,fuelType,ppDeltaT,metAndWaterData,waterFlow,envRegMaxT,streamAvailFrac)
    hourlyCurtailments = np.minimum(hourlyCurtailments,hourlyCurtailmentsRegs)
    return hourlyCurtailments
################################################################################

###### CALCULATE CURTAILMENTS FROM REGRESSIONS FOR 1 GENERATOR #################
#Output 1d list w/ hourly curtailment (MW) for generator
def runCurtailRegression(metAndWaterData,coeffs,incCurtailments,pt,coolType,
            ptCurtailed):
    hrlyCurts = np.array([1 for i in range(metAndWaterData.shape[0])])
    if incCurtailments == True and pt in ptCurtailed: hrlyCurts = setAvivaCurtailments(coeffs,metAndWaterData)
    return list(hrlyCurts)

#Use Aviva's regressions to curtail generation. Note that curtailment values
#are fraction of total capacity (val of 0.3 means max capac is 30% of total capac)
def setAvivaCurtailments(coeffs,metAndWaterData):
    hrlyCurtsAviva = np.array([0 for i in range(metAndWaterData.shape[0])])
    for param in coeffs:
        if param != 'intercept': 
            hrlyCurtsAviva = np.add(hrlyCurtsAviva,coeffs[param]*np.array(metAndWaterData[param]))
        else:
            hrlyCurtsAviva = np.add(hrlyCurtsAviva,coeffs[param])
    hrlyCurtsAviva[hrlyCurtsAviva>1] = 1 #regressions can result in val >1; constrain to 1.
    return hrlyCurtsAviva
################################################################################

###### LOAD REGRESSION COEFFICIENTS ############################################
#dict of cooling type: reg coeffs
def loadRegCoeffs(runLoc): 
    regCoeffs,coalCoeffs,ngccCoeffs = dict(),dict(),dict()
    coalCoeffs['once through'] = {90:{'waterF':-.024,'intercept':2.22},
                                100:{'waterF':-.016,'intercept':1.97}} 
    coalCoeffs['recirculating'] = {90:{'airF':-.015,'rh':-.0058,'intercept':2.54},
                                 100:{'airF':-.002,'rh':-.0006,'intercept':1.15}}
    coalCoeffs['dry cooling'] = {90:{'airF':-.00177,'intercept':.36},
                                100:{'airF':-.00177,'intercept':.36}}
    ngccCoeffs['once through'] = {100:{'airF':-.00132,'intercept':.06}} 
    ngccCoeffs['recirculating'] = {100:{'airF':-.0013,'intercept':.0628}}
    ngccCoeffs['dry cooling'] = {100:{'airF':-.0017,'rh':5.38E-6,'intercept':.1}}
    regCoeffs['Coal Steam'] = coalCoeffs
    regCoeffs['Combined Cycle'] = ngccCoeffs
    return regCoeffs
################################################################################

###### GET COEFFICIENTS FOR PARTICULAR GEN OR TECH #############################
#If plant not eligible for curtailment or cool type not included in reg coeffs (e.g., wind), 
#assume no curtailment. See loadRegCoeffs for strucutre of coeffs.
#Inputs: parameters for getting reg coeffs, and regCoeffs (dict of cool type: cool design T:
#parameter:coefficient)
def getCoeffsForGenOrTech(plantType,hr,fuelAndCoalType,coolType,fgdType,plantTypesCurtailed,
                    regCoeffs,coolDesignT):
    if (plantType not in plantTypesCurtailed or coolType not in list(regCoeffs[plantType].keys())
            or coolDesignT not in list(regCoeffs[plantType][coolType].keys())): 
        coeffs = None
    else:
        coeffs = regCoeffs[plantType][coolType][coolDesignT]
    return coeffs
################################################################################

###### GET KEY PARAMETERS FOR CURTAILMENT OF EXISTING GENERATORS ###############
#Parameters that affect curtailment: coal steam vs NGCC, bit vs. subbit vs. lignite,
#HR, once through vs. recirc vs. dry cooling, wet FGD vs. lime spray dryer
def getKeyCurtailParams(gen,genFleet):
    genRow = getGenRowInFleet(gen,genFleet)
    heads = genFleet[0]
    (plantTypeCol,hrCol,coalTypeCol,coolTypeCol,so2ControlCol,stateCol,
        capacCol) = (heads.index('PlantType'),heads.index('Heat Rate (Btu/kWh)'),
        heads.index('Modeled Fuels'),heads.index('Cooling Tech'),
        heads.index('Wet/DryScrubber'),heads.index('State Name'),heads.index('Capacity (MW)'))
    (plantType,hr,state,capac) = (genRow[plantTypeCol],genRow[hrCol],genRow[stateCol],genRow[capacCol])
    fuelAndCoalType = isolateFirstFuelType(genRow[coalTypeCol])
    coolType = getCoolType(genRow[coolTypeCol])
    fgdType = getSO2Control(genRow[so2ControlCol])
    return (plantType,hr,fuelAndCoalType,coolType,fgdType,state,capac)

def getGenRowInFleet(gen,genFleet):
    (orisID,unitID) = separateGenSymbol(gen)
    (orisCol,unitIdCol) =  (genFleet[0].index('ORIS Plant Code'),genFleet[0].index('Unit ID')) 
    return [row for row in genFleet if str(row[orisCol]) == orisID and row[unitIdCol] == unitID][0] #return 1d list

def isolateFirstFuelType(fuel):
    multiFuelDivider = '&' #some plants have multiple modeled fuels divided by &
    if multiFuelDivider in fuel: fuel = fuel[:fuel.index(multiFuelDivider)]
    return fuel

def getCoolType(coolingType):
    possibleCoolingTypes = ['once through','dry cooling','recirculating']
    finalCoolType = 'NA'
    for possibleCoolingType in possibleCoolingTypes:
        if possibleCoolingType in coolingType: finalCoolType = possibleCoolingType
    return finalCoolType
    
def getSO2Control(so2Control):
    possibleFGDTypes = ['wet','dry']
    fgdType = 'NA'
    for possibleFGDType in possibleFGDTypes:
        if possibleFGDType in so2Control.lower(): fgdType = possibleFGDType
    return fgdType
################################################################################

####### GET KEY PARAMS FOR REGRESSION ##########################################
def getKeyCurtailParamsNewTechs(newTechsCE,techRow):
    heads = newTechsCE[0]
    (plantTypeCol,hrCol,fuelTypeCol,coolTypeCol,fgdCol) = (heads.index('TechnologyType'),
                    heads.index('HR(Btu/kWh)'),heads.index('FuelAndCoalType'),
                    heads.index('Cooling Tech'),heads.index('SO2 Scrubber'))
    (plantType,hr,coolType,fgdType) = (techRow[plantTypeCol],techRow[hrCol],
                                        techRow[coolTypeCol],techRow[fgdCol])
    fuelAndCoalType = techRow[fuelTypeCol]
    return (plantType,hr,fuelAndCoalType,coolType,fgdType)
################################################################################
