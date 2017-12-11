#Michael Craig
#October 4, 2016
#Functions for adding parameters to GAMS database. Used for CE & UC models.

import copy, math
from GAMSAuxFuncs import *
from CalculateOpCost import calcOpCostsTech,calcOpCosts
from TransmissionLineFuncs import getLineSourceAndSink
from AuxFuncs import convertCostToTgtYr

################################################################################
##################### CE & UC PARAMETERS #######################################
################################################################################
##### ADD HOURLY DEMAND (dict of zone:hourly demand)
def addDemandParam(db,demandCEZonal,zoneSet,hourSet,hourSymbols,ipmZones,ipmZoneNums,scaleMWtoGW): 
    demandCENumZoneSymbs = createDictIndexedByZone(demandCEZonal,ipmZones,ipmZoneNums)
    demandDict = getHourly2dParamDict(demandCENumZoneSymbs,hourSymbols,1/scaleMWtoGW)
    (demandName,demandDescrip) = ('pDemand','hourly zonal demand (GWh)')
    demandParam = add2dParam(db,demandDict,zoneSet,hourSet,demandName,demandDescrip)

#Inputs: dictionary of zone:list, optional scalar. Outputs: dictionary of zone symbol:list.
def createDictIndexedByZone(dataDict,ipmZones,ipmZoneNums,*args):
    zoneDict = dict()
    for zone in dataDict: 
        if len(args)>0: zoneDict[createZoneSymbol(ipmZoneNums[ipmZones.index(zone)])] = dataDict[zone]*args[0]
        else: zoneDict[createZoneSymbol(ipmZoneNums[ipmZones.index(zone)])] = dataDict[zone]
    return zoneDict

##### ADD EXISTING GENERATOR PARAMETERS: CO2 ems rate, zone
def addEguParams(db,genFleet,genSet,genSymbols,ipmZones,ipmZoneNums,scaleLbToShortTon,scaleMWtoGW):
    #Heat rate
    scalarHrToMmbtuPerMwh = 1/1000
    hrDict = getEguParamDict(genFleet,'Heat Rate (Btu/kWh)',scalarHrToMmbtuPerMwh*scaleMWtoGW)
    (hrName,hrDescrip) = ('pHr','heat rate (MMBtu/GWh)')
    hrParam = add1dParam(db,hrDict,genSet,genSymbols,hrName,hrDescrip)
    #Emissions rate
    emRateDict = getEguParamDict(genFleet,'CO2EmRate(lb/MMBtu)',1/scaleLbToShortTon)
    (emRateName,emRateDescrip) = ('pCO2emrate','emissions rate (short ton/MMBtu)')
    emRateParam = add1dParam(db,emRateDict,genSet,genSymbols,emRateName,emRateDescrip)
    #Zone
    zoneDict = getEguParamZoneDict(genFleet,'Region Name',ipmZones,ipmZoneNums)
    (zoneName,zoneDesc) = ('pEguzones','zone for each egu')
    zoneParam = add1dParam(db,zoneDict,genSet,genSymbols,zoneName,zoneDesc)

#Return dict of genSymbol:zone num
def getEguParamZoneDict(genFleet,zoneCol,ipmZones,ipmZoneNums):
    zoneCol = genFleet[0].index(zoneCol)
    zoneDict = dict()
    for row in genFleet[1:]:
        zoneDict[createGenSymbol(row,genFleet[0])] = ipmZoneNums[ipmZones.index(row[zoneCol])]
    return zoneDict    

#Add op cost parameter for existing gens
def addEguOpCostParam(db,genFleet,genSet,genSymbols,scaleLbToShortTon,scaleMWtoGW,scaleDollarsToThousands,*co2Price):
    ocDict = getEguOpCostDict(genFleet,scaleLbToShortTon,scaleMWtoGW,scaleDollarsToThousands,co2Price)
    (ocName,ocDescrip) = ('pOpcost','op cost (thousand$/GWh)')
    ocParam = add1dParam(db,ocDict,genSet,genSymbols,ocName,ocDescrip)
   
##### ADD EXISTING GENERATOR HOURLY CAPACITIES
#Add hourly HR & capac param
def addEguHourlyParams(db,hourlyCapacsCE,genSet,hourSet,hourSymbols,scaleMWtoGW):
    capacDict = getHourly2dParamDict(hourlyCapacsCE,hourSymbols,1/scaleMWtoGW)
    (capacName,capacDescrip) = ('pCapac','hourly capacity (GW)')
    capacParam = add2dParam(db,capacDict,genSet,hourSet,capacName,capacDescrip)   
    # hrDict = getHourly2dParamDict(hourlyHrsCE,hourSymbols,1)
    # (hrName,hrDesc) = ('pHr','hourly hr (mmbtu/gwh)')
    # hrParam = add2dParam(db,hrDict,genSet,hourSet,hrName,hrDesc)   

#Creates dictionary of (key,hourSymbol):param   (for capac or HR, key=gen; for demand
#or wind & solar gen, key=zone)
def getHourly2dParamDict(hourlyParam,hourSymbols,scalar):
    paramDict = dict()
    for key in hourlyParam:
        for idx in range(len(hourSymbols)):
            hourSymbol = hourSymbols[idx]
            paramDict[(key,hourSymbol)] = hourlyParam[key][idx]*scalar
    return paramDict
    
##### ADD EXISTING RENEWABLE COMBINED MAXIMUM GENERATION VALUES
#Converts 1d list of param vals to hour-indexed dicts, then adds dicts to GAMS db
def addExistingRenewableMaxGenParams(db,zoneSet,ipmZones,ipmZoneNums,hourSet,hourSymbols,
                                hourlySolarGenCEZonal,hourlyWindGenCEZonal,scaleMWtoGW):    
    hourlySolarGenCEZonalSymbs = createDictIndexedByZone(hourlySolarGenCEZonal,ipmZones,ipmZoneNums)
    hourlyWindGenCEZonalSymbs = createDictIndexedByZone(hourlyWindGenCEZonal,ipmZones,ipmZoneNums)
    maxSolarGenDict = getHourly2dParamDict(hourlySolarGenCEZonalSymbs,hourSymbols,1/scaleMWtoGW)
    maxWindGenDict = getHourly2dParamDict(hourlyWindGenCEZonalSymbs,hourSymbols,1/scaleMWtoGW) 
    (maxSolarGenName,maxSolarGenDescrip) = ('pMaxgensolar','max combined gen by existing solar')
    solarParam = add2dParam(db,maxSolarGenDict,zoneSet,hourSet,maxSolarGenName,maxSolarGenDescrip)
    (maxWindGenName,maxWindGenDescrip) = ('pMaxgenwind','max combined gen by existing wind')
    windParam = add2dParam(db,maxWindGenDict,zoneSet,hourSet,maxWindGenName,maxWindGenDescrip)

#Stores set of values into dictionary keyed by hour
#Inputs: set of param values (1d list), hour symbols (1d list), optional scalar
#Outputs: dictionary of (hour symbol:param val)
def getParamIndexedByHourDict(paramVals,hourSymbols,*scalar):
    paramIndexedByHourDict = dict()
    for idx in range(len(hourSymbols)): paramIndexedByHourDict[hourSymbols[idx]] = paramVals[idx]*scalar[0]
    return paramIndexedByHourDict

##### ADD ZONE AND LINE CONSTRAINTS
#Add parameter mapping lines to zone sources & sinks (pLinesources(l),pLinesinks(l))
def addLineSourceAndSink(db,lineSet,lines,ipmZones,ipmZoneNums):
    lineSources,lineSinks = dict(),dict()
    for line in lines:
        lineSource,lineSink = getLineSourceAndSink(line)
        lineSources[line] = ipmZoneNums[ipmZones.index(lineSource)] #convert zone from name to 1+ number
        lineSinks[line] = ipmZoneNums[ipmZones.index(lineSink)]
    sourceName,sourceDesc = 'pLinesources','source zone for each line'
    sinkName,sinkDesc = 'pLinesinks','sink zone for each line'
    sourceParam = add1dParam(db,lineSources,lineSet,lines,sourceName,sourceDesc)
    sinkParam = add1dParam(db,lineSinks,lineSet,lines,sinkName,sinkDesc)

#Add line capacities (pLinecapacs(l)). lineCapacs already a dict of line:capac;
#convert to GW
def addLineCapacs(db,lineCapacs,lineSet,lines,scaleMWtoGW):
    lineCapacsGW = dict()
    for key in lineCapacs: lineCapacsGW[key] = lineCapacs[key]*1/scaleMWtoGW
    name,desc = 'pLinecapacs','line capacs (GW)'
    param = add1dParam(db,lineCapacsGW,lineSet,lines,name,desc)

##### ADD PUMPED HYDRO PARAMETERS
def addPumpHydroParams(db,genFleetForCE,phEff,phMaxSoc,phInitSoc,pumpHydroGenSet,pumpHydroGenSymbols,scaleMWtoGW):
    effDict = dict()
    for symb in pumpHydroGenSymbols: effDict[symb] = phEff
    (effname,effdesc) = ('pEfficiency','efficiency')
    effParam = add1dParam(db,effDict,pumpHydroGenSet,pumpHydroGenSymbols,effname,effdesc)
    socDict = getSocDict(genFleetForCE,phMaxSoc,pumpHydroGenSymbols,scaleMWtoGW)
    (socname,socdesc) = ('pMaxsoc','max state of charge (GWh)')
    socparam = add1dParam(db,socDict,pumpHydroGenSet,pumpHydroGenSymbols,socname,socdesc)
    initSocDict = dict()
    for gen in socDict: initSocDict[gen] = socDict[gen]*phInitSoc #given as fraction of max SOC
    (socname,socdesc) = ('pInitsoc','initial state of charge (GWh)')
    initsocparam = add1dParam(db,initSocDict,pumpHydroGenSet,pumpHydroGenSymbols,socname,socdesc)
    
#Get state of charge for pumped hydro. phMaxSoc equals multiple of capacity
def getSocDict(genFleetForCE,phMaxSoc,pumpHydroGenSymbols,scaleMWtoGW):
    capacCol = genFleetForCE[0].index('Capacity (MW)')
    socDict = dict()
    for row in genFleetForCE[1:]:
        if createGenSymbol(row,genFleetForCE[0]) in pumpHydroGenSymbols: 
            socDict[createGenSymbol(row,genFleetForCE[0])] = float(row[capacCol])/scaleMWtoGW*phMaxSoc
    return socDict

################################################################################
################################################################################
################################################################################

################################################################################
##################### CAPACITY EXPANSION PARAMETERS ############################
################################################################################
##### ADD NEW TECH PARAMS FOR CE
def addTechParams(db,newTechsCE,techSet,techSymbols,hourSet,hourSymbols,
                    scaleMWtoGW,scaleDollarsToThousands,scaleLbToShortTon,ptCurtailed):
    #Nameplate capacity (for cost calculations)
    capacDict = getTechParamDict(newTechsCE,techSymbols,'Capacity(MW)',ptCurtailed,
                                                                    1/scaleMWtoGW)
    (capacName,capacDescrip) = ('pCapactech','capacity (GW) of techs')
    techCapacParam = add1dParam(db,capacDict,techSet,techSymbols,capacName,capacDescrip)
    #Heat rate
    scalarHrToMmbtuPerMwh = 1/1000
    hrDict = getTechParamDict(newTechsCE,techSymbols,'HR(Btu/kWh)',ptCurtailed,
                                            scalarHrToMmbtuPerMwh*scaleMWtoGW)
    (hrName,hrDescrip) = ('pHrtech','heat rate (MMBtu/GWh)')
    techHrParam = add1dParam(db,hrDict,techSet,techSymbols,hrName,hrDescrip)
    #Op cost
    ocDict = getTechOpCostDict(newTechsCE,ptCurtailed,scaleMWtoGW/scaleDollarsToThousands)
    (ocName,ocDescrip) = ('pOpcosttech','op cost for tech (thousand$/GWh)')
    ocParam = add1dParam(db,ocDict,techSet,techSymbols,ocName,ocDescrip)
    #Fixed O&M
    fixedomDict = getTechParamDict(newTechsCE,techSymbols,'FOM(2012$/MW/yr)',ptCurtailed,
                                                scaleMWtoGW*1/scaleDollarsToThousands)
    for tech in fixedomDict: fixedomDict[tech] = convertCostToTgtYr('fom',fixedomDict[tech])
    (fixedomName,fixedomDescrip) = ('pFom','fixed O&M (thousand$/GW/yr)')
    techFixedomParam = add1dParam(db,fixedomDict,techSet,techSymbols,fixedomName,fixedomDescrip)
    #Overnight capital cost
    occDict = getTechParamDict(newTechsCE,techSymbols,'CAPEX(2012$/MW)',ptCurtailed,
                                scaleMWtoGW*1/scaleDollarsToThousands)
    for tech in occDict: occDict[tech] = convertCostToTgtYr('occ',occDict[tech])
    (occName,occDescrip) = ('pOcc','overnight capital cost (thousand$/GW)')
    techOccParam = add1dParam(db,occDict,techSet,techSymbols,occName,occDescrip)
    #Emissions rate
    emRateDict = getTechParamDict(newTechsCE,techSymbols,'CO2EmissionsRate(lb/MMBtu)',ptCurtailed,
                                                                            1/scaleLbToShortTon)
    (emRateName,emRateDescrip) = ('pCO2emratetech','co2 emissions rate (short ton/MMBtu)')
    techEmRateParam = add1dParam(db,emRateDict,techSet,techSymbols,emRateName,emRateDescrip)
    #Lifetime
    lifetimeDict = getTechParamDict(newTechsCE,techSymbols,'Lifetime(years)',ptCurtailed)
    (lifetimeName,lifetimeDescrip) = ('pLife','years')
    techLifetimeParam = add1dParam(db,lifetimeDict,techSet,techSymbols,lifetimeName,lifetimeDescrip)

#Creates dict of (techSymbol:paramVal) for given parameter name
def getTechParamDict(newTechsCE,techSymbols,paramColName,ptCurtailed,*scalar):
    techCol = newTechsCE[0].index('TechnologyType')
    paramCol = newTechsCE[0].index(paramColName)
    techSymbolsInNewTechsCE = [createTechSymbol(row,newTechsCE[0],ptCurtailed) for row in newTechsCE]
    paramDict = dict()
    for techSymbol in techSymbols:
        rowIdx = techSymbolsInNewTechsCE.index(techSymbol)
        if len(scalar)>0: paramDict[techSymbol] = float(newTechsCE[rowIdx][paramCol])*scalar[0] 
        else: paramDict[techSymbol] = float(newTechsCE[rowIdx][paramCol])
    return paramDict

#Takes in techs and returns dictionary of (tech:opCost)
def getTechOpCostDict(newTechs,ptCurtailed,scalar):
    opCosts = calcOpCostsTech(newTechs)
    paramDict = dict()
    for idx in range(1,len(newTechs)):
        paramDict[createTechSymbol(newTechs[idx],newTechs[0],ptCurtailed)] = opCosts[idx-1]*scalar #op costs = 1d list of vals, so offset by 1
    return paramDict

#Add hourly curtailed tech capac.
#Input capacs: (tech,loc):[capacs]. Output dict added to db: idxed by (loc,tech,hr)
def addTechCurtailedHourlyCapac(db,hourlyCurtailedTechCapacsCE,cellSet,techCurtailedSet,hourSet,hourSymbols,scaleMWtoGW):
    capacDict = getHourlyTechParamDict(hourlyCurtailedTechCapacsCE,hourSymbols,1/scaleMWtoGW)
    (cName,cDes) = ('pCapactechcurtailed','hourly curtailed capac (GW)')
    capParam = add3dParam(db,capacDict,cellSet,techCurtailedSet,hourSet,cName,cDes)

#Add hourly HR params to tech.
#Input HRs: (tech,loc):[HR]. Output dict added to db: idxed by (loc,tech,hr)
# def addTechHourlyHR(db,hourlyTechHrsCE,cellSet,techSet,hourSet,hourSymbols,scaleMWtoGW):
#     hrDict = getHourlyTechParamDict(hourlyTechHrsCE,hourSymbols,1)
#     (hrName,hrDescrip) = ('pHrtech','heat rate (MMBtu/GWh)')
#     techHrParam = add3dParam(db,hrDict,cellSet,techSet,hourSet,hrName,hrDescrip)
    
#Creates dictionary of (locSymbol,techSymbol,hourSymbol):capac or HR
def getHourlyTechParamDict(hourlyParam,hourSymbols,scalar):
    paramDict = dict()
    for (techSymbol,locSymbol) in hourlyParam:
        for idx in range(len(hourSymbols)):
            hourSymbol = hourSymbols[idx]
            paramDict[(locSymbol,techSymbol,hourSymbol)] = hourlyParam[(techSymbol,locSymbol)][idx]*scalar
    return paramDict

##### ADD MAP FROM CELLS TO ZONES
def addCellsToZones(db,cellSet,cellsToZones,ipmZones,ipmZoneNums):
    cellsToZoneNums,cellSymbols = dict(),list()
    for cell in cellsToZones: 
        cellsToZoneNums[cell] = ipmZoneNums[ipmZones.index(cellsToZones[cell])]
        cellSymbols.append(cell)
    (name,desc) = ('pCellzones','zone each cell is in')
    cToZParam = add1dParam(db,cellsToZoneNums,cellSet,cellSymbols,name,desc)

##### ADD PLANNING RESERVE MARGIN FRACTION PARAMETER
#Add zonal planning reserve
def addPlanningReserveParam(db,planningReserveZonal,ipmZones,ipmZoneNums,zoneSet,zoneSymbols,scaleMWtoGW): 
    reserveDict = createDictIndexedByZone(planningReserveZonal,ipmZones,ipmZoneNums,1/scaleMWtoGW)
    planName,planDesc = 'pPlanningreserve','planning reserve'
    add1dParam(db,reserveDict,zoneSet,zoneSymbols,planName,planDesc)

#Add map of peak hour to zone
def addPeakHourToZoneParam(db,peakDemandHourZonal,peakHourSet,peakHrSymbols,ipmZones,ipmZoneNums):
    peakDict = dict()
    for zone in peakDemandHourZonal: 
        peakDict[createHourSymbol(peakDemandHourZonal[zone])] = ipmZoneNums[ipmZones.index(zone)]
    (peakName,peakDesc) = ('pPeakhtozone','map peak hours to zones')
    zoneParam = add1dParam(db,peakDict,peakHourSet,peakHrSymbols,peakName,peakDesc)

##### ADD DISCOUNT RATE PARAMETER
def addDiscountRateParam(db,discountRate):
    add0dParam(db,'pR','discount rate',discountRate)

##### ADD FIRM FRACTION FOR EXISTING GENERATORS
#Firm fraction goes towards meeting planning reserve margin
def addExistingPlantFirmFractions(db,genFleet,genSet,genSymbols,firmCapacityCreditsExistingGens):
    firmCreditDict = getFirmCreditExistingGenDict(genFleet,firmCapacityCreditsExistingGens)
    (firmCreditName,firmCreditDescrip) = ('pFirmcapacfractionegu','firm capacity fraction')
    firmCreditExistingGenParam = add1dParam(db,firmCreditDict,genSet,genSymbols,firmCreditName,firmCreditDescrip)

#Returns dict of (genSymbol:capacCredit) based on plant type of each generator    
def getFirmCreditExistingGenDict(genFleet,firmCapacityCreditsExistingGens):
    plantTypeCol = genFleet[0].index('PlantType')
    firmCapacityCreditsExistingGensDict = dict()
    for row in genFleet[1:]:
        capacCredit = firmCapacityCreditsExistingGens[row[plantTypeCol]]
        firmCapacityCreditsExistingGensDict[createGenSymbol(row,genFleet[0])] = capacCredit
    return firmCapacityCreditsExistingGensDict

##### ADD HOURLY CAPACITY FACTORS FOR NEW RENEWABLE TECHS
#Add CFs for new renew builds. Input: CFs as zone:[CF]
#Output: dict added to GAMS file as (z,tech,h):[CFs]
def addRenewTechCFParams(db,renewTechSet,renewTechSymbols,zoneSet,hourSet,hourSymbols,
                        newWindCFsCEZonal,newSolarCFsCEZonal,ipmZones,ipmZoneNums):
    renewtechCfDict = dict()
    for renewtech in renewTechSymbols:
        if renewtech == 'Wind': relevantCfs = copy.deepcopy(newWindCFsCEZonal)
        elif renewtech == 'Solar PV': relevantCfs = copy.deepcopy(newSolarCFsCEZonal)
        for zone in relevantCfs:
            for idx in range(len(hourSymbols)): 
                renewtechCfDict[(createZoneSymbol(ipmZoneNums[ipmZones.index(zone)]),
                                renewtech,hourSymbols[idx])] = relevantCfs[zone][idx]
    (renewtechCFName,renewtechCFDescrip) = ('pCf','capacity factors for new wind and solar')
    renewtechCfParam = add3dParam(db,renewtechCfDict,zoneSet,renewTechSet,hourSet,renewtechCFName,renewtechCFDescrip)
    
##### ADD CO2 EMISSIONS CAP
def addCppEmissionsCap(db,co2CppSercCurrYearLimit):
    add0dParam(db,'pCO2emcap','CPP co2 emissions cap [short tons]',co2CppSercCurrYearLimit)

##### ADD WEIGHTS TO SCALE REPRESENTATIVE SEASONAL DEMAND UP
def addSeasonDemandWeights(db,seasonDemandWeights):
    for season in seasonDemandWeights:
        add0dParam(db,'pWeight' + season,'weight on rep. seasonal demand',seasonDemandWeights[season])

##### ADD LIMIT ON MAX NUMBER OF NEW BUILDS PER TECH BY ZONE
def addMaxNumNewBuilds(db,newTechsCE,zoneSet,ipmZones,ipmZoneNums,techSet,maxAddedZonalCapacPerTech,ptCurtailed):
    techMaxNewBuildsDict = getMaxNumBuilds(newTechsCE,ipmZones,ipmZoneNums,maxAddedZonalCapacPerTech,ptCurtailed)
    (maxBuildName,maxBuildDescrip) = ('pNmax','max num builds per tech')
    maxBuildParam = add2dParam(db,techMaxNewBuildsDict,zoneSet,techSet,maxBuildName,maxBuildDescrip)

#Returns dict of (zone,tech):maxCapac
def getMaxNumBuilds(newTechsCE,ipmZones,ipmZoneNums,maxAddedZonalCapacPerTech,ptCurtailed):
    capacCol = newTechsCE[0].index('Capacity(MW)')
    techCol = newTechsCE[0].index('TechnologyType')
    techMaxNewBuildsDict = dict()
    for row in newTechsCE[1:]:
        for zone in ipmZones:
            techMaxNewBuildsDict[(createZoneSymbol(ipmZoneNums[ipmZones.index(zone)]),
                createTechSymbol(row,newTechsCE[0],ptCurtailed))] = math.ceil(maxAddedZonalCapacPerTech/int(row[capacCol]))
    return techMaxNewBuildsDict

###### ADD MAX HYDRO GEN PER TIME BLOCK
#hydroPotPerSeason is a dict of season:genSymbol:gen (MWh), so just need to scale it.
#Note that sesaon can be 'special'!
def addHydroMaxGenPerSeason(db,hydroGenSet,hydroGenSymbols,hydroPotPerSeason,scaleMWtoGW):
    for season in hydroPotPerSeason:
        maxGenDict = dict()
        for symb,pot in hydroPotPerSeason[season].items(): maxGenDict[symb] = pot/scaleMWtoGW
        name,desc = 'pMaxhydrogen' + season[:3],'max gen by each hydro unit' #ex: pMaxhydrogensum
        maxGenParam = add1dParam(db,maxGenDict,hydroGenSet,hydroGenSymbols,name,desc)

################################################################################
################################################################################
################################################################################

################################################################################
##################### UNIT COMMITMENT PARAMETERS ###############################
################################################################################
#Add UC parameters
def addEguUCParams(db,fleetUC,genSet,genSymbols,scaleMWtoGW,scaleDollarsToThousands):
    #Min load
    minLoadDict = getEguParamDict(fleetUC,'MinLoad(MW)',1/scaleMWtoGW)
    (minLoadName,minLoadDescrip) = ('pMinload','min load (GW)')
    minLoadParam = add1dParam(db,minLoadDict,genSet,genSymbols,minLoadName,minLoadDescrip)
    #Ramp rate
    rampDict = getEguParamDict(fleetUC,'RampRate(MW/hr)',1/scaleMWtoGW)
    (rampName,rampDescrip) = ('pRamprate','ramp rate (GW/hr)')
    rampParam = add1dParam(db,rampDict,genSet,genSymbols,rampName,rampDescrip)
    #Start up fixed cost
    startCostDict = getEguParamDict(fleetUC,'StartCost($)',1/scaleDollarsToThousands)
    (startName,startDescrip) = ('pStartupfixedcost','startup fixed cost (thousand$)')
    startCostParam = add1dParam(db,startCostDict,genSet,genSymbols,startName,startDescrip)
    #Min down time
    minDownDict = getEguParamDict(fleetUC,'MinDownTime(hrs)',1)
    (minDownName,minDownDescrip) = ('pMindowntime','min down time (hrs)')
    minDownParam = add1dParam(db,minDownDict,genSet,genSymbols,minDownName,minDownDescrip)

#Add non-time-varying capacity param for generators
# def addEguCapacParam(db,genFleet,genSet,genSymbols,scaleMWtoGW):
#     capacDict = getEguParamDict(genFleet,'Capacity (MW)',1/scaleMWtoGW)
#     (capacName,capacDescrip) = ('pCapac','capacity (GW)')
#     capacParam = add1dParam(db,capacDict,genSet,genSymbols,capacName,capacDescrip)

##### RESERVE PARAMETERS
#Add reg reserve parameters
def addRegReserveParameters(db,regUp,regDown,rrToRegTime,hourSet,hourSymbols,scaleMWtoGW,modelName):
    rampToRegParam = db.add_parameter('pRampratetoregreservescalar',0,'convert ramp rate to reg timeframe')
    rampToRegParam.add_record().value = rrToRegTime
    #Add hourly reg reserves; in CE model, increases w/ built wind, hence diff
    #name than UC model.
    regUpDict = getParamIndexedByHourDict(regUp,hourSymbols,1/scaleMWtoGW)
    if modelName == 'UC': regParamName = 'pRegupreserves'
    elif modelName == 'CE': regParamName = 'pRegupreserveinitial'
    (regUpName,regUpDescr) = (regParamName,'hourly reg up reserves (GWh)')
    regParam = add1dParam(db,regUpDict,hourSet,hourSymbols,regUpName,regUpDescr)
    regDownDict = getParamIndexedByHourDict(regDown,hourSymbols,1/scaleMWtoGW)
    if modelName == 'UC': regParamName = 'pRegdownreserves'
    elif modelName == 'CE': regParamName = 'pRegdownreserveinitial' #not used in current project
    (regDownName,regDownDescr) = (regParamName,'hourly reg down reserves (GWh)')
    regParam = add1dParam(db,regDownDict,hourSet,hourSymbols,regDownName,regDownDescr)

#Add reserve parameter quantities
def addFlexReserveParameters(db,flexRes,rrToFlexTime,hourSet,hourSymbols,scaleMWtoGW,modelName):
    rampToFlexParam = db.add_parameter('pRampratetoflexreservescalar',0,'convert ramp rate to flex timeframe')
    rampToFlexParam.add_record().value = rrToFlexTime
    flexDict = getParamIndexedByHourDict(flexRes,hourSymbols,1/scaleMWtoGW)
    if modelName == 'UC': regParamName = 'pFlexreserves'
    elif modelName == 'CE': regParamName = 'pFlexreserveinitial'
    (flexName,flexDesc) = (regParamName,'hourly flex reserves (GWh)')
    flexParam = add1dParam(db,flexDict,hourSet,hourSymbols,flexName,flexDesc)
    
def addContReserveParameters(db,contRes,rrToContTime,hourSet,hourSymbols,scaleMWtoGW):
    rampToContParam = db.add_parameter('pRampratetocontreservescalar',0,'convert ramp rate to cont timeframe')
    rampToContParam.add_record().value = rrToContTime
    contDict = getParamIndexedByHourDict(contRes,hourSymbols,1/scaleMWtoGW)
    (contName,contDesc) = ('pContreserves','hourly cont reserves (GWh)')
    contParam = add1dParam(db,contDict,hourSet,hourSymbols,contName,contDesc)

#Add initial conditions
def addEguInitialConditions(db,genSet,genSymbols,fleetUC,onOffInitial,genAboveMinInitial,
                            mdtCarriedInitial,scaleMWtoGW):
    onOffInitialDict = getInitialCondsDict(fleetUC,onOffInitial,1)
    (onOffInitialName,onOffInitialDescrip) = ('pOnoroffinitial','whether initially on (1) or off (0) based on last UC')
    onOffInitialParam = add1dParam(db,onOffInitialDict,genSet,genSymbols,onOffInitialName,onOffInitialDescrip)
    mdtCarryDict = getInitialCondsDict(fleetUC,mdtCarriedInitial,1)
    (mdtCarryName,mdtCarryDescrip) = ('pMdtcarriedhours','remaining min down time hrs from last UC (hrs))')
    mdtCarryParam = add1dParam(db,mdtCarryDict,genSet,genSymbols,mdtCarryName,mdtCarryDescrip)
    genAboveMinDict = getInitialCondsDict(fleetUC,genAboveMinInitial,scaleMWtoGW) 
    (genAboveMinName,genAboveMinDescrip) = ('pGenabovemininitial','initial gen above min load based on last UC (GW)')
    genAboveMinParam = add1dParam(db,genAboveMinDict,genSet,genSymbols,genAboveMinName,genAboveMinDescrip)

def getInitialCondsDict(fleetUC,initialCondValues,*scalar):
    initCondsDict = dict()
    for rowNum in range(1,len(fleetUC)):
        initCondsDict[createGenSymbol(fleetUC[rowNum],fleetUC[0])] = initialCondValues[rowNum-1] * scalar[0]
    return initCondsDict

##### WHICH GENERATORS ARE ELIGIBLE TO PROVIDE RESERVES
#Add parameter for which existing generators can provide flex, cont, or reg reserves
def addEguEligibleToProvideRes(db,fleetUC,genSet,genSymbols,*stoMarket):
    fleetOrTechsFlag = 'fleet'
    eligibleFlexDict = getEligibleSpinDict(fleetUC,fleetOrTechsFlag,stoMarket)
    (eligFlexName,eligFlexDesc) = ('pFlexeligible','egu eligible to provide flex (1) or not (0)')
    eligFlexParam = add1dParam(db,eligibleFlexDict,genSet,genSymbols,eligFlexName,eligFlexDesc)
    eligContDict = getEligibleSpinDict(fleetUC,fleetOrTechsFlag,stoMarket)
    (eligContName,eligContDesc) = ('pConteligible','egu eligible to provide cont (1) or not (0)')
    eligContParam = add1dParam(db,eligContDict,genSet,genSymbols,eligContName,eligContDesc)
    eligibleRegDict = getEguParamDict(fleetUC,'RegOfferElig',1)
    (eligRegName,eligRegDescrip) = ('pRegeligible','egu eligible to provide reg res (1) or not (0)')
    eligibleRegParam = add1dParam(db,eligibleRegDict,genSet,genSymbols,eligRegName,eligRegDescrip)

#Returns dict of whether units can provide spin reserves or not based on the plant type
def getEligibleSpinDict(fleetOrTechsData,fleetOrTechsFlag,*stoMktOrPtCurt):
    (windPlantType,solarPlantType) = getWindAndSolarPlantTypes()
    plantTypesNotProvideRes = {windPlantType,solarPlantType}
    if fleetOrTechsFlag == 'fleet':
        if len(stoMktOrPtCurt) > 0:
            if len(stoMktOrPtCurt[0]) > 0 and stoMktOrPtCurt[0][0] == 'energy': 
                plantTypesNotProvideRes.add('Storage') 
        plantTypeCol = fleetOrTechsData[0].index('PlantType')
    elif fleetOrTechsFlag=='techs': 
        if len(stoMktOrPtCurt) > 0 and len(stoMktOrPtCurt[0]) > 0: 
            ptCurtailed = stoMktOrPtCurt[0][0]
            print('check this eligible spin dict in GAMSAddParam!')
        plantTypeCol = fleetOrTechsData[0].index('TechnologyType')
    eligibleSpinDict = dict()
    for rowNum in range(1,len(fleetOrTechsData)):
        plantType = fleetOrTechsData[rowNum][plantTypeCol]
        if plantType in plantTypesNotProvideRes: provideSpin = 0
        else: provideSpin = 1
        if fleetOrTechsFlag=='fleet': symbol = createGenSymbol(fleetOrTechsData[rowNum],fleetOrTechsData[0])
        elif fleetOrTechsFlag=='techs': symbol = createTechSymbol(fleetOrTechsData[rowNum],fleetOrTechsData[0],
                                                                                            ptCurtailed)
        eligibleSpinDict[symbol] = provideSpin
    return eligibleSpinDict

def getWindAndSolarPlantTypes():
    return ('Wind','Solar PV')

#Add cost of CNSE
def addCostNonservedEnergy(db,cnse,scaleMWtoGW,scaleDollarsToThousands):
    add0dParam(db,'pCnse','cost of non-served energy (thousand$/GWh)',
               cnse*scaleMWtoGW*1/scaleDollarsToThousands)

#Add CO2 price
def addCo2Price(db,co2Price,scaleDollarsToThousands):
    add0dParam(db,'pCO2price','co2 emissions price (thousand$/short ton)',
               co2Price*1/scaleDollarsToThousands)
################################################################################
################################################################################
################################################################################

################################################################################
############ GENERIC FUNCTIONS TO ADD PARAMS TO GAMS DB ########################
################################################################################
def add0dParam(db,paramName,paramDescrip,paramValue):
    addedParam = db.add_parameter(paramName,0,paramDescrip)
    addedParam.add_record().value = paramValue

def add1dParam(db,paramDict,idxSet,setSymbols,paramName,paramDescrip):
    addedParam = db.add_parameter_dc(paramName,[idxSet],paramDescrip)
    for idx in setSymbols:
        addedParam.add_record(idx).value = paramDict[idx]
    return addedParam

def add2dParam(db,param2dDict,idxSet1,idxSet2,paramName,paramDescrip):
    addedParam = db.add_parameter_dc(paramName,[idxSet1,idxSet2],paramDescrip)
    for k,v in iter(param2dDict.items()):
        addedParam.add_record(k).value = v
    return addedParam    

def add3dParam(db,param3dDict,idxSet1,idxSet2,idxSet3,paramName,paramDescrip):
    addedParam = db.add_parameter_dc(paramName,[idxSet1,idxSet2,idxSet3],paramDescrip)
    for k,v in iter(param3dDict.items()):
        addedParam.add_record(k).value = v
    return addedParam        

#Takes in gen fleet and param col name, and returns a dictionary of (genSymbol:paramVal) 
#for each row.
def getEguParamDict(genFleet,paramColName,*scalar):
    paramCol = genFleet[0].index(paramColName)
    paramDict = dict()
    for row in genFleet[1:]:
        paramDict[createGenSymbol(row,genFleet[0])] = float(row[paramCol])*scalar[0]
    return paramDict

#Takes in gen fleet and returns dictionary of (genSymbol:opCost)
def getEguOpCostDict(genFleet,scaleLbToShortTon,scaleMWtoGW,scaleDollarsToThousands,*co2Price):
    if len(co2Price[0])>0: (opCosts,hrs) = calcOpCosts(genFleet,scaleLbToShortTon,co2Price[0][0]) #thousand $/GWh
    else: (opCosts,hrs) = calcOpCosts(genFleet,scaleLbToShortTon) #thousand $/GWh
    paramDict = dict()
    for idx in range(1,len(genFleet)):
        genSymb = createGenSymbol(genFleet[idx],genFleet[0])
        paramDict[genSymb] = opCosts[idx-1]*scaleMWtoGW/scaleDollarsToThousands #op costs = 1d list of vals, so offset by 1
    return paramDict
################################################################################
################################################################################
################################################################################

