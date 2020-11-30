"""
Michael Craig
October 4, 2016
Converts annual CO2 cap into a CO2 price. Two methods: run ED model in GAMS
or use generation stack.
"""

from operator import *
from DemandFuncs import *
from renewables.GetHydroMaxGenPotential import importHydroPotentialGen,getMonthsPotential
try:
    from gams import *
except ImportError:
    print('gams module not found! GAMS functions will not work.')

from GAMSUtil.GAMSAddParamToDatabaseFuncs import *
from GAMSUtil.GAMSAddSetToDatabaseFuncs import *
from GAMSUtil.GAMSAuxFuncs import *
from TrimDemandREGenAndResForUC import getDemandAndREGenForUC


def convertCo2CapToPrice(genFleet,hourlyWindGen,hourlySolarGen,demandScaled,co2Cap,
        scaleMWtoGW,scaleDollarsToThousands,scaleLbToShortTon,dataRoot,currYear,useGAMS=False):
    """MAIN FUNCTION

    Inputs: gen fleet (2d list), 1d lists of hourly wind & solar gen & total demand,
    annual co2 cap (short tons), scalars, default flag indicating whether use GAMS
    or gen stack approach.
    Outputs: co2 price ($/ton)


    :param genFleet:
    :param hourlyWindGen:
    :param hourlySolarGen:
    :param demandScaled:
    :param co2Cap:
    :param scaleMWtoGW:
    :param scaleDollarsToThousands:
    :param scaleLbToShortTon:
    :param dataRoot:
    :param currYear:
    :param useGAMS:
    :return:
    """
    if useGAMS == False:
        return convertCo2CapToPriceWithGenStack(genFleet,hourlyWindGen,hourlySolarGen,
                                                demandScaled,co2Cap,scaleLbToShortTon,currYear)
    else:
        return convertCo2CapToPriceWithGAMS(genFleet,hourlyWindGen,hourlySolarGen,demandScaled,co2Cap,
                        scaleMWtoGW,scaleDollarsToThousands,dataRoot)


def convertCo2CapToPriceWithGenStack(genFleet,hourlyWindGen,hourlySolarGen,
                                demandScaled,co2Cap,scaleLbToShortTon,currYear):
    """CALCULATE CO2 PRICE TO COMPLY WITH CO2 LIMIT USING GEN STACK

    Calculate net demand, then for each hour of year, calculate op cost of all gens
    for given co2 price, arrange in increasing cost order, and meet demand using
    gen stack. Based on gen stack, calculate hourly emissions, repeat for whole year,
    and compare to annual CO2 cap.
    Inputs: gen fleet (2d list), hourly wind & solar gen & demand (1d lists),
    annual co2 cap (tons).
    Outputs: co2 price ($/ton)

    :param genFleet:
    :param hourlyWindGen:
    :param hourlySolarGen:
    :param demandScaled:
    :param co2Cap:
    :param scaleLbToShortTon:
    :param currYear:
    :return:
    """
    netDemand = calcNetDemand(demandScaled,hourlyWindGen,hourlySolarGen)
    reFuelTypes = ['Wind','Solar']
    fuelCol = genFleet[0].index('Modeled Fuels')
    genFleetNoRE = [row for row in genFleet if row[fuelCol] not in reFuelTypes]
    netDemandNoHydro = removeHydroFromDemand(genFleetNoRE,currYear,netDemand)
    genFleetNoREOrH = [row for row in genFleetNoRE if ((row[fuelCol] != 'Hydro') 
                                                        or (row[fuelCol] != 'Pumped Storage'))]
    print('fleet:',genFleetNoREOrH)
    (genBaseOpCosts,genHrs) = calcOpCosts(genFleetNoREOrH,scaleLbToShortTon)
    (genCapacities,genCo2Ems) = getGenCapacsAndCo2Ems(genFleetNoREOrH,genHrs,scaleLbToShortTon) #MW, ton/MWh
    (annualCo2Ems,co2Price) = (co2Cap*10,-1) #ton, $/ton
    while annualCo2Ems > co2Cap:
        co2Price += 1
        annualCo2Ems = calculateAnnualCo2Ems(genCo2Ems,genCapacities,netDemand,
                                        co2Price,genFleetNoREOrH,scaleLbToShortTon)
        print(annualCo2Ems)
    return co2Price


def removeHydroFromDemand(genFleetNoRE, currYear, netDemand):
    """Net out hydro generation from demand by loading monthly potentials of hydro units,combining potentials, and
    then assumign load follow.

    :param genFleetNoRE:
    :param currYear:
    :param netDemand:
    :return:
    """
    daysPerMonth = {1:31,2:28,3:31,4:30,5:31,6:30,7:31,8:31,9:30,10:31,11:30,12:31}
    plantCol,orisCol = genFleetNoRE[0].index('PlantType'),genFleetNoRE[0].index('ORIS Plant Code')
    capacCol = genFleetNoRE[0].index('Capacity (MW)')

    # 2d list of monthly gen by ORIS in current year; col 1 = ID
    hydroPots = importHydroPotentialGen(currYear,dataRoot)

    hydroUnits = [row for row in fleet[1:] if row[plantCol] == 'Hydro']
    monthlyPotential,netDemandNetHydro,lastHr = [0 for idx in range(12)],list(),0
    for month in range(1,13):
        capacNoData,capacData = 0,0
        for row in hydroUnits:
            oris, capac = row[orisCol],row[capacCol]
            potential,hasData = getMonthsPotential(oris,hydroPots,month)
            print('**', oris, capac, potential)
            if hasData == False: capacNoData += capac
            else: 
                monthlyPotential[month-1] += potential
                capacData += capac
            print('monthly pot:',monthlyPotential)
        monthlyPotential[month-1] += monthlyPotential[month-1]*capacNoData/capacData #fill in missing data
        print('monthly potential w/ filled missing data:',monthlyPotential,capacNoData,capacData)
        monthlyNetDemand,lastHr = netDemand[lastHr:daysPerMonth[month]*24-1],lastHr + daysPerMonth[month]*24
        print('monthly net demand:',monthlyNetDemand)
        shareMoDemand = [val/sum(monthlyNetDemand) for val in monthlyNetDemand]
        print(shareMoDemand,sum(shareMoDemand),'*HEY')
        hourlyHydroGen = [val * monthlyPotential for val in shareMoDemand]
        print(sum(hourlyHydroGen),',and,',hourlyHydroGen)
        netDemandNetHydro.extend([monthlyNetDemand[idx]-hourlyHydroGen[idx] for idx in range(len(monthlyNetDemand))])
        print(netDemandNetHydro)
    return netDemandNetHydro


def getGenCapacsAndCo2Ems(genFleetNoRE,genHrs,scaleLbToShortTon):
    """Returns 1d list of gen fleet's emissions (ton/MWh) & capacities (MW)

    :param genFleetNoRE:
    :param genHrs:
    :param scaleLbToShortTon:
    :return:
    """
    capacCol = genFleetNoRE[0].index('Capacity (MW)')
    co2Col = genFleetNoRE[0].index('CO2EmRate(lb/MMBtu)')
    genCapacities = [float(row[capacCol]) for row in genFleetNoRE[1:]]
    genCo2Ems = [float(row[co2Col])/scaleLbToShortTon for row in genFleetNoRE[1:]] #ton/MMBtu
    genCo2Ems = list(map(mul,genCo2Ems,genHrs)) #ton/MWh
    return (genCapacities, genCo2Ems)


def calculateAnnualCo2Ems(genCo2Ems,genCapacities,netDemand,co2Price,genFleetNoRE,scaleLbToShortTon):
    """

    Calculates annual CO2 ems using gen stack.

    Inputs: 1d lists of op costs (HR*FC+VOM) ($/MWh), CO2 ems (ton/MWh), capacities (MW),
    and net demand (MW), plus a CO2 price ($/ton)
    Outputs: annual co2 ems (tons)

    :param genCo2Ems:
    :param genCapacities:
    :param netDemand:
    :param co2Price:
    :param genFleetNoRE:
    :param scaleLbToShortTon:
    :return:
    """
    (genNewOpCosts,genHrs) = calcOpCosts(genFleetNoRE,scaleLbToShortTon,co2Price)
    costCapacEms = [[genNewOpCosts[idx],genCapacities[idx],genCo2Ems[idx]] 
                            for idx in range(len(genNewOpCosts))]
    sortedCostCapacEms = sorted(costCapacEms)
    sortedCapacs = [row[1] for row in sortedCostCapacEms]
    cumCapacs = [sum(sortedCapacs[:idx+1]) for idx in range(len(sortedCapacs))]
    hourlyCo2Ems = [getHourCo2Ems(hourNetDemand, cumCapacs, sortedCostCapacEms) for hourNetDemand in netDemand]
    return sum(hourlyCo2Ems)


def getHourCo2Ems(hourNetDemand,cumCapacs,sortedCostCapacEms):
    """For given hourly net demand value, calculate total system CO2 ems using gen stack.

    Inputs: net demand for 1 hour, cumulative capacity of stack in increasing cost
    order including co2 price (1d list), 2d list of cost & capacity & co2 ems rate
    in increasing cost order including co2 price
    Outputs: total system co2 ems in 1 hour

    :param hourNetDemand:
    :param cumCapacs:
    :param sortedCostCapacEms:
    :return:
    """
    supplyGap = [hourNetDemand - cumCapac  for cumCapac in cumCapacs]
    supplyGapLessThanZero = [le(gap,0) for gap in supplyGap] #le is <=
    idxMeetDemand = supplyGapLessThanZero.index(True)
    hourCo2Ems = sum([sortedCostCapacEms[idx][2]*sortedCostCapacEms[idx][1] 
                        for idx in range(idxMeetDemand)]) #excludes last unit
    hourCo2Ems += (sortedCostCapacEms[idxMeetDemand][2] *
                   (sortedCostCapacEms[idxMeetDemand][1]+supplyGap[idxMeetDemand])) #include last unit
    return hourCo2Ems


def convertCo2CapToPriceWithGAMS(fleetUC, hourlyWindGen, hourlySolarGen, demandScaled, co2Cap, scaleMWtoGW,
                                 scaleDollarsToThousands, dataRoot):
    """CALCULATE CO2 PRICE TO COMPLY WITH CO2 LIMIT USING GAMS

    Converts CO2 mass limit to CO2 price [$/ton]

    :param fleetUC:
    :param hourlyWindGen:
    :param hourlySolarGen:
    :param demandScaled:
    :param co2Cap:
    :param scaleMWtoGW:
    :param scaleDollarsToThousands:
    :param dataRoot:
    :return:
    """
    co2Emissions = co2Cap + 1
    (co2Price, co2PriceIncrement) = (0,1)

    while co2Emissions>co2Cap:
        co2Emissions = runEdAndGetCo2Ems(fleetUC, hourlyWindGen, hourlySolarGen, demandScaled, co2Price, scaleMWtoGW,
                                         scaleDollarsToThousands, dataRoot)

        print('CO2 mass limit:', co2Cap,', CO2 ems:', co2Emissions,', CO2 price:', co2Price)

        if co2Emissions>co2Cap:
            co2Price += co2PriceIncrement

        elif co2Emissions < co2Cap and co2Price == 0:
            return co2Price

        elif co2Emissions < co2Cap:
            while co2Emissions<co2Cap:
                co2Price-=1
                co2Emissions = runEdAndGetCo2Ems(fleetUC, hourlyWindGen, hourlySolarGen, demandScaled, co2Price,
                                                 scaleMWtoGW, scaleDollarsToThousands)
                print(co2Emissions)

                if co2Emissions > co2Cap:
                    print('Final CO2 price:', co2Price, '$/ton')

                    return co2Price+1


def runEdAndGetCo2Ems(fleetUC,hourlyWindGen,hourlySolarGen,demandScaled,co2Price,
                    scaleMWtoGW,scaleDollarsToThousands,dataRoot):
    fuelCol = fleetUC[0].index('Modeled Fuels')
    fleetED = [row for row in fleetUC if (row[fuelCol] != 'Wind' and row[fuelCol] != 'Solar')]
    (yearCO2Ems,dayIncrement,daysInYear) = (0,5,365)
    for day in range(1,daysInYear+1,dayIncrement):
        (demandED,hourlyWindGenED,hourlySolarGenED,hoursForED) = getDemandAndREGenForUC(day,
                                                                dayIncrement,demandScaled,
                                                                hourlyWindGen,hourlySolarGen)
        netDemandED = calcNetDemand(demandED,hourlyWindGenED,hourlySolarGenED)
        edModel = callEconDispatch(fleetED,netDemandED,hoursForED,co2Price,scaleMWtoGW,scaleDollarsToThousands,dataRoot)
        yearCO2Ems += extractDailyCO2Ems(edModel)
    return yearCO2Ems


def callEconDispatch(fleetED, netDemandED, hoursForED, co2Price, scaleMWtoGW, scaleDollarsToThousands, dataRoot):

    currDir = os.getcwd()

    gamsFileDir = os.path.join(dataRoot, 'GAMS')
    gamsSysDir = 'C:\\GAMS\\win64\\24.7'

    wsED = GamsWorkspace(working_directory=gamsFileDir, system_directory=gamsSysDir)
    dbED = wsED.add_database()

    (genSet,genSymbols,hourSet,hourSymbols) = addSetsEconDispatch(dbED, fleetED, hoursForED)
    addParametersToEconDispatch(dbED,hourSet, hourSymbols, genSet, genSymbols, netDemandED,
                                fleetED,co2Price,scaleMWtoGW,scaleDollarsToThousands)

    edFile = 'EconDispatch25June2016.gms'
    edModel = wsED.add_job_from_file(edFile)
    optED = GamsOptions(wsED)
    optED.defines['gdxincname'] = dbED.name
    edModel.run(optED, databases=dbED)
    return edModel


def addSetsEconDispatch(dbED,fleetED,hoursForED):
    (hourSet,hourSymbols) = addHourSet(dbED,hoursForED)
    genSymbols = isolateGenSymbols(fleetED,'')
    (genSetName,genSetDescription,genSetDimension) = ('egu','existing generators',1)
    genSet = addSet(dbED,genSymbols,genSetName,genSetDescription,genSetDimension) 
    return (genSet,genSymbols,hourSet,hourSymbols)


def addParametersToEconDispatch(dbED,hourSet,hourSymbols,genSet,genSymbols,netDemandED,
                                fleetED,co2Price,scaleMWtoGW,scaleDollarsToThousands):
    addDemandParam(dbED,netDemandED,hourSet,hourSymbols,scaleMWtoGW) 
    addEguParams(dbED,fleetED,genSet,genSymbols,scaleMWtoGW,scaleDollarsToThousands)
    addEguCapacParam(dbED,fleetED,genSet,genSymbols,scaleMWtoGW)
    addCo2Price(dbED,co2Price,scaleDollarsToThousands) 


def extractDailyCO2Ems(edModel):
    return extract0dVarResultsFromGAMSModel(edModel,'vCO2ems')