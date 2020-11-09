#Michael Craig
#October 20, 2016
#Save UC model results for a day into 2d list for all days

# from GAMSAuxFuncs import extract2dVarResultsIntoDictNoLA
from GAMSUtil.GAMSAuxFuncs import createHourSymbol
import pandas as pd

############ SAVE HOURLY X PLANT RESULTS #######################################
#Adds results from curr UC run to input 2d lists
def saveHourlyResultsByPlant(genByPlant,regUpByPlant,regDownByPlant,flexByPlant,contByPlant, turnonByPlant,
                             turnoffByPlant,onOffByPlant,genToRow,hourToColPlant,ucModel,ucDay,daysOpt):
    saveHourByPlantVar(genByPlant,genToRow,hourToColPlant,ucModel,ucDay,daysOpt,'vGen')
    saveHourByPlantVar(regUpByPlant,genToRow,hourToColPlant,ucModel,ucDay,daysOpt,'vRegup')
    # saveHourByPlantVar(regDownByPlant,genToRow,hourToColPlant,ucModel,ucDay,daysOpt,'vRegdown')
    saveHourByPlantVar(flexByPlant,genToRow,hourToColPlant,ucModel,ucDay,daysOpt,'vFlex')
    saveHourByPlantVar(contByPlant,genToRow,hourToColPlant,ucModel,ucDay,daysOpt,'vCont')
    saveHourByPlantVar(turnonByPlant,genToRow,hourToColPlant,ucModel,ucDay,daysOpt,'vTurnon')
    saveHourByPlantVar(turnoffByPlant,genToRow,hourToColPlant,ucModel,ucDay,daysOpt,'vTurnoff')
    saveHourByPlantVar(onOffByPlant,genToRow,hourToColPlant,ucModel,ucDay,daysOpt,'vOnoroff')

#Add entries to 2d list after extracting results only for optimization horizon.
#Inputs: 2d list w/ gen by plants, uc model (gams obj), first day of uc run,
#num days running in optimizaiton horizon (not inc. LA)
def saveHourByPlantVar(varHourByPlantList,genToRow,hourToColPlant,ucModel,ucDay,daysOpt,varName):
    hoursForOptSet = getHoursInOptimHorizon(ucDay,daysOpt)
    for rec in ucModel.out_db[varName]:
        (gen,hour) = (rec.key(0),rec.key(1)) #Vars are indexed as egu,h
        if hour in hoursForOptSet: 
            (rowIdx,colIdx) = (genToRow[gen],hourToColPlant[hour])
            varHourByPlantList[rowIdx][colIdx] = rec.level

#Inputs: curr UC day, # days in optimization horizon. Outputs: set of hours in curr optimization horizon
def getHoursInOptimHorizon(day,daysOpt):
    (firstHour,lastHour) = ((day-1)*24+1,((day-1)+daysOpt)*24) 
    return set([createHourSymbol(hr) for hr in range(firstHour,lastHour+1)])
################################################################################    


def saveHourlyPumpedHydroResults(pumphydroSoc, pumphydroCharge, ucModel, ucDay, daysOpt):
    hoursForOptSet = getHoursInOptimHorizon(ucDay, daysOpt)

    # get ordered hours
    hoursColTot = [h for h in pumphydroSoc[0][1:]]

    #get ordered pumped hydro genIds
    phIDs = [row[0] for row in pumphydroSoc[1:]]

    for rec in ucModel.out_db['vSoc']:
        (gen,hour) = (rec.key(0), rec.key(1)) #Vars are indexed as egu,h
        if hour in hoursForOptSet:
            (rowIdx, colIdx) = (phIDs.index(gen), hoursColTot.index(hour))
            pumphydroSoc[rowIdx+1][colIdx+1] = rec.level

    for rec in ucModel.out_db['vCharge']:
        (gen,hour) = (rec.key(0), rec.key(1)) #Vars are indexed as egu,h
        if hour in hoursForOptSet:
            (rowIdx, colIdx) = (phIDs.index(gen), hoursColTot.index(hour))
            pumphydroCharge[rowIdx+1][colIdx+1] = rec.level


def saveHourlySystemResults(sysResults, ucModel, ucDay, daysOpt):
    """ SAVE HOURLY SYSTEM RESULTS TO PANDAS DATA FRAME

    :param sysResults: pandas data frame with systems results
    :param ucModel: GAMS model object
    :param ucDay: day of UC simulation
    :param daysOpt: number of days in UC simulation
    """

    df = pd.DataFrame(sysResults)

    resultLabelToEqnName = {'nse':'vNse','mcGen':'meetdemand','mcRegup':'meetregupreserves',
                            'mcFlex':'meetflexreserves','mcCont':'meetcontreserves'} #'mcRegdown':'meetregdownreserves'

    hoursForOptSet = getHoursInOptimHorizon(ucDay, daysOpt)

    for result in resultLabelToEqnName:
        varName = resultLabelToEqnName[result]
        for rec in ucModel.out_db[varName]:
            zone = rec.key(0)
            hour = rec.key(1)
            if hour in hoursForOptSet:
                if 'mc' in result:
                    value = rec.marginal
                else:
                    value = rec.level

                # append row with result to data frame
                df = df.append({'zone': zone, 'hour': hour, 'variable': result, 'value': value}, ignore_index=True)

    return df
