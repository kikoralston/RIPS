#Michael Craig
#October 4, 2016
#Helper functions for Python-GAMS interaction

#CREATE SYMBOLS
def createGenSymbol(row,headers):
    (orisCol,unitCol) = (headers.index('ORIS Plant Code'),headers.index('Unit ID'))
    return str(row[orisCol]) + '+' + str(row[unitCol])

def createTechSymbol(row,headers,ptCurtailed):
    techCol = headers.index('TechnologyType')
    coolTechAbbrevs = {'once through':'OT','recirculating':'RC','dry cooling':'DC'}
    if row[techCol] in ptCurtailed:
        coolCol = headers.index('Cooling Tech')
        symbol = row[techCol] + '+' + coolTechAbbrevs[row[coolCol]]
    else: 
        symbol = row[techCol]
    return symbol

def separateGenSymbol(genSymbol):
    separatorIdx = genSymbol.index('+')
    return (genSymbol[:separatorIdx],genSymbol[separatorIdx+1:])

def getTechFromTechSymbol(techSymbol):
    return techSymbol.split('+')[0]

def getTechAndCoolFromTechSymbol(techSymbol):
    coolTechAbbrevsRev = {'OT':'once through','RC':'recirculating','DC':'dry'}
    return techSymbol.split('+')[0],coolTechAbbrevsRev[techSymbol.split('+')[1]]

def createZoneSymbol(zoneNum):
    return 'z' + str(zoneNum)

def getZoneNumFromSymbol(symbol):
    return int(symbol[1])

def createHourSymbol(hour):
    return 'h' + str(hour)

#If modify this function, modify next function!
def createTechAndLocLabel(tech,loc):
    return tech+'-'+loc

#If modify this function, modify prior function!
def splitTechAndLocLabel(label):
    splitLabel = label.split('-')
    return splitLabel[0],splitLabel[1]

#EXTRACT RESULTS
def extract0dVarResultsFromGAMSModel(gamsModel,varName):
    for rec in gamsModel.out_db[varName]: varValue = rec.level
    return varValue

def extract1dVarResultsFromGAMSModel(modelResults,varName):
    varResults = []
    for rec in modelResults.out_db[varName]:
        varResults.append((rec.key(0),rec.level))
    return varResults

def extract2dVarResultsIntoDict(modelResults,varName):
    varResultsDict = dict()
    for rec in modelResults.out_db[varName]:
        varResultsDict[(rec.key(0),rec.key(1))] = rec.level
    return varResultsDict

def extract2dVarResultsIntoList(modelResults,varName):
    varResults = list()
    for rec in modelResults.out_db[varName]:
        varResults.append([(rec.key(0),rec.key(1)),rec.level])
    return varResults

def extractTechBuildResults(modelResults,varName):
    varResults = list()
    for rec in modelResults.out_db[varName]:
        loc,tech = rec.key(0),rec.key(1)
        varResults.append([createTechAndLocLabel(tech,loc),rec.level])
    return varResults

def extract2dVarResultsIntoDictNoLA(modelResults,varName,hoursOpt):
    varResultsDict = dict()
    hoursOptSet = set(hoursOpt)
    for rec in modelResults.out_db[varName]:
        (gen,hour) = (rec.key(0),rec.key(1)) #Vars are indexed as egu,h
        if hour in hoursOptSet: varResultsDict[(gen,hour)] = rec.level
    return varResultsDict
