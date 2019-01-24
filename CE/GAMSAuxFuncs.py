"""
Michael Craig
October 4, 2016
Helper functions for Python-GAMS interaction
"""


def createGenSymbol(row, headers):
    """CREATE GEN SYMBOLS

    :param row: row of genFleet 2d list with data of single generator
    :param headers: list of headers
    :return: strinf with gen code
    """

    (orisCol, unitCol) = (headers.index('ORIS Plant Code'),headers.index('Unit ID'))
    return str(row[orisCol]) + '+' + str(row[unitCol])


def createTechSymbol(row, headers, ptCurtailed):
    """CREATE COOLING TECH SYMBOLS

    :param row: row of genFleet 2d list with data of single generator
    :param headers: list of headers
    :param ptCurtailed: set of curtailed technologies
    :return:
    """
    techCol = headers.index('TechnologyType')
    coolTechAbbrevs = {'once through':'OT','recirculating':'RC','dry cooling':'DC'}

    if row[techCol] in ptCurtailed:
        coolCol = headers.index('Cooling Tech')
        symbol = row[techCol] + '+' + coolTechAbbrevs[row[coolCol]]
    else: 
        symbol = row[techCol]

    return symbol


def separateGenSymbol(genSymbol):
    """ Split generator symbol into plant code and unit code

    :param genSymbol: string with gen symbol
    :return: tuple with strings of plant code and unit code
    """
    separatorIdx = genSymbol.index('+')

    return (genSymbol[:separatorIdx], genSymbol[separatorIdx+1:])


def getTechFromTechSymbol(techSymbol):
    """ Get tech type from tech symbol

    :param techSymbol: string with tech symbol
    :return: string with tech type
    """
    return techSymbol.split('+')[0]


def getTechAndCoolFromTechSymbol(techSymbol):
    """ Get tech and cooling types from tech symbol

    :param techSymbol: string with tech symbol
    :return: tuple with string of tech name and cooling type
    """
    coolTechAbbrevsRev = {'OT':'once through', 'RC':'recirculating', 'DC':'dry'}

    return techSymbol.split('+')[0], coolTechAbbrevsRev[techSymbol.split('+')[1]]


def createZoneSymbol(zoneNum):
    """ create zone symbol for GAMS (z1, z2,, ...)

    :param zoneNum: number of zone
    :return: string with zone symbol
    """
    return 'z' + str(zoneNum)


def getZoneNumFromSymbol(symbol):
    """ Get zone number from zone symbol

    NOTE THAT THIS FUNCTION ONLY WORKS FOR UP TO 9 ZONES (IT ASSUMES ZONE NUMBER IS ONE DIGIT)

    :param symbol: string with zone symbol (z1, z2, ...)
    :return: integer with number of zone
    """
    return int(symbol[1])


def createHourSymbol(hour):
    """ create hour of the year symbol for GAMS (h1, h2, h744, h8760, ...)

    :param hour: integer with hour of the year
    :return: string with hour symbol
    """
    return 'h' + str(hour)


def createTechAndLocLabel(tech,loc):
    """ create tech and location symbol for GAMS

    IF MODIFY THIS FUNCTION, MODIFY NEXT FUNCTION!

    :param tech:
    :param loc:
    :return:
    """
    return tech+'-'+loc


def splitTechAndLocLabel(label):
    """split string with label of tech and location into tuple

    IF MODIFY THIS FUNCTION, MODIFY PRIOR FUNCTION!

    :param label: label if tech and location candidate
    :return: tuple with strings of tech and location
    """
    splitLabel = label.split('-')

    return splitLabel[0],splitLabel[1]


def extract0dVarResultsFromGAMSModel(gamsModel,varName):
    """Extract results from GAMS of variable of dimension 0 (scalar)

    Reads resulting GAMS model and extracts values of variable varName

    :param gamsModel: GAMS model object
    :param varName: string with variable name
    :return: value of variable
    """

    for rec in gamsModel.out_db[varName]: varValue = rec.level

    return varValue


def extract1dVarResultsFromGAMSModel(modelResults,varName):
    """Extract results from GAMS of variable of dimension 1 (vector)

    Reads resulting GAMS model and extracts values of variable varName

    :param modelResults: GAMS model object
    :param varName: string with variable name
    :return: list with value of variable
    """
    varResults = []
    for rec in modelResults.out_db[varName]:
        varResults.append((rec.key(0),rec.level))
    return varResults


def extract2dVarResultsIntoDict(modelResults,varName):
    """Extract results from GAMS of variable of dimension 2 (2d array)

    Reads resulting GAMS model and extracts values of variable varName and saves to dictionary

    :param modelResults: GAMS model object
    :param varName: string with variable name
    :return: dict with value of variable
    """

    varResultsDict = dict()

    for rec in modelResults.out_db[varName]:
        varResultsDict[(rec.key(0),rec.key(1))] = rec.level

    return varResultsDict


def extract2dVarResultsIntoList(modelResults,varName):
    """Extract results from GAMS of variable of dimension 2 (2d array)

    Reads resulting GAMS model and extracts values of variable varName and saves to 2d list

    :param modelResults: GAMS model object
    :param varName: string with variable name
    :return: 2d list with value of variable
    """

    varResults = list()

    for rec in modelResults.out_db[varName]:
        varResults.append([(rec.key(0),rec.key(1)),rec.level])

    return varResults


def extractTechBuildResults(modelResults,varName):
    """ Extracts Tech build results

    :param modelResults: GAMS model object
    :param varName: string with variable name
    :return: 2d list with value of variable
    """
    varResults = list()

    for rec in modelResults.out_db[varName]:
        loc,tech = rec.key(0),rec.key(1)
        varResults.append([createTechAndLocLabel(tech,loc),rec.level])

    return varResults


def extract2dVarResultsIntoDictNoLA(modelResults,varName,hoursOpt):
    """ EXTRACT RESULTS

    :param modelResults:
    :param varName:
    :param hoursOpt:
    :return:
    """
    varResultsDict = dict()
    hoursOptSet = set(hoursOpt)
    for rec in modelResults.out_db[varName]:
        (gen,hour) = (rec.key(0),rec.key(1)) #Vars are indexed as egu,h
        if hour in hoursOptSet: varResultsDict[(gen,hour)] = rec.level
    return varResultsDict
