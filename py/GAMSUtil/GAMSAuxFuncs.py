
# Michael Craig
# October 4, 2016
# Helper functions for Python-GAMS interaction


def add0dParam(db, paramName, paramDescrip, paramValue):
    """Generic function to add a scalar (0-dimension) parameter to database

    :param db: gams database object
    :param paramName: (string) name of parameter
    :param paramDescrip: (string) description of parameter
    :param paramValue: (float) value of parameter
    """
    addedParam = db.add_parameter(paramName, 0, paramDescrip)
    addedParam.add_record().value = paramValue


def add1dParam(db, paramDict, idxSet, setSymbols, paramName, paramDescrip):
    """Generic function to add a 1-dimensional (vector) parameter to database

    :param db: gams database object
    :param paramDict: simple dictionary with data.
    :param idxSet: gams set object of domain of parameter array
    :param setSymbols: (list) list with symbol names
    :param paramName: (string) name of parameter
    :param paramDescrip: (string) description of parameter
    :return: gams parameter object created
    """
    addedParam = db.add_parameter_dc(paramName, [idxSet], paramDescrip)
    for idx in setSymbols:
        addedParam.add_record(idx).value = paramDict[idx]
    return addedParam


def add2dParam(db, param2dDict, idxSet1, idxSet2, paramName, paramDescrip):
    """Generic function to add a 2-dimensional (matrix) parameter to database

    :param db: gams database object
    :param param2dDict: dictionary with data. keys are a tuple (key1, key2)
    :param idxSet1: gams set object of first domain of parameter array
    :param idxSet2: gams set object of second domain of parameter array
    :param paramName: (string) name of parameter
    :param paramDescrip: (string) description of parameter
    :return: gams parameter object created
    """
    addedParam = db.add_parameter_dc(paramName, [idxSet1, idxSet2], paramDescrip)
    for k, v in iter(param2dDict.items()):
        addedParam.add_record(k).value = v
    return addedParam


def add_NdParam(db, paramDict, list_idxSet, paramName, paramDescrip):
    """ Generic function to add a N-dimensional parameter to database

    paramDict is a simple dictionary (it MUST NOT BE a nested dictionary). See function 'nested_dict_to_dict()'
    to convert a nested dictionary into a simple dictionary where the key is a tuple with combination of nested keys

    :param db: gams database object
    :param paramDict: simple dictionary with data. keys must be a tuple with N values
    :param list_idxSet: 1d list of size N with GAMS sets for domain of parameters (must be in the correct orders)
    :param paramName: (string) name of parameter
    :param paramDescrip: (string) description of parameter
    :return: gams parameter object
    """
    addedParam = db.add_parameter_dc(paramName, list_idxSet, paramDescrip)
    for k, v in iter(paramDict.items()):
        addedParam.add_record(k).value = v
    return addedParam


def add3dParam(db, param3dDict, idxSet1, idxSet2, idxSet3, paramName, paramDescrip):
    addedParam = db.add_parameter_dc(paramName, [idxSet1, idxSet2, idxSet3], paramDescrip)
    for k, v in iter(param3dDict.items()):
        addedParam.add_record(k).value = v
    return addedParam

def createGenSymbol(row, headers):
    """Create gen symbols

    :param row: row of genFleet 2d list with data of single generator
    :param headers: list of headers
    :return: string with gen code
    """

    (orisCol, unitCol) = (headers.index('ORIS Plant Code'),headers.index('Unit ID'))
    return str(row[orisCol]) + '+' + str(row[unitCol])


def createTechSymbol(row, headers, ptCurtailed):
    """Create cooling tech symbols

    :param row: row of genFleet 2d list with data of single generator
    :param headers: 1-d list with header
    :param ptCurtailed: 1-d list with curtailed technologies
    :return: string with symbol
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

    note that this function only works for up to 9 zones (it assumes zone number is one digit)

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

    **if modify this function, modify next function!**

    :param tech:
    :param loc:
    :return:
    """
    return tech+'-'+loc


def splitTechAndLocLabel(label):
    """split string with label of tech and location into tuple

    **if modify this function, modify prior function!**

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
    """ Extract results of 2d variable

    :param modelResults: GAMS model object
    :param varName: string with variable name
    :param hoursOpt: list with hours of optimization simulation
    :return: dictionary with results {(gen,hour): value}
    """
    varResultsDict = dict()
    hoursOptSet = set(hoursOpt)
    for rec in modelResults.out_db[varName]:
        (gen,hour) = (rec.key(0),rec.key(1)) #Vars are indexed as egu,h
        if hour in hoursOptSet: varResultsDict[(gen,hour)] = rec.level
    return varResultsDict


# Takes in gen fleet and param col name, and returns a dictionary of (genSymbol:paramVal)
# for each row.
def getEguParamDict(genFleet, paramColName, *scalar):
    paramCol = genFleet[0].index(paramColName)
    paramDict = dict()

    for row in genFleet[1:]:
        paramDict[createGenSymbol(row, genFleet[0])] = float(row[paramCol]) * scalar[0]

    return paramDict


# Takes in gen fleet and returns dictionary of (genSymbol:opCost)
def getEguOpCostDict(genFleet, scaleLbToShortTon, scaleMWtoGW, scaleDollarsToThousands, *co2Price):

    if len(co2Price[0]) > 0:
        (opCosts, hrs) = calcOpCosts(genFleet, scaleLbToShortTon, co2Price[0][0])  # thousand $/GWh
    else:
        (opCosts, hrs) = calcOpCosts(genFleet, scaleLbToShortTon)  # thousand $/GWh

    paramDict = dict()
    for idx in range(1, len(genFleet)):
        genSymb = createGenSymbol(genFleet[idx], genFleet[0])
        paramDict[genSymb] = opCosts[idx - 1] * scaleMWtoGW / scaleDollarsToThousands  # op costs = 1d list of vals, so offset by 1

    return paramDict
