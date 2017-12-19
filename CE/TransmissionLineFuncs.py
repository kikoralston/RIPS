# Michael Craig, 12 May 2017
# Transmission line functions

import os
from AuxFuncs import *


def setLines(runLoc, capac):
    if runLoc == 'pc':
        dataDir = 'C:\\Users\\mtcraig\\Desktop\\EPP Research\\Databases\\ERCOTZonalDemandAndXmission'
    else:
        dataDir = 'Data'
    # Load data
    lines = readCSVto2dList(os.path.join(dataDir, 'ZonalTranmissionERCOT22Aug17.csv'))
    sourceCol, sinkCol = lines[0].index('Source'), lines[0].index('Sink')
    capacCol = lines[0].index('Limit')
    # Create 2 empty lists to store lines and capacities
    lineList, lineCapacs = list(), dict()
    for row in lines[1:]:
        lineName = createLineName(row[sourceCol], row[sinkCol])
        lineList.append(lineName)
        lineCapacs[lineName] = capac if capac != 'stclair' else float(row[capacCol])
    return lineList, lineCapacs


# If modify this, modify getLineSourceAndSink
def createLineName(zoneFrom, zoneTo):
    return zoneFrom + '_to_' + zoneTo


def getLineSourceAndSink(line):
    temp = line.split('_to_')
    return temp[0], temp[1]
