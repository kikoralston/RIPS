#Michael Craig
#August 19, 2016
#Import and process water availability data

import os, csv, datetime
from AuxFuncs import *
from ModifyGeneratorCapacityWithWaterTData18Aug2016 import getGenToCellAndCellToGenDictionaries

################################################################################
####### SET PARAMETERS  ########################################################
################################################################################
def setParams():
    #Water avail params
    processAvailData = False
    baseDir = 'C:\\Users\\mtcraig\\Desktop\\EPP Research\\UWData\\DataFromUW\\DummyDataWithScripts9Aug2016'
    dataDir = os.path.join(baseDir,'RBMRawWaterAvailabilityData19Aug2016')
    networkFilename = 'bcc-csm1-1-m_rcp45_r1i1p1_Network'
    flowFilename = 'bcc-csm1-1-m_rcp45_r1i1p1.DA_flow'
    flowBaseFilename = flowFilename[:flowFilename.index('.')]
    outputFolder = os.path.join(baseDir,'RBMProcessedWaterAvail26Aug2016',flowBaseFilename)
    locPrecision = 4
    tgtYr = 2015
    #Fleet params
    genFleetDir = 'C:\\Users\\mtcraig\\Desktop\\EPP Research\\PythonCCPaper'
    genFleetName = 'genFleetBeforeCE.csv'
    genFleet = readCSVto2dList(os.path.join(genFleetDir,genFleetName))
    return (processAvailData,dataDir,networkFilename,flowFilename,outputFolder,locPrecision,tgtYr,genFleet)
################################################################################
################################################################################
################################################################################

################################################################################
####### MASTER FUNCTION  #######################################################
################################################################################
def masterFunction():
    (processAvailData,dataDir,networkFilename,flowFilename,outputFolder,locPrecision,
                            tgtYr,genFleet) = setParams()
    if processAvailData==True:
        processWaterAvailDataForAllCells(dataDir,networkFilename,flowFilename,outputFolder,locPrecision)
    genToWaterAvailInTgtYr = getWaterAvailDataForFleetInTgtYear(genFleet,tgtYr,outputFolder,locPrecision)
################################################################################
################################################################################
################################################################################

################################################################################
####### GET AND SAVE  WATER AVAILABILITY DATA FOR ALL CELLS  ###################
################################################################################
def processWaterAvailDataForAllCells(dataDir,networkFilename,flowFilename,outputFolder,locPrecision):
    (mapNodesToLatLong,mapLatLongToNode) = loadNetworkData(dataDir,networkFilename)
    if not os.path.exists(outputFolder): os.makedirs(outputFolder)
    numCellsToProcessSimult = 4
    dayToDates = createDatesList() #1d list of datetimes from 1950/1/1 - 2100/12/31
    (cellCtr,totalCellCtr,tgtCellLatLongs) = (0,0,set())
    for (lat,lon) in mapLatLongToNode:
        cellCtr += 1
        totalCellCtr += 1
        tgtCellLatLongs.add((lat,lon))
        if cellCtr == numCellsToProcessSimult or totalCellCtr == len(mapLatLongToNode):    
            # print('Processing cells:',tgtCellLatLongs)       
            processAndSaveWaterAvailData(dataDir,outputFolder,flowFilename,
                                        mapNodesToLatLong,locPrecision,tgtCellLatLongs,dayToDates)   
            (cellCtr,tgtCellLatLongs) = (0,set())
            print('Cells processed:',totalCellCtr)
    print('Num cells processed:',totalCellCtr)

#Load network data, returning dicts of nodes to cell (lat,lon) and vice versa.
def loadNetworkData(dataDir,networkFilename):
    (mapNodesToLatLong,mapLatLongToNode) = (dict(),dict())
    networkFile = os.path.join(dataDir,networkFilename)
    f = open(networkFile, 'r')
    while 1:
        line = f.readline().rstrip("\n")
        if line=='': break
        lineSplit = line.split()
        if 'Node' in lineSplit: #if not Node row, skip it
            #Node, lat, long all given in column after label
            node = int(lineSplit[lineSplit.index('Node') + 1])
            lat = float(lineSplit[lineSplit.index('Lat') + 1])
            lon = float(lineSplit[lineSplit.index('Long') + 1])
            mapNodesToLatLong[node] = (lat,lon)
            if (lat,lon) in mapLatLongToNode: mapLatLongToNode[(lat,lon)].append(node)
            else: mapLatLongToNode[(lat,lon)] = [node]
    f.close()
    print('Num cells with water availability data:',len(mapLatLongToNode))
    return (mapNodesToLatLong,mapLatLongToNode)

#Get empty 1d list of datetime dates from 1950 to 2100
def createDatesList():
    startSim = datetime.date(1950,1,1)
    endSim = datetime.date(2100,12,31)
    (startOrdinal,endOrdinal) = (startSim.toordinal(),endSim.toordinal())
    numDays = endOrdinal - startOrdinal + 1
    datesList = [startSim + datetime.timedelta(days=x) for x in range(0,numDays)]
    datesListAsDay = [day for day in range(1,numDays+1)]
    dayToDates = dict()
    for idx in range(len(datesListAsDay)): dayToDates[datesListAsDay[idx]] = datesList[idx]
    return dayToDates

#Process then save water avail (Q_out) data for set of target cells
def processAndSaveWaterAvailData(dataDir,outputFolder,flowFilename,
                                mapNodesToLatLong,locPrecision,tgtCellLatLongs,dayToDates):
    latLongToWaterAvail = processWaterAvailData(dataDir,flowFilename,mapNodesToLatLong,
                                                    tgtCellLatLongs,dayToDates)
    writeWaterAvailData(latLongToWaterAvail,outputFolder,locPrecision)

#Flow file cols: nd (day of simulation), node, Q_in (cfs), Q_out (cfs), Q_diff (cfs), 
#depth (ft), width (ft), u (ft/sec)
def processWaterAvailData(dataDir,flowFilename,mapNodesToLatLong,tgtCellLatLongs,dayToDates):
    latLongToWaterAvail = dict()
    for (lat,lon) in tgtCellLatLongs: latLongToWaterAvail[(lat,lon)] = [['Datetime','QOut(cfs)']]
    (dayCol,nodeCol,qOutCol) = (0,1,3)
    flowFile = os.path.join(dataDir,flowFilename)
    f = open(flowFile, 'r')
    while 1:
        line = f.readline().rstrip('\n')
        if line == '': break
        lineSplit = line.split()
        (node,qOut,day) = (int(lineSplit[nodeCol]),float(lineSplit[qOutCol]),int(lineSplit[dayCol]))
        (nodeLat,nodeLon) = mapNodesToLatLong[node]
        if (nodeLat,nodeLon) in tgtCellLatLongs:
            currDate = dayToDates[day]
            latLongToWaterAvail[(nodeLat,nodeLon)].append([currDate,qOut])
    f.close()
    return latLongToWaterAvail

#Write availability data
def writeWaterAvailData(latLongToWaterAvail,outputFolder,locPrecision):
    for (lat,lon) in latLongToWaterAvail:
        waterAvail = latLongToWaterAvail[(lat,lon)]
        cellFilename = createWaterAvailFilenameToReadOrWrite(locPrecision,lat,lon)
        write2dListToCSV(waterAvail,os.path.join(outputFolder,cellFilename))

def createWaterAvailFilenameToReadOrWrite(locPrecision, inputLat, inputLong):
    return '%.*f_%.*fWaterAvail.csv' %(locPrecision, inputLat, locPrecision, inputLong)
################################################################################
################################################################################
################################################################################


################################################################################
####### MAP GENERATORS TO AND FROM RBM GRID CELLS  #############################
################################################################################
def getWaterAvailDataForFleetInTgtYear(genFleet,tgtYr,waterAvailDir,locPrecision):
    (genToWaterAvailInTgtYr,genToWaterAvailInTgtYrList) = (dict(),[])
    (genToCellLatLongsDictValues,genToCellLatLongsDict,cellLatLongToGenDict,
                genToCellLatLongsList) = getGenToCellAndCellToGenDictionaries(genFleet)
    for (cellLat,cellLon) in cellLatLongToGenDict:
        gensInCell = cellLatLongToGenDict[(cellLat,cellLon)]
        cellWaterAvailTgtYr = loadCellWaterAvailForTgtYr(cellLat,cellLon,locPrecision,
                                                            waterAvailDir,tgtYr) #1d list [qout(cfs)]
        for gen in gensInCell:
            genToWaterAvailInTgtYr[gen] = cellWaterAvailTgtYr #1d list
            genToWaterAvailInTgtYrList.append([gen] + cellWaterAvailTgtYr) #1d list
    write2dListToCSV(genToWaterAvailInTgtYrList,'waterAvailByGenInYr' + str(tgtYr) + '.csv')    
    return genToWaterAvailInTgtYr

def loadCellWaterAvailForTgtYr(cellLat,cellLon,locPrecision,waterAvailDir,tgtYr):
    cellWaterAvail = loadCellWaterAvail(cellLat,cellLon,locPrecision,waterAvailDir) #2d list: [datetime,qout(cfs)]
    print('**str or datetime',cellWaterAvail[1][0],type(cellWaterAvail[1][0]))
    cellWaterAvailTgtYr = isolateCellWaterAvailInTgtYr(cellWaterAvail,tgtYr) #1d list [qout] for tgt yr
    return cellWaterAvailTgtYr

def loadCellWaterAvail(cellLat,cellLon,locPrecision,waterAvailDir):
    cellFilename = createWaterAvailFilenameToReadOrWrite(locPrecision,cellLat,cellLon)
    return readCSVto2dList(os.path.join(waterAvailDir,cellFilename))

def isolateCellWaterAvailInTgtYr(cellWaterAvail,tgtYr):
    (dateCol,qoutCol) = (cellWaterAvail[0].index('Datetime'),cellWaterAvail[0].index('QOut(cfs)'))
    cellWaterAvailTgtYr = [row[qoutCol] for row in cellWaterAvail[1:] if str(tgtYr) in row[dateCol]]
    return cellWaterAvailTgtYr #1d list [qout]
################################################################################
################################################################################
################################################################################

################################################################################
####### TEST FUNCTION  #########################################################
################################################################################
#Flow file cols: nd (day of simulation), node, Q_in (cfs), Q_out (cfs), Q_diff (cfs), 
#depth (ft), width (ft), u (ft/sec)
def testProcessWaterAvailData():
    dataDir = 'C:\\Users\\mtcraig\\Desktop\\EPP Research\\UWData\\DataFromUW\\DummyDataWithScripts9Aug2016\\RBMRawWaterAvailabilityData19Aug2016'
    flowFilename = 'bcc-csm1-1-m_rcp45_r1i1p1.DA_flow'
    flowFile = os.path.join(dataDir,flowFilename)
    f = open(flowFile, 'r')
    i = 0
    test2dList = []
    while i<200000:
        line = f.readline().rstrip('\n')
        if line == '': break
        lineSplit = line.split()
        test2dList.append(lineSplit)
        i += 1
    f.close()
    write2dListToCSV(test2dList,'SliceOfWaterAvailData.csv')

# testProcessWaterAvailData()
################################################################################
################################################################################
################################################################################

masterFunction()
