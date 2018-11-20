import os
import sys
from AssignCellsToIPMZones import getIPMPolys, assignCellsToIPMZones
from AssignCellsToStates import getStatePolys
from AuxFuncs import *
from CO2CapCalculations import getCo2Cap, interpolateCO2Cap
from ast import literal_eval as lev
import copy


class Generalparameters:
    """ Generalparameters class

    This class contains all general parameters used in the simulation. I aggregated all of them inside this class
    to make it easier to pass them as argument to the functions.

    The method load reads a formatted txt file with the values of the parameters (be careful with the order and format)

    """

    def __init__(self):
        self.dataRoot = ''
        self.pathSysGams = ''
        self.runCE = True
        self.runUC = False
        self.runFirstUCYear = False
        self.startYear = 2015
        self.endYear = 2046
        self.yearStepCE = 10
        self.daysPerSeason = 10
        self.analysisArea = 'TVA'
        self.useLineLimits = False     # True: Tries to read file with transmission limits. False: limits = +Inf

        self.incCurtailments = True  # whether to model curtailments,whether to model env regs on water T
        self.incRegs = True  # whether to model curtailments,whether to model env regs on water T
        self.coolDesignT = 100  # design temperature of cooling techs

        # PTs curtailed via regression (ptCurtailed) and via enviro regulations (ptCurtailedRegs)
        self.ptCurtailed = {'Coal Steam', 'Combined Cycle'}
        self.ptCurtailedRegs = {'Coal Steam', 'Coal Steam CCS', 'Combined Cycle', 'Combined Cycle CCS', 'Nuclear'}
        self.ptCurtailedAll = self.ptCurtailed | self.ptCurtailedRegs

        # Determines which cells new plants can be sited in; can be all all cells or those w/ gens already
        self.cellsEligibleForNewPlants = 'withGens'  # 'all' (all cells) or 'withGens' (only cells with gen inside)
        # Of those cells given above parameter, whether input all of them or just the one w/ max wtr T per zone.
        self.cellNewTechCriteria = 'all'  # 'all' or 'maxWaterT' (maxWaterT: only include cell w/ max water T in CE.)

        self.compressFleet = True
        self.co2CapScenario = 'none'
        self.scenario = 'normal'

        self.co2CapEndYr, self.co2CapEnd = None, None
        self.fuelPricesTimeSeries = []

        self.resultsDir = ''

        self.processRBMData = False  # set to True if first time improting set of RBM data

        self.tzAnalysis = 'EST'
        self.projectName = 'rips'
        self.windGenDataYr = 2009  # picked somewhat arbitrarily - average CF of all years in WIND data

        self.capacExpFilename = 'CERIPS24July2017PS.gms'
        self.maxAddedZonalCapacPerTech = 10000
        self.incITC = True
        self.retirementCFCutoff = .3  # retire units w/ CF lower than given value
        self.ptEligRetCF = ['Coal Steam']
        self.selectCurtailDays = True
        self.planningReserve = 0.15  # fraction of peak demand
        self.discountRate = 0.07  # fraction
        self.allowCoalWithoutCCS = True
        self.onlyNSPSUnits = False
        self.permitOncethru = True

        self.ucFilename = 'UCRIPS18April2017.gms'
        self.calculateCO2Price = True
        self.daysOpt = 1
        self.daysLA = 1
        self.ocAdderMin = 0
        self.ocAdderMax = 0.05  # $/MWh

        self.phEff = 0.81
        self.phMaxSoc = 5
        self.phInitSoc = .5  # max soc as multiple of capacity; init SOC as fraction of max soc

        self.scaleMWtoGW = 1000
        self.scaleDollarsToThousands = 1000
        self.scaleLbToShortTon = 2000

        self.states = []
        self.statesAbbrev = []
        self.ipmZones = []

        self.fipsToZones = None
        self.fipsToPolys = None
        self.statePolys = None

        # CREATE LIST W/ IDXS THAT PARALLEL ZONES - NEED TO PRESERVE ORDER SO DON'T USE A DICT
        self.ipmZoneNums = []

        # TRANSMISSION SYSTEM DATA
        self.lineList = None
        self.lineCapacs = None

        # OTHER
        self.ncores_py = 1      # number of cores to use for parallel simulation in python
        self.ncores_gams = 1    # number of cores to use for parallel simulation in gams

        # OLD VARIABLES
        self.testModel = False  # use dummy test system; not currently working

    def __str__(self):

        outstr = '#\n# ------ PARAMETER FILE --------\n#\n# Description:\n#\n# End Description.\n#\n'
        outstr = outstr + 'datadir = {}\n'.format(self.dataRoot)
        outstr = outstr + 'pathSysGams = {}\n'.format(self.pathSysGams)
        outstr = outstr + 'runCE = {}\n'.format(self.runCE)
        outstr = outstr + 'runUC = {}\n'.format(self.runUC)
        outstr = outstr + 'runFirstUCYear = {}\n'.format(self.runFirstUCYear)
        outstr = outstr + 'startYear = {}\n'.format(self.startYear)

        outstr = outstr + 'endYear = {}\n'.format(self.endYear)
        outstr = outstr + 'yearStepCE = {}\n'.format(self.yearStepCE)
        outstr = outstr + 'daysPerSeason = {}\n'.format(self.daysPerSeason)
        outstr = outstr + 'analysisArea = {}\n'.format(self.analysisArea)
        outstr = outstr + 'useLineLimits = {}\n'.format(self.useLineLimits)

        outstr = outstr + '#\n# -------- KEY CURTAILMENT PARAMETERS --------\n#\n'
        outstr = outstr + 'incCurtailments = {} \t # whether to model curtailments,whether to model env regs ' \
                          'on water T\n'.format(self.incCurtailments)
        outstr = outstr + 'incRegs = {} \t # whether to model curtailments,whether to model env regs on ' \
                          'water T\n'.format(self.incRegs)
        outstr = outstr + 'coolDesignT = {} \t # design temperature of cooling techs\n'.format(self.coolDesignT)
        outstr = outstr + '# PTs curtailed via regression (ptCurtailed) and via enviro regulations (ptCurtailedRegs)\n'
        outstr = outstr + 'ptCurtailed = {}\n'.format(self.list2string(list(self.ptCurtailed)))
        outstr = outstr + 'ptCurtailedRegs = {}\n'.format(self.list2string(list(self.ptCurtailedRegs)))
        outstr = outstr + '# Determines which cells new plants can be sited in; can be all all cells or those w/ ' \
                          'gens already\n'
        outstr = outstr + 'cellsEligibleForNewPlants = {} \t # \'all\' (all cells) or \'withGens\' (only cells ' \
                          'with gen inside)\n'.format(self.cellsEligibleForNewPlants)
        outstr = outstr + '# Of those cells given above parameter, whether input all of them or just the one w/ max ' \
                          'wtr T per zone.\n'
        outstr = outstr + 'cellNewTechCriteria = {} \t # \'all\' or \'maxWaterT\' (maxWaterT: only include cell w/ ' \
                          'max water T in CE.\n'.format(self.cellNewTechCriteria)

        outstr = outstr + '#\n# -------- GENERAL PARAMETERS --------\n#\n'
        outstr = outstr + 'compressFleet = {}\n'.format(self.compressFleet)
        outstr = outstr + 'co2CapScenario = {}\n'.format(self.co2CapScenario)
        outstr = outstr + 'scenario = {}\n'.format(self.scenario)

        outstr = outstr + '#\n# -------- RESULTS DIRECTORY --------\n#\n'
        outstr = outstr + 'resultsDir = {}\n'.format(self.resultsDir)

        outstr = outstr + '#\n# -------- THERMAL CURTAILMENT PARAMETERS --------\n#\n'
        outstr = outstr + 'processRBMData = {} \t # set to True if first time improting ' \
                          'set of RBM data\n'.format(self.processRBMData)

        outstr = outstr + '#\n# -------- RENEWABLE CAPACITY FACTOR PARAMETERS --------\n#\n'
        outstr = outstr + 'tzAnalysis = {}\n'.format(self.tzAnalysis)
        outstr = outstr + 'projectName = {}\n'.format(self.projectName)
        outstr = outstr + 'windGenDataYr = {} \t # picked somewhat arbitrarily - average CF of all years ' \
                          'in WIND data\n'.format(self.windGenDataYr)

        outstr = outstr + '#\n# -------- CAPACITY EXPANSION PARAMETERS --------\n#\n'
        outstr = outstr + 'capacExpFilename = {}\n'.format(self.capacExpFilename)
        outstr = outstr + 'maxAddedZonalCapacPerTech = {}\n'.format(self.writeMaxCapacTechParam())
        outstr = outstr + 'incITC = {}\n'.format(self.incITC)
        outstr = outstr + 'retirementCFCutoff = {} \t # retire units w/ CF lower than given ' \
                          'value\n'.format(self.retirementCFCutoff)
        outstr = outstr + 'ptEligRetCF = {}\n'.format(self.list2string(self.ptEligRetCF))
        outstr = outstr + 'selectCurtailDays = {}\n'.format(self.selectCurtailDays)
        outstr = outstr + 'planningReserve = {} \t # fraction of peak demand\n'.format(self.planningReserve)
        outstr = outstr + 'discountRate = {}\n'.format(self.discountRate)
        outstr = outstr + 'allowCoalWithoutCCS = {}\n'.format(self.allowCoalWithoutCCS)
        outstr = outstr + 'onlyNSPSUnits = {}\n'.format(self.onlyNSPSUnits)
        outstr = outstr + 'permitOncethru = {}\n'.format(self.permitOncethru)

        outstr = outstr + '#\n# -------- UNIT COMMITMENT PARAMETERS --------\n#\n'
        outstr = outstr + 'ucFilename = {}\n'.format(self.ucFilename)
        outstr = outstr + 'calculateCO2Price = {}\n'.format(self.calculateCO2Price)
        outstr = outstr + 'daysOpt = {}\n'.format(self.daysOpt)
        outstr = outstr + 'daysLA = {}\n'.format(self.daysLA)
        outstr = outstr + 'ocAdderMin = {}\n'.format(self.ocAdderMin)
        outstr = outstr + 'ocAdderMax = {} \t # $/MWh\n'.format(self.ocAdderMax)
        outstr = outstr +  '#\n# -------- PUMPED HYDRO PARAMETERS --------\n#\n'
        outstr = outstr + 'phEff = {}\n'.format(self.phEff)
        outstr = outstr + 'phMaxSoc = {}\n'.format(self.phMaxSoc)
        outstr = outstr + 'phInitSoc = {} \t # max soc as multiple of capacity; init SOC as fraction ' \
                          'of max soc\n'.format(self.phInitSoc)
        outstr = outstr + '#\n# -------- CONVERSION PARAMETERS --------\n#\n'
        outstr = outstr + 'scaleMWtoGW = {}\n'.format(self.scaleMWtoGW)
        outstr = outstr + 'scaleDollarsToThousands = {}\n'.format(self.scaleDollarsToThousands)
        outstr = outstr + 'scaleLbToShortTon = {}\n'.format(self.scaleLbToShortTon)
        outstr = outstr + '#\n# -------- OTHER PARAMETERS --------\n#\n'
        outstr = outstr + 'ncores_py = {} \t    # number of cores to use for parallel simulation in python'.format(self.ncores_py)
        outstr = outstr + 'ncores_gams = {} \t  # number of cores to use for parallel simulation in gams'.format(self.ncores_gams)

        return outstr

    @staticmethod
    def importFuelPrices(dataRoot, scenario):
        fuelPriceDir = os.path.join(dataRoot, 'FuelPricesCapacityExpansion')

        if scenario == 'ng':
            fuelFileName = 'FuelPriceTimeSeries2Aug2016LowNG.csv'
        else:
            fuelFileName = 'FuelPriceTimeSeries2Aug2016.csv'

        return readCSVto2dList(os.path.join(fuelPriceDir, fuelFileName))

    @staticmethod
    def list2string(ll):

        a = str(ll)
        b = ((a.replace('[', '')).replace(']', '')).replace('\'', '')

        return b

    def load(self, fname):
        """
        Reads parameter file and allocates to object

        :param fname: string with path to parameter file
        """

        i = 0
        data = []
        with open(os.path.expanduser(fname), 'r') as csvfile:
            spamreader = csv.reader(csvfile, delimiter='=')

            for row in spamreader:
                if len(row) > 0:
                    if row[0][0] != '#':
                        # remove inline comment and trim
                        if row[1].find('#') > 0:
                            row[1] = (row[1][0:(row[1].find('#'))]).strip()
                        else:
                            row[1] = row[1].strip()
                        data.append(row)
                        # print('{0:3d} : {1}'.format(i, row))
                        # i = i + 1

        # lev is an abbreviation for ast.literal_eval()

        self.dataRoot = data[0][1]
        self.pathSysGams = data[1][1]
        self.runCE = lev(data[2][1])
        self.runUC = lev(data[3][1])
        self.runFirstUCYear = lev(data[4][1])
        self.startYear = lev(data[5][1])
        self.endYear = lev(data[6][1])
        self.yearStepCE = lev(data[7][1])
        self.daysPerSeason = lev(data[8][1])
        self.analysisArea = data[9][1]
        self.useLineLimits = data[10][1]

        self.incCurtailments = lev(data[11][1])
        self.incRegs = lev(data[12][1])
        self.coolDesignT = data[13][1]

        self.ptCurtailed = set(map(str.strip, data[14][1].split(',')))  # set
        self.ptCurtailedRegs = set(map(str.strip, data[15][1].split(',')))  # set
        self.ptCurtailedAll = self.ptCurtailed | self.ptCurtailedRegs

        self.cellsEligibleForNewPlants = data[16][1]
        self.cellNewTechCriteria = data[17][1]

        self.compressFleet = lev(data[18][1])
        self.co2CapScenario = data[19][1]
        self.scenario = data[20][1]

        self.co2CapEndYr, self.co2CapEnd = getCo2Cap(self.co2CapScenario)
        self.fuelPricesTimeSeries = self.importFuelPrices(self.dataRoot, self.scenario)

        self.resultsDir = data[21][1]
        folderName = ('Area' + self.analysisArea + 'Cells' + self.cellNewTechCriteria +
                      ('Curtail' if self.incCurtailments == True else 'NoCurtail') +
                      ('EnvRegs' if self.incRegs == True else 'NoRegs') +
                      'C' + self.co2CapScenario + 'S' + self.scenario[:3])
        self.resultsDir = os.path.join(self.resultsDir, folderName)

        if not os.path.exists(self.resultsDir):
            os.makedirs(self.resultsDir)

        self.processRBMData = lev(data[22][1])

        self.tzAnalysis = data[23][1]
        self.projectName = data[24][1]
        self.windGenDataYr = lev(data[25][1])

        self.capacExpFilename = data[26][1]
        #self.maxAddedZonalCapacPerTech = lev(data[27][1])
        self.maxAddedZonalCapacPerTech = self.readMaxCapacTechParam(data[27][1])
        self.incITC = lev(data[28][1])
        self.retirementCFCutoff = lev(data[29][1])
        self.ptEligRetCF = list(map(str.strip, data[30][1].split(',')))  # list
        self.selectCurtailDays = lev(data[31][1])
        self.planningReserve = lev(data[32][1])
        self.discountRate = lev(data[33][1])
        self.allowCoalWithoutCCS = lev(data[34][1])
        self.onlyNSPSUnits = lev(data[35][1])
        self.permitOncethru = lev(data[36][1])

        self.ucFilename = data[37][1]
        self.calculateCO2Price = lev(data[38][1])
        self.daysOpt = lev(data[39][1])
        self.daysLA = lev(data[40][1])
        self.ocAdderMin = lev(data[41][1])
        self.ocAdderMax = lev(data[42][1])

        self.phEff = lev(data[43][1])
        self.phMaxSoc = lev(data[44][1])
        self.phInitSoc = lev(data[45][1])

        self.scaleMWtoGW = lev(data[46][1])
        self.scaleDollarsToThousands = lev(data[47][1])
        self.scaleLbToShortTon = lev(data[48][1])

        self.ncores_py = lev(data[49][1])
        self.ncores_gams = lev(data[50][1])

        self.setstates()

        # get IPM polygons
        self.fipsToZones, self.fipsToPolys = getIPMPolys(self.dataRoot, self.ipmZones)

        # get state polygons
        self.statePolys = getStatePolys(self.dataRoot, self.states)

        # CREATE LIST W/ IDXS THAT PARALLEL ZONES - NEED TO PRESERVE ORDER SO DON'T USE A DICT
        self.ipmZoneNums = [idx for idx in range(1, len(self.ipmZones) + 1)]
        write2dListToCSV([['Zone', 'ZoneNum']] + rotate([self.ipmZones, self.ipmZoneNums]),
                         os.path.join(self.resultsDir, 'zoneNamesToNumbers.csv'))

        # IMPORT TRANSMISSION SYSTEM DATA
        self.setLines()

    def setstates(self):
        if self.analysisArea == 'coreSERC':
            self.states = ['North Carolina', 'South Carolina', 'Georgia', 'Mississippi', 'Alabama',
                           'Kentucky', 'Tennessee']
            self.statesAbbrev = ['NC', 'SC', 'GA', 'MS', 'AL', 'KY', 'TN']
            self.ipmZones = ['S_SOU', 'S_VACA', 'S_C_KY', 'S_C_TVA']
        elif self.analysisArea == 'allSERC':
            self.states = ['North Carolina', 'South Carolina', 'Georgia', 'Mississippi', 'Alabama',
                           'Kentucky', 'Tennessee', 'Missouri', 'Virginia', 'Illinois', 'Louisiana']
            self.statesAbbrev = ['NC', 'SC', 'GA', 'MS', 'AL', 'KY', 'TN', 'MO', 'VA', 'IL', 'LA']
            self.ipmZones = ['S_SOU', 'S_VACA', 'S_C_KY', 'S_C_TVA', 'S_D_WOTA', 'S_D_N_AR', 'S_D_AMSO', 'S_D_REST']
        elif self.analysisArea == 'onlyTN':
            self.states = ['Tennessee']
            self.statesAbbrev = ['TN']
            self.ipmZones = ['S_C_TVA']
        elif self.analysisArea == 'TVA':
            self.states = ['North Carolina', 'Georgia', 'Mississippi', 'Alabama', 'Kentucky', 'Tennessee']
            self.statesAbbrev = ['NC', 'GA', 'MS', 'AL', 'KY', 'TN']
            self.ipmZones = ['S_C_TVA']
        elif self.analysisArea == 'S_C_KY':
            self.states = ['Kentucky']
            self.statesAbbrev = ['KY']
            self.ipmZones = ['S_C_KY']
        elif self.analysisArea == 'test':
            self.states = ['North Carolina', 'South Carolina', 'Georgia', 'Mississippi', 'Alabama']
            self.statesAbbrev = ['NC', 'SC', 'GA', 'MS', 'AL']
            self.ipmZones = ['S_SOU', 'S_VACA']
        else:
            print('------------------------------------------------------------------------------------------')
            print('ERROR!!!')
            print('Analysis area is set to {0} which is not a valid option!'.format(self.analysisArea))
            print('Valid analysis areas: coreSERC, allSERC, onlyTN, test, TVA')
            print('Please change analysis area to a valid option in the general parameters file.')
            print('------------------------------------------------------------------------------------------')
            sys.exit()

    def writefile(self, fname):

        with open(os.path.expanduser(fname), 'w') as csvfile:
            csvfile.write(self.__str__())

    def setLines(self):

        # set path with transmission data (dataRoot/Transmission)
        dataDir = os.path.join(self.dataRoot, 'Transmission')

        # Create 2 empty lists to store lines and capacities
        self.lineList, self.lineCapacs = list(), dict()

        # Load data
        if self.useLineLimits:
            if os.path.isfile(os.path.join(dataDir, 'ZonalTransmission.csv')):
                lines = readCSVto2dList(os.path.join(dataDir, 'ZonalTransmission.csv'))

                sourceCol, sinkCol = lines[0].index('Source'), lines[0].index('Sink')
                capacCol = lines[0].index('Limit')

                for row in lines[1:]:
                    lineName = self.createLineName(row[sourceCol], row[sinkCol])
                    self.lineList.append(lineName)

                    self.lineCapacs[lineName] = float(row[capacCol])

            else:
                print()
                print('File ZonalTransmission.csv not found in {}'.format(dataDir))
                print('Transmission upper bounds will be set to +INF.')
                print('Changing useLineLimits parameter to FALSE.')
                print()
                self.useLineLimits = False

        if not self.useLineLimits:
            for s in self.ipmZones:
                sinks = copy.deepcopy(self.ipmZones)
                sinks.remove(s)

                for t in sinks:
                    lineName = self.createLineName(s, t)
                    self.lineList.append(lineName)
                    self.lineCapacs[lineName] = float('inf')

    def writeMaxCapacTechParam(self):

        if isinstance(self.maxAddedZonalCapacPerTech, dict):
            # data is a dictionary {type: maxCapacValue}
            s = Generalparameters.dict2string(self.maxAddedZonalCapacPerTech)
        else:
            # data is a single value
            s = str(self.maxAddedZonalCapacPerTech)

        return s

    @staticmethod
    def createLineName(zoneFrom, zoneTo):
        # If modify this, modify getLineSourceAndSink
        return zoneFrom + '_to_' + zoneTo

    @staticmethod
    def getLineSourceAndSink(line):
        temp = line.split('_to_')
        return temp[0], temp[1]

    @staticmethod
    def readMaxCapacTechParam(s):

        if s.find(':') > -1:
            # data is a dictionary {type: maxCapacValue}
            d = Generalparameters.string2dict(s)
        else:
            # data is a single value
            d = float(s)

        return d

    @staticmethod
    def string2dict(s):

        a = list(map(str.strip, s.split(',')))
        b = dict()

        for x in a:
            c = x.split(':')
            k = c[0].strip()  # key
            v = float(c[1].strip())  # value

            b.update({k: v})

        return b

    @staticmethod
    def dict2string(d):

        a = str(d)
        b = ((a.replace('{', '')).replace('}', '')).replace('\'', '')

        return b







