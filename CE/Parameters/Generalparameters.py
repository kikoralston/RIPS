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

        self.referenceCase = False     # True: Run case without curtailments or effects of cc in demand

        self.incCurtailments = True  # whether to model curtailments,whether to model env regs on water T
        self.incRegs = True  # whether to model curtailments,whether to model env regs on water T
        self.coolDesignT = 100  # design temperature of cooling techs

        # PTs curtailed via regression (ptCurtailed) and via enviro regulations (ptCurtailedRegs)
        self.ptCurtailed = {'Coal Steam', 'Combined Cycle', 'Nuclear'}
        self.ptCurtailedRegs = {'Coal Steam', 'Coal Steam CCS', 'Combined Cycle', 'Combined Cycle CCS', 'Nuclear'}
        self.ptCurtailedAll = self.ptCurtailed | self.ptCurtailedRegs

        # Determines which cells new plants can be sited in; can be all all cells or those w/ gens already
        self.cellsEligibleForNewPlants = 'withGens'  # 'all' (all cells) or 'withGens' (only cells with gen inside)
        # Of those cells given above parameter, whether input all of them or just the one w/ max wtr T per zone.
        self.cellNewTechCriteria = 'all'  # 'all' or 'maxWaterT' (maxWaterT: only include cell w/ max water T in CE.)

        self.compressFleet = True
        self.co2CapScenario = 'none'
        self.scenario = 'FuelPriceTimeSeries.csv'

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
        self.retireOldPlants = True

        self.ucFilename = 'UCRIPS18April2017.gms'
        self.calculateCO2Price = True
        self.daysOpt = 1                    # number of days being optimized in each individual UC run
        self.daysLA = 1                     # number of "Look Ahead Days"
        self.ocAdderMin = 0
        self.ocAdderMax = 0.05              # $/MWh
        self.ucDayInitial = 1               # initial day of the year for UC simulation
        self.ucDayEnd = 365                 # final day of the year for UC simulation

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
        self.coldStart = False  # "Cold Start" for CE and UC models. Read files with initial conditions in first run
        self.gcmranking = []    # list with ranking of GCMs that will be chosen in each CE year (e.g. [3, 9, 15])
        self.rcp = ''           # name of rcp being simulated (rcp45 or rcp85)

        # OLD VARIABLES
        self.testModel = False  # use dummy test system; not currently working

    def __str__(self):

        outstr = '#\n# ------ PARAMETER FILE --------\n#\n# Description:\n#\n# End Description.\n#\n'
        outstr = outstr + 'dataRoot = {}\n'.format(self.dataRoot)
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

        outstr = outstr + 'referenceCase = {}\n'.format(self.referenceCase)

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
        outstr = outstr + 'retireOldPlants = {}\n'.format(self.retireOldPlants)

        outstr = outstr + '#\n# -------- UNIT COMMITMENT PARAMETERS --------\n#\n'
        outstr = outstr + 'ucFilename = {}\n'.format(self.ucFilename)
        outstr = outstr + 'calculateCO2Price = {}\n'.format(self.calculateCO2Price)
        outstr = outstr + 'daysOpt = {}\n'.format(self.daysOpt)
        outstr = outstr + 'daysLA = {}\n'.format(self.daysLA)
        outstr = outstr + 'ocAdderMin = {}\n'.format(self.ocAdderMin)
        outstr = outstr + 'ocAdderMax = {} \t # $/MWh\n'.format(self.ocAdderMax)
        outstr = outstr + 'ucDayInitial = {}\n'.format(self.ucDayInitial)
        outstr = outstr + 'ucDayEnd = {}\n'.format(self.ucDayEnd)

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
        outstr = outstr + 'ncores_py = {}   # number of cores to use for parallel simulation in python\n'.format(self.ncores_py)
        outstr = outstr + 'ncores_gams = {} # number of cores to use for parallel simulation in gams\n'.format(self.ncores_gams)
        outstr = outstr + 'coldStart = {}   # "Cold Start" for CE and UC models. Read files with initial conditions in first run\n'.format(self.coldStart)
        outstr = outstr + 'gcmranking = {}  # list with ranking of GCMs that will be chosen in each CE year (e.g. [3, 9, 15])\n'.format(self.list2string(self.gcmranking))
        outstr = outstr + 'rcp = {}  # name of rcp being simulated (rcp45 or rcp85)\n'.format(self.rcp)

        return outstr

    @staticmethod
    def importFuelPrices(dataRoot, scenario):
        fuelPriceDir = os.path.join(dataRoot, 'FuelPricesCapacityExpansion')

        return readCSVto2dList(os.path.join(fuelPriceDir, scenario))

    @staticmethod
    def list2string(ll):

        if ll is None:
            b = ''
        else:
            a = str(ll)
            b = ((a.replace('[', '')).replace(']', '')).replace('\'', '')

        return b

    def load(self, fname):
        """
        Reads parameter file and allocates to object

        :param fname: string with path to parameter file
        """

        i = 0
        data = {}
        with open(os.path.expanduser(fname), 'r') as csvfile:
            csvreader = csv.reader(csvfile, delimiter='=')

            for row in csvreader:
                if len(row) > 0:
                    if row[0][0] != '#':
                        # remove inline comment and trim
                        if row[1].find('#') > 0:
                            row[1] = (row[1][0:(row[1].find('#'))]).strip()
                        else:
                            row[1] = row[1].strip()

                        key = row[0].strip()
                        data[key] = row[1]

        # lev is an abbreviation for ast.literal_eval()

        for key in data:
            #print(key + ': ' + data[key])
            if key in ['ptCurtailed', 'ptCurtailedRegs']:
                setattr(self, key, set(map(str.strip, data[key].split(','))))
            elif key in ['ptEligRetCF', 'gcmranking']:
                setattr(self, key, list(map(str.strip, data[key].split(','))))
            elif key in ['dataRoot', 'resultsDir', 'pathSysGams']:
                setattr(self, key, data[key])
            else:
                try:
                    setattr(self, key, lev(data[key]))
                except (ValueError, SyntaxError):
                    setattr(self, key, data[key])

        # convert ranking of gcms to expected format
        if self.gcmranking == ['']:
            # change empty list to None
            self.gcmranking = None
        else:
            # change string values to integer
            self.gcmranking = list(map(int, self.gcmranking))

#        folderName = ('Area' + self.analysisArea + 'Cells' + self.cellNewTechCriteria +
#                      ('Curtail' if self.incCurtailments == True else 'NoCurtail') +
#                      ('EnvRegs' if self.incRegs == True else 'NoRegs') +
#                      'C' + self.co2CapScenario + 'S' + self.scenario[:3])
#        self.resultsDir = os.path.join(self.resultsDir, folderName)

        self.ptCurtailedAll = self.ptCurtailed | self.ptCurtailedRegs

        self.fuelPricesTimeSeries = self.importFuelPrices(self.dataRoot, self.scenario)

        self.setstates()

        self.maxAddedZonalCapacPerTech = self.readMaxCapacTechParam(self.maxAddedZonalCapacPerTech)

        # get IPM polygons
        self.fipsToZones, self.fipsToPolys = getIPMPolys(self.dataRoot, self.ipmZones)

        # get state polygons
        self.statePolys = getStatePolys(self.dataRoot, self.states)

        # CREATE LIST W/ IDXS THAT PARALLEL ZONES - NEED TO PRESERVE ORDER SO DON'T USE A DICT
        self.ipmZoneNums = [idx for idx in range(1, len(self.ipmZones) + 1)]

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

        if isinstance(s, str):
            if s.find(':') > -1:
                # data is a dictionary {type: maxCapacValue}
                d = Generalparameters.string2dict(s)
            else:
                # data is a single value
                d = float(s)
        else:
            # data is a single numeric value
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
