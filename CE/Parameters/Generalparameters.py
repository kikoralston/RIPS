import os
import sys
from ipmzones.AssignCellsToIPMZones import getIPMPolys, getStatePolys
from AuxFuncs import *
from ast import literal_eval as lev
import copy


class Generalparameters:
    """ Generalparameters class

    This class contains all general parameters used in the simulation. I aggregated all of them inside this class
    to make it easier to pass them as argument to the functions.

    :param dataRoot: (string) path to folder with general input data
    :param pathSysGams: (string) path to gams installation
    :param runCE: (boolean) if True run CE simulation
    :param runUC: (boolean) if True run UCED simulation
    :param runFirstUCYear: (boolean) if True run UCED in initial year of simulation
    :param startYear: (integer) first year of simulation
    :param endYear: (integer) end year of simulation
    :param yearStepCE: (integer) size of steps (in years) in CE simulation
    :param daysPerSeason: (integer) number of days in each season considered in CE simulation
    :param analysisArea: (string) analysis area
    :param useLineLimits: (boolean) if True, tries to read file with transmission limits (dataRoot/Transmission/). If False, limits = +Inf
    :param referenceCase: (boolean) If True, run case without curtailments or effects of climate change in demand
    :param incDemandCC: (boolean) If True, run case WITH effects of climate change in demand
    :param incHydroCC: (boolean) If True, run case WITH effects of climate change in hydropower
    :param incCurtailments: (boolean) whether to model curtailments,whether to model env regs on water T
    :param incRegs: (boolean) whether to model curtailments,whether to model env regs on water T
    :param ptCurtailed: (set) types of power plants for which thermal deratings simulations will be performed (default: {'Coal Steam', 'Combined Cycle', 'Nuclear'})
    :param ptCurtailedRegs: (set) types of power plants for which thermal curtailments due to enviro regulations simulations will be performed (default: {'Coal Steam', 'Coal Steam CCS', 'Combined Cycle', 'Combined Cycle CCS', 'Nuclear'})
    :param ptCurtailedAll: (set) used internally only. should not be set by user. Union of ptCurtailed and ptCurtailedRegs
    :param cellsEligibleForNewPlants: (string) Determines which cells new plants can be sited in; can be all all cells or those w/ gens already. 'all' (all cells) or 'withGens' (only cells with gen inside)
    :param cellNewTechCriteria: (string) Of those cells given above parameter, whether input all of them or just the one w/ max wtr T per zone. 'all' or 'maxWaterT' (maxWaterT: only include cell w/ max water T in CE.)
    :param compressFleet:  (boolean) Compress fleet by agregating small generators?
    :param co2CapPercentage: (float) Cap of co2 emission in final year of analysis **as percentage of emission in base year** (i.e., 100% - %reduction)
    :param scenario: (string) File with fossil fuel cost scenario (e.g. 'FuelPriceTimeSeries.csv')
    :param fuelPricesTimeSeries: (list) used internally only. should not be set by user.
    :param resultsDir: (string) path to folder where output data will be saved
    :param tzAnalysis: (string) timezone of analysis (default: 'EST')
    :param projectName: (string) project name
    :param windGenDataYr: (integer) year of wind database that will be used for wind generation simulations (e.g.: 2009, this was picked somewhat arbitrarily - average CF of all years in WIND data)
    :param capacExpFilename: (string) name of gams file with capacity expansion model (e.g., 'CERIPS24July2017PS.gms')
    :param maxAddedZonalCapacPerTech: (integer) maximum number of new techs (by type) that can be added in each CE xtep
    :param incITC: (boolean) include ITC?
    :param retirementCFCutoff: (float) retire units w/ CF lower than this value
    :param ptEligRetCF: (list) types of power plants that can be retired due to low CF
    :param selectCurtailDays: (boolean) add days with maximum curtailment to special days
    :param planningReserve: (float) planning reserve as fraction of peak demand
    :param discountRate: (float) discount rate used for computation of annualized costs
    :param allowCoalWithoutCCS: (boolean) allow Coal without CCS?
    :param onlyNSPSUnits: (boolean) allow only NSPS Units?
    :param permitOncethru: (boolean) permit adding power plants with Once-through cooling?
    :param retireOldPlants: (boolean) retire Old Plants?
    :param ucFilename: (string) 'UCRIPS18April2017.gms'
    :param calculateCO2Price: (boolean) True
    :param daysOpt: (integer) number of days being optimized in each individual UC run
    :param daysLA: (integer) number of "Look Ahead Days" in each individual UC run
    :param ocAdderMin: (float) 0
    :param ocAdderMax: (float) 0.05 $/MWh
    :param ucDayInitial: (integer) initial day of the year for complete UC annual simulation
    :param ucDayEnd: (integer) final day of the year for complete UC annual simulation
    :param costnse: (float) cost of non supplied energy (cnse) in $/MWh
    :param phEff: (float) pumped storage efficiency
    :param phMaxSoc: (float) maximum state of charge of pumped hydro. Max soc as multiple of capacity
    :param phInitSoc: (float) initial state of charge of pumped hydro. Init SOC as fraction of max soc
    :param scaleMWtoGW: 1000
    :param scaleDollarsToThousands: 1000
    :param scaleLbToShortTon: 2000
    :param states: (list) used internally only. should not be set by user.
    :param statesAbbrev: (list) used internally only. should not be set by user.
    :param ipmZones: (list) used internally only. should not be set by user.
    :param fipsToZones: (list) used internally only. should not be set by user.
    :param fipsToPolys: (list) used internally only. should not be set by user.
    :param statePolys: (list) used internally only. should not be set by user.
    :param ipmZoneNums: (list) used internally only. should not be set by user.
    :param lineList: (list) used internally only. should not be set by user.
    :param lineCapacs: (list) used internally only. should not be set by user.
    :param ncores_py: (integer) number of cores to use for parallel simulation in python
    :param ncores_gams: (integer) number of cores to use for parallel simulation in gams
    :param coldStart: (boolean) "Cold Start" for CE and UC models. Read files with initial conditions in first run
    :param gcmranking: (list) list with ranking of GCMs that will be chosen in each CE year (e.g. [3, 9, 15])
    :param rcp: (string) name of rcp being simulated (rcp45 or rcp85)
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

        self.incDemandCC = True  # whether to model CC impacts in demand
        self.incHydroCC = True  # whether to model CC impacts in hydropower

        self.incCurtailments = True  # whether to model curtailments,whether to model env regs on water T
        self.incRegs = True  # whether to model curtailments,whether to model env regs on water T

        # PTs curtailed via regression (ptCurtailed) and via enviro regulations (ptCurtailedRegs)
        self.ptCurtailed = {'Coal Steam', 'Combined Cycle', 'Nuclear'}
        self.ptCurtailedRegs = {'Coal Steam', 'Coal Steam CCS', 'Combined Cycle', 'Combined Cycle CCS', 'Nuclear'}
        self.ptCurtailedAll = self.ptCurtailed | self.ptCurtailedRegs

        # Determines which cells new plants can be sited in; can be all all cells or those w/ gens already
        self.cellsEligibleForNewPlants = 'withGens'  # 'all' (all cells) or 'withGens' (only cells with gen inside)
        # Of those cells given above parameter, whether input all of them or just the one w/ max wtr T per zone.
        self.cellNewTechCriteria = 'all'  # 'all' or 'maxWaterT' (maxWaterT: only include cell w/ max water T in CE.)

        self.compressFleet = True
        self.co2CapPercentage = float("inf")
        self.scenario = 'FuelPriceTimeSeries.csv'

        self.fuelPricesTimeSeries = []

        self.resultsDir = ''

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
        self.costnse = 5000                 # cost of non supplied energy (cnse) in $/MWh

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

        outstr = outstr + 'incDemandCC = {}\n'.format(self.incDemandCC)
        outstr = outstr + 'incHydroCC = {}\n'.format(self.incHydroCC)

        outstr = outstr + '#\n# -------- KEY CURTAILMENT PARAMETERS --------\n#\n'
        outstr = outstr + 'incCurtailments = {} \t # whether to model curtailments,whether to model env regs ' \
                          'on water T\n'.format(self.incCurtailments)
        outstr = outstr + 'incRegs = {} \t # whether to model curtailments,whether to model env regs on ' \
                          'water T\n'.format(self.incRegs)
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
        outstr = outstr + 'co2CapPercentage = {}\n'.format(self.co2CapPercentage)
        outstr = outstr + 'scenario = {}\n'.format(self.scenario)

        outstr = outstr + '#\n# -------- RESULTS DIRECTORY --------\n#\n'
        outstr = outstr + 'resultsDir = {}\n'.format(self.resultsDir)

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
        outstr = outstr + 'costnse = {}\n'.format(self.costnse)

        outstr = outstr + '#\n# -------- PUMPED HYDRO PARAMETERS --------\n#\n'
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
        Reads formatted txt file with the values of the parameters and allocates to object

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

        self.co2CapPercentage = float(self.co2CapPercentage)

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
