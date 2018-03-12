import csv
import os
from TransmissionLineFuncs import *
from AssignCellsToIPMZones import getIPMPolys, assignCellsToIPMZones
from AssignCellsToStates import getStatePolys
from AuxFuncs import *


class Generalparameters:
    """
    Generalparameters class
    """

    def __init__(self):
        self.dataRoot = ''
        self.xsedeRun = False
        self.runCE = True
        self.runUC = False
        self.runFirstUCYear = False
        self.startYear = 2015
        self.endYear = 2046
        self.yearStepCE = 10
        self.daysPerSeason = 10
        self.analysisArea = 'TVA'

        self.incCurtailments = True  # whether to model curtailments,whether to model env regs on water T
        self.incRegs = True  # whether to model curtailments,whether to model env regs on water T
        self.coolDesignT = 100  # design temperature of cooling techs
        # PTs curtailed via regression (ptCurtailed) and via enviro regulations (ptCurtailedRegs)
        self.ptCurtailed = {'Coal Steam', 'Combined Cycle'}
        self.ptCurtailedRegs = {'Coal Steam', 'Coal Steam CCS', 'Combined Cycle', 'Combined Cycle CCS', 'Nuclear'}
        # Determines which cells new plants can be sited in; can be all all cells or those w/ gens already
        self.cellsEligibleForNewPlants = 'withGens'  # 'all' (all cells) or 'withGens' (only cells with gen inside)
        # Of those cells given above parameter, whether input all of them or just the one w/ max wtr T per zone.
        self.cellNewTechCriteria = 'all'  # 'all' or 'maxWaterT' (maxWaterT: only include cell w/ max water T in CE.)

        self.compressFleet = True
        self.co2CapScenario = 'none'
        self.scenario = 'normal'

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

        self.annualDemandGrowth = 0.02  # fraction per year

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

        # OLD VARIABLES
        self.testModel = False  # use dummy test system; not currently working

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

        self.dataRoot = data[0][1]
        self.xsedeRun = data[1][1]
        self.runCE = data[2][1]
        self.runUC = data[3][1]
        self.runFirstUCYear = data[4][1]
        self.startYear = data[5][1]
        self.endYear = data[6][1]
        self.yearStepCE = data[7][1]
        self.daysPerSeason = data[8][1]
        self.analysisArea = data[9][1]

        self.incCurtailments = data[10][1]
        self.incRegs = data[11][1]
        self.coolDesignT = data[12][1]
        self.ptCurtailed = set(list(map(str.strip, data[13][1])))   # set
        self.ptCurtailedRegs = set(list(map(str.strip, data[14][1])))   #set
        self.cellsEligibleForNewPlants = data[15][1]
        self.cellNewTechCriteria = data[16][1]

        self.compressFleet = data[17][1]
        self.co2CapScenario = data[18][1]
        self.scenario = data[19][1]

        self.resultsDir = data[20][1]
        folderName = ('Area' + analysisArea + 'Cells' + cellNewTechCriteria +
                      ('Curtail' if incCurtailments == True else 'NoCurtail') +
                      ('EnvRegs' if incRegs == True else 'NoRegs') +
                      'C' + co2CapScenario + 'S' + scenario[:3])
        self.resultsDir = os.path.join(self.resultsDir, folderName)

        if not os.path.exists(resultsDir): os.makedirs(resultsDir)

        self.processRBMData = data[21][1]

        self.tzAnalysis = data[22][1]
        self.projectName = data[23][1]
        self.windGenDataYr = data[24][1]

        self.capacExpFilename = data[25][1]
        self.maxAddedZonalCapacPerTech = data[26][1]
        self.incITC = data[27][1]
        self.retirementCFCutoff = data[28][1]
        self.ptEligRetCF = list(map(str.strip, data[29][1]))    # list
        self.selectCurtailDays = data[30][1]
        self.planningReserve = data[31][1]
        self.discountRate = data[32][1]
        self.allowCoalWithoutCCS = data[33][1]
        self.onlyNSPSUnits = data[34][1]
        self.permitOncethru = data[35][1]

        self.ucFilename = data[36][1]
        self.calculateCO2Price = data[37][1]
        self.daysOpt = data[38][1]
        self.daysLA = data[39][1]
        self.ocAdderMin = data[40][1]
        self.ocAdderMax = data[41][1]

        self.phEff = data[42][1]
        self.phMaxSoc = data[43][1]
        self.phInitSoc = data[44][1]

        self.scaleMWtoGW = data[45][1]
        self.scaleDollarsToThousands = data[46][1]
        self.scaleLbToShortTon = data[47][1]

        self.annualDemandGrowth = data[48][1]

        self.setstates()

        # get IPM polygons
        self.fipsToZones, self.fipsToPolys = getIPMPolys(self.dataRoot, self.ipmZones)

        # get state polygons
        self.statePolys = getStatePolys(self.dataRoot, self.states)

        # CREATE LIST W/ IDXS THAT PARALLEL ZONES - NEED TO PRESERVE ORDER SO DON'T USE A DICT
        self.ipmZoneNums = [idx for idx in range(1, len(self.ipmZones) + 1)]
        write2dListToCSV([['Zone', 'ZoneNum']] + rotate([self.ipmZones, self.ipmZoneNums]),
                         os.path.join(resultsDir, 'zoneNamesToNumbers.csv'))

        # IMPORT TRANSMISSION SYSTEM DATA
        self.lineList, self.lineCapacs = setLines(self.ipmZones)

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
        elif self.analysisArea == 'test':
            self.states = ['North Carolina', 'South Carolina', 'Georgia', 'Mississippi', 'Alabama']
            self.statesAbbrev = ['NC', 'SC', 'GA', 'MS', 'AL']
            self.ipmZones = ['S_SOU', 'S_VACA']
        elif self.analysisArea == 'TVA':
            self.states = ['North Carolina', 'Georgia', 'Mississippi', 'Alabama', 'Kentucky', 'Tennessee']
            self.statesAbbrev = ['NC', 'GA', 'MS', 'AL', 'KY', 'TN']
            self.ipmZones = ['S_C_TVA']


