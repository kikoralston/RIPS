import csv
import os


class Reserveparameters:
    """
    Reserveparameter class

    This class contains parameters to configure simulation of regulated reserves.

    :param regLoadFrac: (float) fraction of hourly load in regulated up & down reserves
    :param contLoadFrac: (float) fraction of hourly load in contingency reserves
    :param regErrorPercentile: (float) percentile of 10-m wind & 5-m solar forecast errors used in regulated reserves
    :param flexErrorPercentile: (float) percentile of hourly wind & solar forecast errors used in regulated reserves
    :param regUpCostCoeffs: (dict) dictionary with regulated reserve costs by type of plant {'Combined Cycle': 6, 'Combined Cycle CCS': 6, 'O/G Steam': 4, 'Coal Steam': 10, 'Coal Steam CCS': 10}  # $/MWh
    :param regReserveMinutes: (float) timeframe of regulated reserves. Regulated reserves must be provided w/in 5 mins
    :param flexReserveMinutes: (float) timeframe of spinning reserves. Spin reserves must be provided w/in 10 minutes
    :param contingencyReserveMinutes: (float) timeframe of contingency res (must be provided w/in 30 minutes)
    :param minutesPerHour: (float) minutes per hour
    :param rampRateToRegReserveScalar: (float) computed internally
    :param rampRateToFlexReserveScalar: (float) computed internally
    :param rampRateToContReserveScalar: (float) computed internally
    """

    def __init__(self):
        # Requirement parameters - based on WWSIS Phase 2
        self.regLoadFrac = 0.0              # frac of hourly load in reg up & down
        self.contLoadFrac = 0.0             # frac of hourly load in contingency
        self.regErrorPercentile = 0         # percentile of 10-m wind & 5-m solar forecast errors used in reg reserves
        self.flexErrorPercentile = 0        # percentile of hourly wind & solar forecast errors used in reg reserves
        # Cost coeff - from Denholm et al. 2013, val of E sto in grid apps
        self.regUpCostCoeffs = {'Combined Cycle': 6, 'Combined Cycle CCS': 6, 'O/G Steam': 4,
                                'Coal Steam': 10, 'Coal Steam CCS': 10}  # $/MWh
        # Timeframes
        self.regReserveMinutes = 0          # reg res must be provided w/in 5 mins
        self.flexReserveMinutes = 0         # spin reserves must be provided w/in 10 minutes
        self.contingencyReserveMinutes = 0  # contingency res must be provided w/in 30 minutes
        self.minutesPerHour = 60
        self.rampRateToRegReserveScalar = self.regReserveMinutes / self.minutesPerHour          # ramp rate in MW/hr
        self.rampRateToFlexReserveScalar = self.flexReserveMinutes / self.minutesPerHour        # ramp rate in MW/hr
        self.rampRateToContReserveScalar = self.contingencyReserveMinutes / self.minutesPerHour

    def __str__(self):

        strout = '#\n# ------ RESERVE PARAMETER FILE --------\n#\n# Description:\n#\n# End Description.\n#\n'

        strout = strout + '# ----- Requirement parameters (based on WWSIS Phase 2) ------\n'
        strout = strout + 'regLoadFrac = {}  # frac of hourly load in reg up & down\n'.format(self.regLoadFrac)
        strout = strout + 'contLoadFrac = {}  # frac of hourly load in contingency\n'.format(self.contLoadFrac)
        strout = strout + 'regErrorPercentile = {}  # percentile of 10-m wind & 5-m solar forecast errors ' \
                          'used in reg reserves\n'.format(self.regErrorPercentile)
        strout = strout + 'flexErrorPercentile = {}  # percentile of hourly wind & solar forecast errors ' \
                          'used in reg reserves\n'.format(self.flexErrorPercentile)
        strout = strout + '# ---- Cost coeff (from Denholm et al. 2013, val of E sto in grid apps) ----\n'

        strout = strout + 'regUpCostCoeffs = {} # $/MWh\n'.format(self.dict2string(self.regUpCostCoeffs))
        strout = strout + '# ---- Timeframes -----\n'
        strout = strout + 'regReserveMinutes = {}  # reg res must be provided w/in 5 mins\n'.format(self.regReserveMinutes)
        strout = strout + 'flexReserveMinutes = {}  # spin reserves must be provided w/in 10 minutes\n'.format(self.flexReserveMinutes)
        strout = strout + 'contingencyReserveMinutes = {}  # contingency res must be provided w/in 30 minutes\n'.format(self.contingencyReserveMinutes)
        strout = strout + 'minutesPerHour = {}\n'.format(self.minutesPerHour)

        return strout

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

    def load(self, fname):
        """
        Reads formatted txt file with the values of the parameters and allocates to object

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
                        #print('{0:3d} : {1}'.format(i, row))
                        i = i + 1

        # Requirement parameters - based on WWSIS Phase 2
        self.regLoadFrac = float(data[0][1])
        self.contLoadFrac = float(data[1][1])
        self.regErrorPercentile = float(data[2][1])
        self.flexErrorPercentile = float(data[3][1])
        # Cost coeff - from Denholm et al. 2013, val of E sto in grid apps
        self.regUpCostCoeffs = self.string2dict(data[4][1])
        # Timeframes
        self.regReserveMinutes = float(data[5][1])
        self.flexReserveMinutes = float(data[6][1])
        self.contingencyReserveMinutes = float(data[7][1])
        self.minutesPerHour = float(data[8][1])
        self.rampRateToRegReserveScalar = self.regReserveMinutes / self.minutesPerHour
        self.rampRateToFlexReserveScalar = self.flexReserveMinutes / self.minutesPerHour
        self.rampRateToContReserveScalar = self.contingencyReserveMinutes / self.minutesPerHour

    def writefile(self, fname):

        with open(os.path.expanduser(fname), 'w') as csvfile:
            csvfile.write(self.__str__())
