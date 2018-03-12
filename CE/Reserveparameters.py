import csv
import os


class Reserveparameters:
    """
    Reserveparameter class
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

    def load(self, fname):
        """
        Reads reserve parameter file and allocates to object

        :param fname: string with path to reserve parameter file
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
                        print('{0:3d} : {1}'.format(i, row))
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
