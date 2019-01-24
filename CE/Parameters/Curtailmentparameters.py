import csv
import os
from ast import literal_eval as lev


class Curtailmentparameters:
    """
    Curtailmentparameters class
    """

    def __init__(self):
        # Env reg parameters: dict of state : max water T (deg C)
        self.envRegMaxT = dict()
        self.maxFracFlow = dict()
        # RBM processing parameters
        self.outputHeaders = []
        self.locPrecision = 0
        self.numCellsToProcess = 0
        self.tempAndSpatFilename = ''
        self.nsegFilename = ''
        # Water T directories
        self.rbmRootDir = './'
        self.rbmDataDir = './'      # name of sub folder inside rbmRootDir
        self.rbmOutputDir = './'    # name of sub folder inside rbmRootDir
        # ----- GCMs ---------
        self.listgcms = list()
        self.basenamemeteo = 'forcing_maca_hourly_{0}_{1:4d}0101-{1:4d}1231.nc'
        self.basenamestreamT = 'serc.{}.stream_T.nc'
        self.basenameflow = 'serc.{}.KW.flow.MACA.regulated.nc'

    def __str__(self):

        strout = '#\n#  ------ CURTAILMENT PARAMETER FILE --------\n#\n# Description:\n#\n# End Description.\n#\n'

        strout = strout + '# ----- Environ regul parameters ----------\n'
        strout = strout + '# python dictionay format (with quotes and curly braces) state : max water T (deg C)\n'
        strout = strout + 'envRegMaxT = {}\n'.format(self.dict2string(self.envRegMaxT))
        strout = strout + 'python dictionary format (without quotes and curly braces) month : max fractions\n'
        strout = strout + 'maxFracFlow = {}\n'.format(self.dict2string(self.maxFracFlow))

        strout = strout + '# ----- RBM processing parameters ----\n'
        strout = strout + 'outputHeaders = {}\n'.format(self.outputHeaders)
        strout = strout + 'locPrecision = {}\n'.format(self.locPrecision)
        strout = strout + 'numCellsToProcess = {}\n'.format(self.numCellsToProcess)
        strout = strout + 'tempAndSpatFilename = {}\n'.format(self.tempAndSpatFilename)
        strout = strout + 'nsegFilename = {}\n'.format(self.nsegFilename)

        strout = strout + '# ----- Water T directories ------\n'
        strout = strout + 'rbmRootDir = {}\n'.format(self.rbmRootDir)
        strout = strout + 'rbmDataDir = {}\n'.format(self.rbmDataDir)
        strout = strout + 'rbmOutputDir = {}\n'.format(self.rbmOutputDir)
        strout = strout + '# ----- GCMs ---------\n'
        strout = strout + 'listgcms = {}\n'.format(self.listgcms)
        strout = strout + 'basenamemeteo = {}\n'.format(self.basenamemeteo)
        strout = strout + 'basenamestreamT = {}\n'.format(self.basenamestreamT)
        strout = strout + 'basenameflow = {}\n'.format(self.basenameflow)

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
        Reads curtailment parameter file and allocates to object

        :param fname: string with path to curtailment parameter file
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

        # Env regulations parameters: dict of state : max water T (deg C)
        self.envRegMaxT = self.string2dict(data[0][1])    # dict
        self.maxFracFlow = self.string2dict(data[1][1])  # dict

        # RBM processing parameters
        self.outputHeaders = list(map(str.strip, data[2][1].split(',')))  # list
        self.locPrecision = lev(data[3][1])
        self.numCellsToProcess = lev(data[4][1])
        self.tempAndSpatFilename = data[5][1]
        self.nsegFilename = data[6][1]
        # Water T directories
        self.rbmRootDir = data[7][1]
        self.rbmDataDir = data[8][1]
        self.rbmDataDir = os.path.join(self.rbmRootDir, self.rbmDataDir)
        self.rbmOutputDir = data[9][1]
        self.rbmOutputDir = os.path.join(self.rbmRootDir, self.rbmOutputDir, self.tempAndSpatFilename)
        # ----- GCMs ---------
        self.listgcms = list(map(str.strip, data[10][1].split(',')))  # list
        self.basenamemeteo = data[11][1]
        self.basenamestreamT = data[12][1]
        self.basenameflow = data[13][1]

    def writefile(self, fname):

        with open(os.path.expanduser(fname), 'w') as csvfile:
            csvfile.write(self.__str__())



