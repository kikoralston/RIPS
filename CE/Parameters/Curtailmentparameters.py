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

    def __str__(self):

        strout = '#\n#  ------ CURTAILMENT PARAMETER FILE --------\n#\n# Description:\n#\n# End Description.\n#\n'

        strout = strout + '# ----- Environ regul parameters ----------\n'
        strout = strout + '# python dictionay format (with quotes and curly braces) state : max water T (deg C)\n'
        strout = strout + 'envRegMaxT = {}\n'.format(self.envRegMaxT)

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

        # Env reg parameters: dict of state : max water T (deg C)

        self.envRegMaxT = self.string2dict(data[0][1])    # dict
        # RBM processing parameters
        self.outputHeaders = list(map(str.strip, data[1][1].split(',')))  # list
        self.locPrecision = lev(data[2][1])
        self.numCellsToProcess = lev(data[3][1])
        self.tempAndSpatFilename = data[4][1]
        self.nsegFilename = data[5][1]
        # Water T directories
        self.rbmRootDir = data[6][1]
        self.rbmDataDir = data[7][1]
        self.rbmDataDir = os.path.join(self.rbmRootDir, self.rbmDataDir)
        self.rbmOutputDir = data[8][1]
        self.rbmOutputDir = os.path.join(self.rbmRootDir, self.rbmOutputDir, self.tempAndSpatFilename)

    def writefile(self, fname):

        with open(os.path.expanduser(fname), 'w') as csvfile:
            csvfile.write(self.__str__())



