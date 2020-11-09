#Michael Craig
#November 16, 2016

from AuxFuncs import write2dListToCSV
import os


def writeHourlyResultsByPlant(genByPlant, regUpByPlant, regDownByPlant, flexByPlant, contByPlant, turnonByPlant,
                              turnoffByPlant, onOffByPlant, resultsDir, year, modelName, plantOrTech):
    """WRITE HOURLY RESULTS BY PLANT

    :param genByPlant:
    :param regUpByPlant:
    :param regDownByPlant:
    :param flexByPlant:
    :param contByPlant:
    :param turnonByPlant:
    :param turnoffByPlant:
    :param onOffByPlant:
    :param resultsDir:
    :param year:
    :param modelName:
    :param plantOrTech:
    """
    write2dListToCSV(genByPlant,os.path.join(resultsDir,'genBy' + plantOrTech + modelName + str(year) + '.csv'))
    write2dListToCSV(regUpByPlant,os.path.join(resultsDir,'regupBy' + plantOrTech + modelName + str(year) + '.csv'))
    write2dListToCSV(regDownByPlant,os.path.join(resultsDir,'regdownBy' + plantOrTech + modelName + str(year) + '.csv'))
    write2dListToCSV(flexByPlant,os.path.join(resultsDir,'flexBy' + plantOrTech + modelName + str(year) + '.csv'))
    write2dListToCSV(contByPlant,os.path.join(resultsDir,'contBy' + plantOrTech + modelName + str(year) + '.csv'))
    write2dListToCSV(turnonByPlant,os.path.join(resultsDir,'turnonBy' + plantOrTech + modelName + str(year) + '.csv')) 
    write2dListToCSV(turnoffByPlant,os.path.join(resultsDir,'turnoffBy' + plantOrTech + modelName + str(year) + '.csv')) 
    write2dListToCSV(onOffByPlant,os.path.join(resultsDir,'onOffBy' + plantOrTech + modelName + str(year) + '.csv')) 


def writeHourlyStoResults(chargeBySto, socBySto, resultsDir, year):
    """WRITE HOURLY RESULTS BY STORAGE UNITS

    :param chargeBySto:
    :param socBySto:
    :param resultsDir:
    :param year:
    """
    write2dListToCSV(chargeBySto,os.path.join(resultsDir,'chargeByStoUC' + str(year) + '.csv'))
    write2dListToCSV(socBySto,os.path.join(resultsDir,'socByStoUC' + str(year) + '.csv'))