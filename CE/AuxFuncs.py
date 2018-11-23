# Michael Craig
# October 4, 2016
# Define several auxiliary functions that are used by a range of scripts.

import csv
import time
import numpy as np


################ READ AND WRITE CSV LISTS
# Read CSV to 2d list
# Input: full file name including dir (str)
# Output: 2d list
def readCSVto2dList(fileNameWithDir):
    with open(fileNameWithDir, 'r') as f:
        f = csv.reader(f)
        f = list(f)
    return f


def write2dListToCSV(list2d, fileNameWithDir):
    """Write 2d list to CSV file

    :param list2d: a 2d list
    :param fileNameWithDir: string with full file name including dir & .csv
    """

    fullFileName = fileNameWithDir
    with open(fullFileName, 'w', newline='') as csvfile:
        w = csv.writer(csvfile)
        w.writerows(list2d)


def convert_dict2list2d(dict_in):
    """ Convert dictionary to 2d list

    :param dict_in: a dictionary (possibly nested). A n-level nested dictionary will have all its keys in
                    the first n columns

    :return: a 2d list with the keys in the first n columns (where n is the "level" of the nested dictionary)
    """

    list2d = []

    for key in dict_in:
        if isinstance(dict_in[key], dict):
            aux = convert_dict2list2d(dict_in[key])
            list2d = list2d + [[key] + row for row in aux]
        elif isinstance(dict_in[key], list):
            list2d = list2d + [[key] + dict_in[key]]
        elif isinstance(dict_in[key], np.ndarray):
            list2d = list2d + [[key] + list(dict_in[key])]
        else:
            list2d = list2d + [[key] + [dict_in[key]]]

    return list2d


def dict_tuples_to_list2d(dict_in):
    """ Convert a dictionary that has tuples as keys to 2d list where the first columns are the separated keys

    :param dict_in: a dictionary (not nested) that has keys that are tuples

    :return: a 2d list with the combined keys (in the tuples) in the first n columns
            (where n is the length of the tuples)
    """

    dict_keys = list(dict_in.keys())
    len_list = len(dict_keys)

    list2d = [0]*len_list

    for (i, key) in enumerate(dict_keys):
        list2d[i] = list(key) + [dict_in[key]]

    return list2d


def nested_dict_to_dict(dict_in):
    """ Convert nested dictionary to simple dictionary (all keys are combined into a tuple as the resulting key)

    :param dict_in: a dictionary (possibly nested).

    :return: a simple dictionary where keys from original nested dictionary are combined in a tuple
    """

    dict_out = dict()

    for key in dict_in:
        if isinstance(dict_in[key], dict):
            aux = nested_dict_to_dict(dict_in[key])
            for key2 in aux:
                if isinstance(key2, tuple):
                    keyfinal = tuple([key] + [k for k in key2])
                else:
                    keyfinal = (key, key2)

                dict_out[keyfinal] = aux[key2]
        else:
            dict_out[key] = dict_in[key]

    return dict_out


def writeDictToCSV(dictWrite, fileNameWithDir):
    """ Writes a dictionary to a csv file

    Works with a 2 level nested dictionary. **Will not work with deeper dictionaries**. For 1 level dictionaries,
    the first column are the keys. For 2 level nested dictionaries, the first two columns are the keys.

    :param dictWrite: dictionary to write
    :param fileNameWithDir: string with file name
    """
    list2d = convert_dict2list2d(dictWrite)

    write2dListToCSV(list2d, fileNameWithDir)


################ CONVERT DOLLAR YEARS
# Convert dollar years to 2012 dollars
# CPI from Minneapolis Fed, https://www.minneapolisfed.org/community/teaching-aids/
# cpi-calculator-information/consumer-price-index-and-inflation-rates-1913.
# Inputs: name of parameter of dollar year to convert, parameter value (cost)
# Outputs: cost in 2012 dollars
def convertCostToTgtYr(paramName, cost):
    paramDollarYears = {'startup': 2011, 'vom': 2012, 'fom': 2012, 'occ': 2012, 'fuel': 2015, 'tgt': 2012}
    targetDollarYear = 2012
    cpiValues = {2015: 237, 2014: 236.7, 2013: 233, 2012: 229.6, 2011: 224.9, 2010: 218.1,
                 2009: 214.5, 2008: 215.3, 2007: 207.3, 2006: 201.6, 2005: 195.3}
    return doConversion(paramName, cost, paramDollarYears, targetDollarYear, cpiValues)


# Convert dollar year
def doConversion(paramName, cost, paramDollarYears, targetDollarYear, cpiValues):
    paramDollarYear = paramDollarYears[paramName]
    (cpiTgtYear, cpiParamYear) = (cpiValues[targetDollarYear], cpiValues[paramDollarYear])
    return cost * cpiTgtYear / cpiParamYear


################

################ ROTATE 2D LIST FROM HORIZ TO VERT OR VICE VERSA
def rotate(list2d):

    if list2d == []:
        list2drotated = []
    elif list2d is not None:
        list2drotated = list()
        for col in range(len(list2d[0])):
            list2drotated.append([row[col] for row in list2d])
    else:
        list2drotated = None

    return list2drotated
################


def str_elapsedtime(start_time):

    delta = time.time() - start_time
    m = int(delta / 60)
    s = int(delta % 60)
    str_out = '{0:d}m:{1:02d}s'.format(m, s)

    return str_out
