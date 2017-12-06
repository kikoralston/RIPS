import numpy as np


def createTempComponents(temp, temp_breaks=[-10, 0, 10, 20, 30]):
    """Computes temperature components for piecewise linear model

    Keyword arguments:
    temp -- vector of temperatures (degrees celsius)
    temp_breaks -- temperature break points of piecewise linear function
    
    Returns:
    matrix with components in columns
    """

    temp = np.array(temp)
    temp_breaks = np.array(temp_breaks)

    nr = np.shape(temp)[0]
    nb = np.shape(temp_breaks)[0]  # number of bounds
    nc = nb + 1  # number of components = number bounds + 1

    temp_comp = np.zeros(shape=(nr, nc))

    temp_comp[(temp < temp_breaks[0]), 0] = temp[temp < temp_breaks[0]]
    temp_comp[temp >= temp_breaks[0], 0] = temp_breaks[0]

    for i in np.arange(1, nb):
        idx_rows = (temp < temp_breaks[i]) & (temp >= temp_breaks[i - 1])
        temp_comp[idx_rows, i] = temp[idx_rows] - temp_breaks[i - 1]
        idx_rows = (temp >= temp_breaks[i])
        temp_comp[idx_rows, i] = temp_breaks[i] - temp_breaks[i - 1]

    idx_rows = (temp > temp_breaks[nb - 1])
    temp_comp[idx_rows, (nc - 1)] = temp[idx_rows] - temp_breaks[nb - 1]

    return temp_comp


def convertDewPoint2RelHum(dew_point, temp):
    """Converts dew point to relative humidity
    uses inverse of formula (8) in:
    http://journals.ametsoc.org/doi/pdf/10.1175/BAMS-86-2-225
    
    Args:
    dew_point: dew point value in Celsius. Can be a vector or an atomic
    temp: air temperature in Celsius. Can be vector or atomic
    
    OBS: if both arguments are vectors, they must be the same length 
    
    Returns:
    relative humidity in %
    """

    # convert to numpy array (in case it is a list) and convert number to vector
    dew_point = np.array([dew_point]).ravel()
    temp = np.array([temp]).ravel()

    if len(dew_point) != len(temp):
        if (len(dew_point) > 1) & (len(temp) > 1):
            print("\nArguments must be:\n"
                  "    1) Both vectors of same length; or\n"
                  "    2) One vector and one atomic value.\n"
                  "Please check arguments!")
            return None

    # definition of constants
    A1 = 17.625  # dimensionless
    B1 = 243.04  # degrees celsius
    C1 = 610.94  # Pascal

    exp_term = A1 * dew_point / (B1 + dew_point) - A1 * temp / (B1 + temp)

    # element-wise minimum
    rh = np.minimum(np.exp(exp_term), 1)

    return 100 * rh


def convertRelHum2DewPoint(rh, temp):
    """Converts relative humidity to dew point
    uses inverse of formula (8) in:
    http://journals.ametsoc.org/doi/pdf/10.1175/BAMS-86-2-225
    
    Args:
    rh: relative humidity value in % (0 to 100). Can be a vector or an atomic
    temp: air temperature in Celsius. Can be vector or atomic
    
    OBS: if both arguments are vectors, they must be the same length
    
    Returns:
    dew point in Celsius
    """

    # convert to numpy array (in case it is a list) and convert number to vector
    rh = np.array([rh]).ravel()
    temp = np.array([temp]).ravel()

    if len(rh) != len(temp):
        if (len(rh) > 1) & (len(temp) > 1):
            print("\nArguments must be:\n"
                  "    1) Both vectors of same length; or\n"
                  "    2) One vector and one atomic value.\n"
                  "Please check arguments!")
            return None

    # definition of constants
    A1 = 17.625  # dimensionless
    B1 = 243.04  # degress celsius
    C1 = 610.94  # Pascal

    dp = B1 * (np.log(rh / 100.) + A1 * temp / (B1 + temp)) / (
        A1 - np.log(rh / 100.) - A1 * temp / (B1 + temp))

    return dp
