from GAMSAuxFuncs import *


def getHourlyCapacitiesForDays(fleetUC, hourlyCapacsAllGens, hoursForUC):
    """GET HOURLY CAPACITY VALUES FOR PARTICULAR UC DAY

    Isolates hourly capacities for generators just to hours included in UC.
    Returns dict of gen:capacs

    :param fleetUC:
    :param hourlyCapacsAllGens:
    :param hoursForUC:
    :return:
    """
    (hourlyCapacsUC, hourlyCapacsUCList) = (dict(), [])

    # In order to treat cases where the hours of the year extend beyond 8760
    # repeat hours of last day of the year
    hoursForUCAux = [((8760 - 24 + (hr % 24 if hr % 24 > 0 else 24)) if hr > 8760 else hr) for hr in hoursForUC]

    for row in fleetUC[1:]:
        genSymbol = createGenSymbol(row, fleetUC[0])

        hourlyCapacsDay = [hourlyCapacsAllGens[genSymbol][hr-1] for hr in hoursForUCAux]

        hourlyCapacsUC[genSymbol] = hourlyCapacsDay
        hourlyCapacsUCList.append([genSymbol] + hourlyCapacsDay)

    # write2dListToCSV(hourlyCapacsUCList,'hourlyCapacsUCDay' + str(hoursForUC[0]//24+1) + '.csv')
    return hourlyCapacsUC
