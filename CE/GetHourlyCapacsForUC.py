from GAMSAuxFuncs import *


def getHourlyCapacitiesForDays(fleetUC, hourlyCapacsCurtailedGens, hoursForUC):
    """GET HOURLY CAPACITY VALUES FOR PARTICULAR UC DAY

    Isolates hourly capacities for generators just to hours included in UC.
    Returns dict of gen:capacs

    :param fleetUC:
    :param hourlyCapacsCurtailedGens:
    :param hoursForUC:
    :return:
    """
    (hourlyCapacsUC, hourlyCapacsUCList) = (dict(),[])
    for row in fleetUC[1:]:
        genSymbol = createGenSymbol(row, fleetUC[0])

        if genSymbol in hourlyCapacsCurtailedGens.keys():
            hourlyCapacsDay = [hourlyCapacsCurtailedGens[genSymbol][hr-1] for hr in hoursForUC]
        else:
            # no simulation of curtailment for this gen (assume capacity = 100%)
            hourlyCapacsDay = [1 for hr in hoursForUC]

        hourlyCapacsUC[genSymbol] = hourlyCapacsDay
        hourlyCapacsUCList.append([genSymbol] + hourlyCapacsDay)

    # write2dListToCSV(hourlyCapacsUCList,'hourlyCapacsUCDay' + str(hoursForUC[0]//24+1) + '.csv')
    return hourlyCapacsUC
