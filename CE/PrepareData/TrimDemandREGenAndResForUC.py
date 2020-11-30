#Michael Craig
#October 4, 2016


def getDemandAndREGenForUC(day, daysOpt, daysLA, demandScaled, hourlyWindGen, hourlySolarGen):
    """Isolate hourly demand and wind & solar gen for hours input to UC model

    :param day: (int) current UC day
    :param daysOpt: (int) number of days to optimize for
    :param daysLA: (int) number of days look ahead
    :param demandScaled: (list) hourly annual demand
    :param hourlyWindGen: (list) hourly wind  gen
    :param hourlySolarGen: (list) hourly solar gen
    :return: lists with hourly demand & wind & solar gen for UC hours, UC hours
    """

    hoursForUC = getUCHours(day, daysOpt, daysLA)

    demandUC, hourlyWindGenUC, hourlySolarGenUC = dict(), dict(), dict()

    if (day + daysOpt + daysLA) <= 366:
        #necessary data doesn't extend beyond end of year

        for zone in demandScaled.keys():
            # -1 b/c hours in year start @ 1, not 0 like Python idx
            demandUC[zone] = [demandScaled[zone][hr-1] for hr in hoursForUC]
            hourlyWindGenUC[zone] = [hourlyWindGen[zone][hr-1] for hr in hoursForUC]
            hourlySolarGenUC[zone] = [hourlySolarGen[zone][hr-1] for hr in hoursForUC]

    else:
        #necessary data extends beyond end of year, so copy data from last day

        for zone in demandScaled.keys():
            demandUC[zone] = getDataPastEndOfYear(day, daysOpt, daysLA, hoursForUC, demandScaled[zone])
            hourlyWindGenUC[zone] = getDataPastEndOfYear(day, daysOpt, daysLA, hoursForUC, hourlyWindGen[zone])
            hourlySolarGenUC[zone] = getDataPastEndOfYear(day, daysOpt, daysLA, hoursForUC, hourlySolarGen[zone])

    return demandUC, hourlyWindGenUC, hourlySolarGenUC, hoursForUC


def getResForUC(day, daysOpt, daysLA, hourlyRegUp, hourlyRegDown, hourlyFlex, hourlyCont):
    """Get reserve requirements for current UC run

    :param day: (int) current UC day
    :param daysOpt: (int) number of days to optimize for
    :param daysLA: (int) number of days look ahead
    :param hourlyRegUp: (list) hourly regulated up reserve requirement
    :param hourlyRegDown: (list) hourly regulated down reserve requirement
    :param hourlyFlex: (list) hourly flexibility reserve requirement
    :param hourlyCont: (list) hourly contigency reserve requirement
    :return: lists with hourly reserve requirements (regulated up & down, flexibility, contigency) for UC hours
    """

    hoursForUC = getUCHours(day, daysOpt, daysLA)

    regUpUC, regDownUC, flexUC, contUC = dict(), dict(), dict(), dict()

    if (day + daysOpt + daysLA) <= 366:

        for zone in hourlyRegUp.keys():
            # -1 b/c hours in year start @ 1, not 0 like Python idx
            regUpUC[zone] = [hourlyRegUp[zone][hr-1] for hr in hoursForUC]
            regDownUC[zone] = [hourlyRegDown[zone][hr-1] for hr in hoursForUC]
            flexUC[zone] = [hourlyFlex[zone][hr-1] for hr in hoursForUC]
            contUC[zone] = [hourlyCont[zone][hr-1] for hr in hoursForUC]
    else:
        for zone in hourlyRegUp.keys():
            regUpUC[zone] = getDataPastEndOfYear(day, daysOpt, daysLA, hoursForUC, hourlyRegUp[zone])
            regDownUC[zone] = getDataPastEndOfYear(day, daysOpt, daysLA, hoursForUC, hourlyRegDown[zone])
            flexUC[zone] = getDataPastEndOfYear(day, daysOpt, daysLA, hoursForUC, hourlyFlex[zone])
            contUC[zone] = getDataPastEndOfYear(day, daysOpt, daysLA, hoursForUC, hourlyCont[zone])

    return regUpUC, regDownUC, flexUC, contUC


def getUCHours(day, daysOpt, daysLA):
    """ Get list with hours of year for current UC model simulation

    :param day: (int) current day of UC simulation
    :param daysOpt: (int) num days to optimize for
    :param daysLA: (int) num days look ahead
    :return: 1d list of hours in UC (1-8760 basis)
    """
    (firstHour, lastHour) = ((day-1)*24+1, ((day-1)+(daysOpt+daysLA))*24)
    hoursForUC = [hr for hr in range(firstHour, int(lastHour)+1)]

    return hoursForUC


def getDataPastEndOfYear(day, daysOpt, daysLA, hoursForUC, dataList):
    """Extend data for number of days past end of year

    :param day: (int) current day of UC simulation
    :param daysOpt: (int) num days to optimize for
    :param daysLA: (int) num days look ahead
    :param hoursForUC: (1-d list) hours of the year for current UC simulation
    :param dataList:
    :return:
    """
    daysExtend = (day + daysOpt + daysLA) - 365
    daysWithData = day - 364
    hoursWithData = hoursForUC[:daysWithData*24]

    return [dataList[hr-1] for hr in hoursWithData] * daysExtend