#Michael Craig
#October 4, 2016


def getDemandAndREGenForUC(day,daysOpt,daysLA,demandScaled,hourlyWindGen,hourlySolarGen):
    """Isolate hourly demand and wind & solar gen for hours input to UC model

    :param day: current UC day
    :param daysOpt: number of days to optimize for
    :param daysLA: number of days LA
    :param demandScaled: hourly annual demand
    :param hourlyWindGen: hourly wind  gen
    :param hourlySolarGen: hourly solar gen
    :return: hourly demand & wind & solar gen for UC hours, UC hours
    """

    hoursForUC = getUCHours(day,daysOpt, daysLA)

    if (day + daysOpt + daysLA) <= 366: #necessary data doesn't extend beyond end of year
        demandUC = [demandScaled[hr-1] for hr in hoursForUC] #-1 b/c hours in year start @ 1, not 0 like Python idx
        hourlyWindGenUC = [hourlyWindGen[hr-1] for hr in hoursForUC] #-1 b/c hours in year start @ 1, not 0 like Python idx
        hourlySolarGenUC = [hourlySolarGen[hr-1] for hr in hoursForUC] #-1 b/c hours in year start @ 1, not 0 like Python idx

    else: #necessary data extends beyond end of year, so copy data from last day
        demandUC = getDataPastEndOfYear(day,daysOpt,daysLA,hoursForUC,demandScaled)
        hourlyWindGenUC = getDataPastEndOfYear(day,daysOpt,daysLA,hoursForUC,hourlyWindGen)
        hourlySolarGenUC = getDataPastEndOfYear(day,daysOpt,daysLA,hoursForUC,hourlySolarGen)

    return (demandUC, hourlyWindGenUC, hourlySolarGenUC, hoursForUC)


def getResForUC(day, daysOpt, daysLA, hourlyRegUp, hourlyRegDown, hourlyFlex, hourlyCont):
    """

    :param day: current UC day
    :param daysOpt: number of days to optimize for
    :param daysLA: number of days LA
    :param hourlyRegUp: hourly annual reg up
    :param hourlyRegDown: hourly annual down req
    :param hourlyFlex:
    :param hourlyCont:
    :return: hourly reg up & down req for UC hours
    """

    hoursForUC = getUCHours(day,daysOpt,daysLA)

    if (day + daysOpt + daysLA) <= 366:
        regUpUC = [hourlyRegUp[hr-1] for hr in hoursForUC] #-1 b/c hours in year start @ 1, not 0 like Python idx
        regDownUC = [hourlyRegDown[hr-1] for hr in hoursForUC] #-1 b/c hours in year start @ 1, not 0 like Python idx
        flexUC = [hourlyFlex[hr-1] for hr in hoursForUC]
        contUC = [hourlyCont[hr-1] for hr in hoursForUC]
    else:
        regUpUC = getDataPastEndOfYear(day,daysOpt,daysLA,hoursForUC,hourlyRegUp)
        regDownUC = getDataPastEndOfYear(day,daysOpt,daysLA,hoursForUC,hourlyRegDown)
        flexUC = getDataPastEndOfYear(day,daysOpt,daysLA,hoursForUC,hourlyFlex)
        contUC = getDataPastEndOfYear(day,daysOpt,daysLA,hoursForUC,hourlyCont)

    return (regUpUC,regDownUC,flexUC,contUC)


def getUCHours(day,daysOpt,daysLA):
    """ Get list with hours of year for UC model

    :param day: first day of UC
    :param daysOpt: num days to optimize for
    :param daysLA: num days LA
    :return: 1d list of hours in UC (1-8760 basis)
    """
    (firstHour, lastHour) = ((day-1)*24+1, ((day-1)+(daysOpt+daysLA))*24)
    hoursForUC = [hr for hr in range(firstHour, int(lastHour)+1)]

    return hoursForUC


def getDataPastEndOfYear(day, daysOpt, daysLA, hoursForUC, dataList):
    """Extend data for # of days past end of year

    :param day:
    :param daysOpt:
    :param daysLA:
    :param hoursForUC:
    :param dataList:
    :return:
    """
    daysExtend = (day + daysOpt + daysLA) - 365
    daysWithData = day - 364
    hoursWithData = hoursForUC[:daysWithData*24]

    return [dataList[hr-1] for hr in hoursWithData] * daysExtend