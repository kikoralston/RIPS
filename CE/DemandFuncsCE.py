"""
Michael Craig
October 4, 2016
Functions that process demand data for CE model - select which days are included in CE, calculate seasonal weights,
get peak demand hour, and calculate planning reserve margin.
"""

import copy, operator, os, sys
from collections import OrderedDict
import numpy as np
from AuxFuncs import *
import pprint


def selectWeeksForExpansion(zonalDemandProfile, zonalNetDemand, zonalHourlyWindGen, zonalHourlySolarGen, daysPerSeason,
                            selectCurtailDays, hrlyCurtailmentsAllGensInTgtYr, currYear, resultsDir, planningMargin):
    """SELECT WEEKS FOR EXPANSION

    :param zonalDemandProfile: nested dictionary {gcm: {zone: [1d list hourly demand]}} with demand for current CE run
    :param zonalNetDemand: nested dictionary {gcm: {zone: [1d list hourly net demand]}} with net demand for CE run.
                           (list with Net demand list already removed additional leap year day)
    :param zonalHourlyWindGen: hourly wind gen (1d lists w/out heads)
    :param zonalHourlySolarGen: hourly solar gen (1d lists w/out heads)
    :param daysPerSeason: num representative days per season (general parameter of simulation)
    :param selectCurtailDays: (boolean) whether to include days to represent thermal curtailments (general parameter)
    :param hrlyCurtailmentsAllGensInTgtYr: hourly curtailments of each generator in current CE year
                                          (nested dictionary {gcm: {genId: [1d list with hourly curtailments]}}
                                          mapping generator to 1d list without headers of hourly curtailments)
    :param currYear: (int) current CE year
    :param resultsDir: (string) path to folder of results
    :param planningMargin: (scalar) parameter with planning reserve margin (general parameter)
    :return: hourly demand, wind and solar values for CE (1d lists w/out headers), and hour numbers for whole CE,
    representative per season, special days, and all other season hours (1d lists, all 1-8760 basis).
    """

    netDemand = OrderedDict()
    for gcm in zonalNetDemand.keys():
        netDemand[gcm] = sumZonalData(zonalNetDemand[gcm])

    demand = OrderedDict()
    for gcm in zonalDemandProfile.keys():
        demand[gcm] = sumZonalData(zonalNetDemand[gcm])

    # Get hour of peak demand in each zone (over all gcms)
    peakDemandHourZonal, planningMarginZonal = getPeakDemandHourAndPlanningMarginCEZonal(zonalDemandProfile,
                                                                                         planningMargin)
    # Get hours for special days

    (peakNetDemandDayHours, netDemandMinusPeak) = getPeakNetDemandDayHours(netDemand)  # 1-8760 basis

    pp = pprint.PrettyPrinter(indent=3, compact=True)
    print('Peak net demand hours:')
    pp.pprint(peakNetDemandDayHours)

    specialDayHours = {gcm: copy.deepcopy(peakNetDemandDayHours[gcm]) for gcm in peakNetDemandDayHours}

    if selectCurtailDays:  # if considering thermal curtailments in selection of special days

        totalHrlyCurtailments = getTotalSystemCurtailments(hrlyCurtailmentsAllGensInTgtYr)

        totalHrlyCurtailments = eliminateLeapYearDay(totalHrlyCurtailments, currYear)

        if len(totalHrlyCurtailments) > 0:
            writeDictToCSV(totalHrlyCurtailments,
                           os.path.join(resultsDir,'curtailmentsHourlyCombined'+str(currYear)+'.csv'))

            peakCurtailmentDayHours = getPeakCurtailmentDayHours(totalHrlyCurtailments)  # 1-8760

            for g in peakNetDemandDayHours:
                if peakCurtailmentDayHours[g] != peakNetDemandDayHours[g]:  # use netDemandMinusPeak if change this
                    specialDayHours[g].extend(peakCurtailmentDayHours[g])

            peakNetDemandAndCurtailmentDayHours = getPeakDemandPlusCurtailmentDayHours(netDemand,
                                                                                       totalHrlyCurtailments)  # 1-8760
            for g in peakNetDemandDayHours:
                if (peakNetDemandAndCurtailmentDayHours[g] != peakCurtailmentDayHours[g] and
                        peakNetDemandAndCurtailmentDayHours[g] != peakNetDemandDayHours[g]):  # use netDemandMinusPeak if change this
                    specialDayHours[g].extend(peakNetDemandAndCurtailmentDayHours[g])

            print('Peak curtailment hours:')
            pp.pprint(peakCurtailmentDayHours)
            print('Peak net demand + curtailment hours:')
            pp.pprint(peakNetDemandAndCurtailmentDayHours)

        else:
            print('List of hourly curtailment is empty!')

    print('Special hours:')
    pp.pprint(specialDayHours)

    # Get representative hours by NLDC
    repSeasonalHours, repHrsBySeason, regHrsBySeason = dict(), dict(), dict()
    for gcm in netDemand:
        (repSeasonalHours[gcm], repHrsBySeason[gcm],
         regHrsBySeason[gcm]) = getRepSeasonalHoursByNLDC(netDemand[gcm], daysPerSeason, specialDayHours[gcm])

    # Create combined dictionary of rep hrs by season and 'special':special hours
    # (used for getting max hydro generation in each time block)
    repAndSpeHoursDict = dict()
    for gcm in netDemand:
        repAndSpeHoursDict[gcm] = copy.deepcopy(repHrsBySeason[gcm])
        repAndSpeHoursDict[gcm]['special'] = specialDayHours[gcm]

    # Combine representative w/ special and peak hours
    hoursForCE = dict()
    for gcm in netDemand:
        hoursForCE[gcm] = copy.copy(specialDayHours[gcm]) + copy.copy(repSeasonalHours[gcm])

    for gcm in peakDemandHourZonal.keys():
        for zone in peakDemandHourZonal[gcm].keys():
            if peakDemandHourZonal[gcm][zone] not in hoursForCE[gcm]:
                hoursForCE[gcm].append(peakDemandHourZonal[gcm][zone])

    demandCE = dict()
    for gcm in netDemand:
        demandCE[gcm] = [demand[gcm][hr - 1] for hr in hoursForCE[gcm]]

    (demandCEZonal, hourlyWindGenCEZonal, hourlySolarGenCEZonal) = isolateDemandAndREGenForCEZonal(hoursForCE,
                                                                                                   zonalDemandProfile,
                                                                                                   zonalHourlyWindGen,
                                                                                                   zonalHourlySolarGen)
    return (demandCE, hoursForCE, repHrsBySeason, specialDayHours, regHrsBySeason, demandCEZonal,
            hourlyWindGenCEZonal, hourlySolarGenCEZonal, peakDemandHourZonal, planningMarginZonal, repAndSpeHoursDict)


def sumZonalData(zonalData):
    """SUM DATA ACROSS ZONES

    :param zonalData: dict zone: 1d list
    :return: 1d list with total accumulated data across zones
    """
    totalData = list()
    for zone in zonalData:
        if totalData == []:
            totalData = copy.deepcopy(zonalData[zone])
        else:
            totalData = list(map(operator.add, totalData, zonalData[zone]))
    return totalData


def getPeakNetDemandDayHours(netDemand):
    """SELECT CE HOURS FOR SPECIAL DAYS

    Get hours of entire day w/ peak net demand

    :param netDemand: dictionary with {gcm: [net demand]}
    :return: hours of day w/ peak net demand (dictionary with {gcm: [hours of days]}, 1-8760 basis),
             net demand without peak demand day (dictionary with {gcm: [net demand]})
    """

    peakNetDemandDayHours = {}
    netDemandMinusPeak = {}

    for gcm in netDemand.keys():
        peakDemandDayHour = netDemand[gcm].index(max(netDemand[gcm]))

        peakNetDemandDayHours[gcm] = getHoursOfDayForGivenHour(peakDemandDayHour)  # 1-8760 basis
        netDemandMinusPeak[gcm] = removeHoursFrom1dList(netDemand[gcm], peakNetDemandDayHours[gcm])

    return peakNetDemandDayHours, netDemandMinusPeak  # 1-8760 basis


def getHoursOfDayForGivenHour(inputHour):
    """Get all hours for day for a given hour of year.

    :param inputHour: single hour (0-8759 basis)
    :return: hours for entire day (Returns value on 1-8760 basis)
    """
    return getHoursOfDayForGivenDay(inputHour // 24)


def getHoursOfDayForGivenDay(numDay):
    """Get all hours for day for a given day

    :param numDay: day in 0-364 basis
    :return: hour in 1-8760 basis.
    """
    hoursPerDay = 24
    (dayStartHour, dayEndHour) = (numDay * hoursPerDay, (numDay + 1) * hoursPerDay)
    return [hr + 1 for hr in range(dayStartHour, dayEndHour)]


def removeHoursFrom1dList(list1d, hoursToRemove):
    """Remove given list of hours from a 1d list, where hours are in 1-8760 format

    :param list1d: 1d list (general)
    :param hoursToRemove: list with hours to be removed
    :return: updated id list
    """
    list1dCopy = copy.copy(list1d)
    for hr in reversed(hoursToRemove): list1dCopy.pop(hr - 1)
    return list1dCopy


def getTotalSystemCurtailments(hrlyCurtailmentsAllGensInTgtYr):
    """Determine total system hourly curtailments by summing curtailments for each generator

    :param hrlyCurtailmentsAllGensInTgtYr: nested dict mapping in each gcm each gen to 2d list of datetime for year of
                                           run to hourly net capacity curtailments (MW)
                                           {gcm: genId: [1d list hourly net capacity curtailments]}
    :return: dict {gcm: [1d list w/ total system hourly thermal curtailments]}
    """
    totalHrlyCurtailments = {}
    for gcm in hrlyCurtailmentsAllGensInTgtYr:
        auxHrlyCurtailments = []
        for gen in hrlyCurtailmentsAllGensInTgtYr[gcm]:
            genHrlyCurtailments = hrlyCurtailmentsAllGensInTgtYr[gcm][gen]

            if len(auxHrlyCurtailments) == 0:
                auxHrlyCurtailments = list(copy.deepcopy(genHrlyCurtailments))
            else:
                auxHrlyCurtailments = list(map(operator.add, auxHrlyCurtailments, genHrlyCurtailments))

        totalHrlyCurtailments[gcm] = auxHrlyCurtailments

    return totalHrlyCurtailments


def eliminateLeapYearDay(totalHrlyCurtailments, currYear):
    """Remove extra leap year day from hourly curtailment data (RBM includes leap year)

    :param totalHrlyCurtailments: total system hourly curtailments in each gcm (dict {gcm: [1d list w/out head]})
    :param currYear: current year
    :return: (dict {gcm: [1d list w/out head]})
    """
    leapYears = [yr for yr in range(2016, 2101, 4)]

    if currYear in leapYears:
        out_dict = OrderedDict()

        for gcm in totalHrlyCurtailments:
            out_dict[gcm] = totalHrlyCurtailments[gcm][:-24]
    else:
        out_dict = totalHrlyCurtailments

    return out_dict


def getPeakCurtailmentDayHours(totalHrlyCurtailments):
    """Get hours for entire day for hour w/ peak system curtailment

    :param totalHrlyCurtailments: total system curtailment dict {gcm: [1d list w/out head]}
    :return: dict {gcm: [1d list w/out head]} of hours of day w/ peak hourly curtailment in each gcm (1-8760 basis)
    """

    peakCurtailmentHour = {g: totalHrlyCurtailments[g].index(min(totalHrlyCurtailments[g]))
                           for g in totalHrlyCurtailments}

    # output in 1-8760 basis
    peakCurtailment24Hours = {g: getHoursOfDayForGivenHour(peakCurtailmentHour[g]) for g in peakCurtailmentHour}

    return peakCurtailment24Hours


def getPeakDemandPlusCurtailmentDayHours(netDemand, totalHrlyCurtailments):
    """Get hours for entire day for day w/ peak net demand + hourly curtailment

    :param netDemand: net demand dict {gcm: [1d list w/out head]}
    :param totalHrlyCurtailments: total system curtailment dict {gcm: [1d list w/out head]}
    :return: dict with 1d list of hours of day (1-8760 basis) for each gcm {gcm: [1d list w/out head]}
    """
    (maxPeakDemandAndCurtailment, maxPeakDemandAndCurtailmentDay) = (None, None)
    assert (netDemand.keys() == totalHrlyCurtailments.keys())

    out_dict = dict()
    for g in netDemand:
        assert (len(netDemand[g]) == len(totalHrlyCurtailments[g]))
        demandPlusCurtail = list(map(operator.add, netDemand[g], totalHrlyCurtailments[g]))

        out_dict[g] = getHoursOfDayForGivenHour(demandPlusCurtail.index(max(demandPlusCurtail)))

    return out_dict


def getRepSeasonalHoursByNLDC(netDemand, daysPerSeason, specialDayHours):
    """SELECT CE HOURS FOR REPRESENTATIVE DAYS PER SEASON

    :param netDemand: net demand (1d list),
    :param daysPerSeason: num representative days per season to select,
    :param specialDayHours: set of hours already incluced as special days in CE (1d list, 1-8760 basis)
    :return: rep hours for each season (1d list no head, 1-8760 basis), dictionaries mapping
    seasons to representative and regular hours
    """
    seasons = ['winter', 'spring', 'summer', 'fall']
    seasonMonths = {'winter': [1, 2, 12], 'spring': [3, 4, 5], 'summer': [6, 7, 8], 'fall': [9, 10, 11]}
    repHrsBySeason, regHrsBySeason, allRepSeasonalHours = dict(), dict(), []

    for season in seasons:
        monthsInSeason = seasonMonths[season]
        (seasonRepHours, otherSeasonHours) = getSeasonRepHoursByNLDC(netDemand,
                                                                     daysPerSeason, monthsInSeason,
                                                                     specialDayHours)  # 1-8760
        allRepSeasonalHours.extend(seasonRepHours)
        repHrsBySeason[season] = seasonRepHours
        regHrsBySeason[season] = otherSeasonHours
    return allRepSeasonalHours, repHrsBySeason, regHrsBySeason


def getSeasonRepHoursByNLDC(netDemand, daysPerSeason, monthsInSeason, specialDayHours):
    """Get representative hours for given season

    :param netDemand: net demand (1d list),
    :param daysPerSeason: rep days per season to select,
    :param monthsInSeason: months in season (1d list, 1-8760 basis)
    :param specialDayHours: hours already included in CE model as special days (1d list, 1-8760 basis)
    :return: rep hours for given season (1d list), all other hours for given season (1d list)
    (both sets of hours start at 1)
    """
    hoursInMonths = getHoursInMonths(monthsInSeason)  # starting @ 1
    hoursInMonthsNotSpecial = [hour for hour in hoursInMonths if hour not in specialDayHours]  # starting at 1
    netDemandInMonths = [netDemand[hr - 1] for hr in hoursInMonthsNotSpecial]  # index backwards so hour 1 = idx 0
    print('Months in season:', monthsInSeason)
    seasonRepHours = selectRepresentativeHours(netDemandInMonths, hoursInMonthsNotSpecial,
                                               daysPerSeason)  # starting at 1
    otherSeasonHours = [hr for hr in hoursInMonthsNotSpecial if hr not in seasonRepHours]  # starting at 1
    return (seasonRepHours, otherSeasonHours)  # starting at 1


def getDaysPerMonth():
    """Returns list of tuples with number of days in each month

    :return: list of tuples with number of days in each month (month, number of days)
    """
    daysPerMonth = [(1, 31), (2, 28), (3, 31), (4, 30), (5, 31), (6, 30), (7, 31), (8, 31), (9, 30), (10, 31), (11, 30),
                    (12, 31)]
    return daysPerMonth


def getHoursInMonths(months):
    """Get 1d list of hours (1-8760 basis) in given list of months

    :param months: 1d list of months
    :return: hours in given months (1d list, hours start at 1 in year)
    """
    daysPerMonth = getDaysPerMonth()

    firstDayInMonthsAsDayInYear = getFirstDayInMonthsAsDayInYear(daysPerMonth)
    daysInMonths = []
    for month in months:
        firstDayInMonth = firstDayInMonthsAsDayInYear[month - 1]
        daysInMonth = [day for day in range(firstDayInMonth, firstDayInMonth + daysPerMonth[month - 1][1])]
        daysInMonths.extend(daysInMonth)
    hoursInMonths = []
    for day in daysInMonths: hoursInMonths.extend([val + 1 for val in range(day * 24, (day + 1) * 24)])
    return hoursInMonths  # starts @ 1


def getFirstDayInMonthsAsDayInYear(daysPerMonth):
    """Get 1d list of first day each month as day in year (starting at 0)

    :param daysPerMonth: num days each month (list of tuples)
    :return: first day in each month (1d list) (days start at 0)
    """
    firstDayInMonthsAsDayInYear = []
    for idx in range(len(daysPerMonth)):
        if idx == 0:
            firstDayInMonthsAsDayInYear.append(0)
        else:
            firstDayInMonthsAsDayInYear.append(firstDayInMonthsAsDayInYear[idx - 1] + daysPriorMonth)
        (lastMonth, daysPriorMonth) = (daysPerMonth[idx][0], daysPerMonth[idx][1])
    return firstDayInMonthsAsDayInYear  # starts at 0


def selectRepresentativeHours(netDemandInMonths, hoursInMonthsNotSpecial, daysPerSeason):
    """Select representiative hours for months in a season

    :param netDemandInMonths: net demand in months (1d list)
    :param hoursInMonthsNotSpecial: hours in months not already included as special hours (1d list, hours start @ 1)
    :param daysPerSeason: num rep days per season to select
    :return: rep hours for season (1d list) (hours start @ 1)
    """
    hoursPerDay = 24
    hoursLowestRmse, lowestRmse = [], sum(netDemandInMonths) ** 2
    for firstHourInDay in range(0, len(netDemandInMonths) - hoursPerDay * (daysPerSeason - 1), hoursPerDay):
        hrInNetDemandSample = [hr for hr in range(firstHourInDay, firstHourInDay + hoursPerDay * daysPerSeason)]
        actualHrsSample = [hoursInMonthsNotSpecial[hr] for hr in hrInNetDemandSample]
        maxDiffHrs = max([actualHrsSample[idx] - actualHrsSample[idx - 1] for idx in range(1, len(actualHrsSample))])
        if maxDiffHrs == 1:
            netDemandInDay = [netDemandInMonths[hr] for hr in hrInNetDemandSample]
            netDemandInDayForAllMonths = netDemandInDay * (len(netDemandInMonths) // len(netDemandInDay))
            truncatedNetDemandInMonths = netDemandInMonths[:len(netDemandInDayForAllMonths)]
            rmse = getRMSE(netDemandInDayForAllMonths, truncatedNetDemandInMonths)
            if rmse < lowestRmse:
                lowestRmse = rmse
                hoursLowestRmse = copy.copy(actualHrsSample)
    if max(netDemandInMonths) == min(netDemandInMonths):
        print('Lowest NRMSE:', lowestRmse / (max(netDemandInMonths) - min(netDemandInMonths) + 0.1))
    else:
        print('Lowest NRMSE:', lowestRmse / (max(netDemandInMonths) - min(netDemandInMonths)))

    return hoursLowestRmse


def getRMSE(sampleData, originalData):
    """Calculate RMSE b/wn 2 sets of data

    :param sampleData:
    :param originalData:
    :return:
    """
    sampleNLDC = sorted(sampleData)
    originalNLDC = sorted(originalData)
    squaredErrors = [(sampleNLDC[idx] - originalNLDC[idx]) ** 2 for idx in range(len(sampleData))]
    rmse = (sum(squaredErrors) / len(squaredErrors)) ** 0.5
    return rmse


def isolateDemandAndREGenForCEZonal(hoursForCE, zonalDemandProfile, zonalHourlyWindGen, zonalHourlySolarGen):
    """ISOLATE DEMAND AND RE GENERATION FOR HOURS FOR CE BY ZONE

    :param hoursForCE: hours for CE model (1-8760 basis, 1d list), dict of zone:hourly demand/wind/solar
    :param zonalDemandProfile: dict of zone:hourly demand for whole year for curr CE year (1d list)
    :param zonalHourlyWindGen: dict of zone:hourly wind for whole year for curr CE year (1d list)
    :param zonalHourlySolarGen: dict of zone:hourly solar for whole year for curr CE year (1d list)
    :return: dict of zones to hourly demand, wind or solar for next CE run
    """
    demandCEZonal, hourlyWindGenCEZonal, hourlySolarGenCEZonal = dict(), dict(), dict()

    for gcm in hoursForCE:
        auxDemand = dict()
        for zone in zonalDemandProfile[gcm]:
            # -1 b/c hours in year start @ 1, not 0 like Python idx
            auxDemand[zone] = [zonalDemandProfile[gcm][zone][hr - 1] for hr in hoursForCE[gcm]]

        demandCEZonal[gcm] = auxDemand

    for gcm in hoursForCE:
        auxDictSolar = dict()
        auxDictWind = dict()
        for zone in zonalHourlyWindGen:
            # -1 b/c hours in year start @ 1, not 0 like Python idx
            auxDictWind[zone] = [zonalHourlyWindGen[zone][hr - 1] for hr in hoursForCE[gcm]]
            auxDictSolar[zone] = [zonalHourlySolarGen[zone][hr - 1] for hr in hoursForCE[gcm]]

        hourlyWindGenCEZonal[gcm] = auxDictWind
        hourlySolarGenCEZonal[gcm] = auxDictSolar

    return demandCEZonal, hourlyWindGenCEZonal, hourlySolarGenCEZonal


def calculateSeasonalWeights(zonalDemandProfile, repHrsBySeason, regHrsBySeason):
    """CALCULATE SEASONAL WEIGHTS TO SCALE REP. DEMAND TO SEASON VALUE

    :param zonalDemandProfile: hourly demand in curr CE year for each zone in each gcm.
    :param repHrsBySeason: 1d list of representative hours per season (1-8760 basis)
    :param regHrsBySeason: 1d list of regular (i.e. non-rep) hours per season (1-8760 basis)
    :return: map of season to weight to scale rep demand to full season demand (scalar)
    """
    seasonDemandWeights, weightsList = dict(), [['Season', 'SeasonWeight']]

    # get list of seasons to iterate over
    seasonsList = repHrsBySeason[list(repHrsBySeason.keys())[0]].keys()

    for season in seasonsList:
        repHourlyDemand, regHourlyDemand = 0, 0
        for gcm in repHrsBySeason:
            demand = sumZonalData(zonalDemandProfile[gcm])
            (repHrs, regHrs) = (repHrsBySeason[gcm][season], regHrsBySeason[gcm][season])
            repHourlyDemand = repHourlyDemand + sum([demand[hr - 1] for hr in repHrs])
            regHourlyDemand = regHourlyDemand + sum([demand[hr - 1] for hr in regHrs])

        # average over all gcms
        seasonWeight = (regHourlyDemand + repHourlyDemand) / repHourlyDemand

        seasonDemandWeights[season] = seasonWeight
        weightsList.append([season, seasonWeight])

    return seasonDemandWeights, weightsList


def getPeakDemandHourAndPlanningMarginCEZonal(demandCEZonal, planningMargin):
    """GET PEAK DEMAND HOUR AND PLANNING RESERVE CAPACITY

    The output dictionary will only keep the largest peak demand value among all gcms (for each zone).
    For example, if we have 2 gcms and 3 zones and the largest peak demand value for z1 occurs in gcm2, for z2 occurs in
    gcm2, and for z3 occurs at gcm1 the output dictionary will look like this:

    {gcm1: {z3: x}, gcm2:{z2: y, z3: w}}

    That is the list of keys in the outer dictionary may not be extensive over all gcms. And the list of keys
    in the inner dictionary may not contain all zones.

    :param demandCEZonal: nested dict of zone to hourly demand in each CE, hours in CE {gcm:zone:[hourly demand]}
    :param planningMargin: margin as fraction
    :return: nested dict of {gcm: {zone: hour of peak demand}} (see observation above)
             dict of {zone: planningMargin} with planning margin (in MW) for each zone
    """
    peakDemandHourZonal, planningMarginZonal = OrderedDict(), OrderedDict()

    gcms = list(demandCEZonal.keys())
    zones = list(demandCEZonal[gcms[0]].keys())

    for z in zones:
        auxDemand = np.array([demandCEZonal[g][z] for g in gcms])
        maxHoursEachGcm = auxDemand.argmax(axis=1)
        peakDemandValuesEachGcm = np.array([auxDemand[i, maxHoursEachGcm[i]] for i in range(len(gcms))])

        peakDemandValue = np.max(peakDemandValuesEachGcm)
        idx_peak_gcm = np.argmax(peakDemandValuesEachGcm)
        maxHourZone = maxHoursEachGcm[idx_peak_gcm] + 1 # +1 b/c idx of pos 0 = 0, want 1

        if gcms[idx_peak_gcm] in peakDemandHourZonal.keys():
            peakDemandHourZonal[gcms[idx_peak_gcm]][z] = maxHourZone
        else:
            aux = OrderedDict()
            aux[z] = maxHourZone
            peakDemandHourZonal[gcms[idx_peak_gcm]] = aux

        planningMarginZonal[z] = peakDemandValue * (1 + planningMargin)

    return peakDemandHourZonal, planningMarginZonal


def get_list_all_hours_ce(hoursForCE):
    """ compiles a 1d list with all hours considered for CE model across all GCM scenarios

    :param hoursForCE: dict {gcm: [hours CE]}
    :return: 1d list with all hours for CE
    """

    a = [hoursForCE[gcm] for gcm in hoursForCE]

    b = []
    for i in a:
        b = b + copy.deepcopy(i)

    b.sort()

    # convert to set (to remove duplicates) and back to list (to be able to sort)
    d = list(set(b))
    d.sort()

    return d


########### TEST FUNCTIONS #####################################################
def testGetRMSE():
    print('testing getRMSE')
    assert (almostEqual(getRMSE([1, 1, 1], [1, 1, 1]), 0))
    assert (almostEqual(getRMSE([1, 2, 10], [1, 2, 10]), 0))
    assert (almostEqual(getRMSE([1, 2, 0], [1, 2, 4]), ((1 + 1 + 4) / 3) ** 0.5))
    assert (almostEqual(getRMSE([1, 1, 10], [5, 5, 5]), ((16 + 16 + 25) / 3) ** 0.5))


def almostEqual(num1, num2):
    return abs(num2 - num1) < 0.005


##### ISOLATE DEMAND AND RE GENERATION FOR HOURS FOR CE
# Inputs: hours for CE model (1-8760 basis, 1d list), hourly demand, wind and solar
# for whole year for curr CE year (1d lists)
# Outputs: hourly demand, wind & solar for next CE run (1d lists)
# def isolateDemandAndREGenForCE(hoursForCE,demandScaled,hourlyWindGen,hourlySolarGen):
#     demandCE = [demandScaled[hr-1] for hr in hoursForCE] #-1 b/c hours in year start @ 1, not 0 like Python idx
#     hourlyWindGenCE = [hourlyWindGen[hr-1] for hr in hoursForCE] #-1 b/c hours in year start @ 1, not 0 like Python idx
#     hourlySolarGenCE = [hourlySolarGen[hr-1] for hr in hoursForCE] #-1 b/c hours in year start @ 1, not 0 like Python idx
#     return (demandCE,hourlyWindGenCE,hourlySolarGenCE)
