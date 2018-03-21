# Michael Craig
# October 4, 2016
# Functions that process demand data for CE model - select which days are
# included in CE, calculate seasonal weights, get peak demand hour, and calculate
# planning reserve margin.

import copy, operator, os
from AuxFuncs import *


########### SELECT WEEKS FOR EXPANSION #########################################
# Inpzuts: demand for current CE run (1d list w/out head), net demand for current CE run (1d list w/out head),
# hourly wind and solar gen (1d lists w/out heads), num representative days per season,
# whether to include days to represent thermal curtailments, hourly curtailments
# of each generator in current CE year (dictionary mapping generator to 2d list w/ headers
# of hourly curtailments), current CE year
# Outputs: hourly demand, wind and solar values for CE (1d lists w/out headers),
# and hour numbers for whole CE, representative per season, special days, and
# all other season hours (1d lists, all 1-8760 basis).
def selectWeeksForExpansion(zonalDemandProfile, zonalNetDemand, zonalHourlyWindGen,
                            zonalHourlySolarGen, daysPerSeason, selectCurtailDays, hrlyCurtailmentsAllGensInTgtYr,
                            currYear, resultsDir, planningMargin):
    netDemand, demand = sumZonalData(zonalNetDemand), sumZonalData(zonalDemandProfile)
    # Get hour of peak demand in each zone
    peakDemandHourZonal, planningMarginZonal = getPeakDemandHourAndPlanningMarginCEZonal(zonalDemandProfile,
                                                                                         planningMargin)
    # Get hours for special days
    specialDayHours = []
    (peakNetDemandDayHours, netDemandMinusPeak) = getPeakNetDemandDayHours(netDemand)  # 1-8760 basis
    specialDayHours.extend(peakNetDemandDayHours)

    if selectCurtailDays:  # if considering thermal curtailments in selection of special days
        totalHrlyCurtailments = getTotalSystemCurtailments(hrlyCurtailmentsAllGensInTgtYr)
        eliminateLeapYearDay(totalHrlyCurtailments, currYear)
        write2dListToCSV([totalHrlyCurtailments], os.path.join(resultsDir,
                                                               'curtailmentsHourlyCombined' + str(currYear) + '.csv'))
        peakCurtailmentDayHours = getPeakCurtailmentDayHours(totalHrlyCurtailments)  # 1-8760
        if peakCurtailmentDayHours != peakNetDemandDayHours:  # use netDemandMinusPeak if change this
            specialDayHours.extend(peakCurtailmentDayHours)
        peakNetDemandAndCurtailmentDayHours = getPeakDemandPlusCurtailmentDayHours(netDemand,
                                                                                   totalHrlyCurtailments)  # 1-8760
        if (peakNetDemandAndCurtailmentDayHours != peakCurtailmentDayHours and
                    peakNetDemandAndCurtailmentDayHours != peakNetDemandDayHours):  # use netDemandMinusPeak if change this
            specialDayHours.extend(peakNetDemandAndCurtailmentDayHours)
        print('Peak net demand hours:', peakNetDemandDayHours)
        print('Peak curtailment hours:', peakCurtailmentDayHours)
        print('Peak net demand + curtailment hours:', peakNetDemandAndCurtailmentDayHours)
        print('Special hours:', specialDayHours)
    # Get representative hours by NLDC
    (repSeasonalHours, repHrsBySeason, regHrsBySeason) = getRepSeasonalHoursByNLDC(netDemand,
                                                                                   daysPerSeason, specialDayHours)
    # Create combined dictionary of rep hrs by season and 'special':special hours
    # (used for getting max hydro generation in each time block)
    repAndSpeHoursDict = copy.deepcopy(repHrsBySeason)
    repAndSpeHoursDict['special'] = specialDayHours
    # Combine representative w/ special and peak hours
    hoursForCE = copy.copy(specialDayHours) + copy.copy(repSeasonalHours)
    for zone in peakDemandHourZonal:
        if peakDemandHourZonal[zone] not in hoursForCE: hoursForCE.append(peakDemandHourZonal[zone])
    demandCE = [demand[hr - 1] for hr in hoursForCE]
    (demandCEZonal, hourlyWindGenCEZonal, hourlySolarGenCEZonal) = isolateDemandAndREGenForCEZonal(hoursForCE,
                                                                                                   zonalDemandProfile,
                                                                                                   zonalHourlyWindGen,
                                                                                                   zonalHourlySolarGen)
    return (demandCE, hoursForCE, repHrsBySeason, specialDayHours, regHrsBySeason, demandCEZonal,
            hourlyWindGenCEZonal, hourlySolarGenCEZonal, peakDemandHourZonal, planningMarginZonal, repAndSpeHoursDict)


##### SUM DATA ACROSS ZONES
def sumZonalData(zonalData):
    totalData = list()
    for zone in zonalData:
        if totalData == []:
            totalData = copy.deepcopy(zonalData[zone])
        else:
            totalData = list(map(operator.add, totalData, zonalData[zone]))
    return totalData


##### SELECT CE HOURS FOR SPECIAL DAYS
# Get hours of entire day w/ peak net demand
# Input: net demand (1d lits w/out head)
# Output: hours of day w/ peak net demand (1d list, 1-8760 basis),
# net demand without peak demand day (1d list)
def getPeakNetDemandDayHours(netDemand):
    peakDemandDayHour = netDemand.index(max(netDemand))
    peakNetDemandDayHours = getHoursOfDayForGivenHour(peakDemandDayHour)  # 1-8760 basis
    netDemandMinusPeak = removeHoursFrom1dList(netDemand, peakNetDemandDayHours)
    return (peakNetDemandDayHours, netDemandMinusPeak)  # 1-8760 basis


# Get all hours for day for a given hour of year.
# Inputs: single hour (0-8759 basis)
# Outputs: hours for entire day (Returns value on 1-8760 basis)
def getHoursOfDayForGivenHour(inputHour):
    return getHoursOfDayForGivenDay(inputHour // 24)


# Get all hours for day for a given day
# Input: day in 0-364 basis. Output: hour in 1-8760 basis.
def getHoursOfDayForGivenDay(numDay):
    hoursPerDay = 24
    (dayStartHour, dayEndHour) = (numDay * hoursPerDay, (numDay + 1) * hoursPerDay)
    return [hr + 1 for hr in range(dayStartHour, dayEndHour)]


# Remove given list of hours from a 1d list, where hours are in 1-8760 format
def removeHoursFrom1dList(list1d, hoursToRemove):
    list1dCopy = copy.copy(list1d)
    for hr in reversed(hoursToRemove): list1dCopy.pop(hr - 1)
    return list1dCopy


# Determine total system hourly curtailments by summing curtailments for each generator
# Input: dict mapping each gen to 2d list of datetime for year of run to hourly net capacity curtailments (MW)
# Output: 1d list w/ total system hourly thermal curtailments
def getTotalSystemCurtailments(hrlyCurtailmentsAllGensInTgtYr):
    totalHrlyCurtailments = None
    for gen in hrlyCurtailmentsAllGensInTgtYr:
        genHrlyCurtailments = hrlyCurtailmentsAllGensInTgtYr[gen]
        if totalHrlyCurtailments == None:
            totalHrlyCurtailments = copy.deepcopy(genHrlyCurtailments)
        else:
            totalHrlyCurtailments = list(map(operator.add, totalHrlyCurtailments, genHrlyCurtailments))
    return totalHrlyCurtailments


# Remove extra leap year day from hourly curtailment data (RBM includes leap year)
# Inputs: total system hourly curtailments (1d list w/out head), curr year
def eliminateLeapYearDay(totalHrlyCurtailments, currYear):
    leapYears = [yr for yr in range(2016, 2101, 4)]
    if currYear in leapYears: totalHrlyCurtailments = totalHrlyCurtailments[:-24]


# Get hours for entire day for hour w/ peak system curtailment
# Inputs: total system curtailment (1d list w/out head)
# Output: 1d list of hours of day w/ peak hourly curtailment (output on 1-8760 basis)
def getPeakCurtailmentDayHours(totalHrlyCurtailments):
    peakCurtailmentDay = totalHrlyCurtailments.index(max(totalHrlyCurtailments))
    return getHoursOfDayForGivenHour(peakCurtailmentDay)  # output in 1-8760 basis


# Get hours for entire day for day w/ peak net demand + hourly curtailment
# Inputs: 1d lists w/out headers of net demand & system curtailments
# Output: 1d list of hours of day (1-8760 basis)
def getPeakDemandPlusCurtailmentDayHours(netDemand, totalHrlyCurtailments):
    (maxPeakDemandAndCurtailment, maxPeakDemandAndCurtailmentDay) = (None, None)
    assert (len(netDemand) == len(totalHrlyCurtailments))
    demandPlusCurtail = list(map(operator.add, netDemand, totalHrlyCurtailments))
    return getHoursOfDayForGivenHour(demandPlusCurtail.index(max(demandPlusCurtail)))


##### SELECT CE HOURS FOR REPRESENTATIVE DAYS PER SEASON
# Inputs: net demand (1d list), num representative days per season to select,
# set of hours already incluced as special days in CE (1d list, 1-8760 basis)
# Outputs: rep hours for each season (1d list no head, 1-8760 basis), dictionaries mapping
# seasons to representative and regular hours
def getRepSeasonalHoursByNLDC(netDemand, daysPerSeason, specialDayHours):
    seasons = ['winter', 'spring', 'summer', 'fall']
    seasonMonths = {'winter': [1, 2, 12], 'spring': [3, 4, 5], 'summer': [6, 7, 8], 'fall': [9, 10, 11]}
    (repHrsBySeason, regHrsBySeason, allRepSeasonalHours) = (dict(), dict(), [])
    for season in seasons:
        monthsInSeason = seasonMonths[season]
        (seasonRepHours, otherSeasonHours) = getSeasonRepHoursByNLDC(netDemand,
                                                                     daysPerSeason, monthsInSeason,
                                                                     specialDayHours)  # 1-8760
        allRepSeasonalHours.extend(seasonRepHours)
        repHrsBySeason[season] = seasonRepHours
        regHrsBySeason[season] = otherSeasonHours
    return (allRepSeasonalHours, repHrsBySeason, regHrsBySeason)


# Get representative hours for given season
# Inputs: net demand (1d list), rep days per season to select, months in season (1d list, 1-8760 basis),
# hours already included in CE model as special days (1d list, 1-8760 basis)
# Outputs: rep hours for given season (1d list), all other hours for given season (1d list)
# (both sets of hours start at 1)
def getSeasonRepHoursByNLDC(netDemand, daysPerSeason, monthsInSeason, specialDayHours):
    hoursInMonths = getHoursInMonths(monthsInSeason)  # starting @ 1
    hoursInMonthsNotSpecial = [hour for hour in hoursInMonths if hour not in specialDayHours]  # starting at 1
    netDemandInMonths = [netDemand[hr - 1] for hr in hoursInMonthsNotSpecial]  # index backwards so hour 1 = idx 0
    print('Months in season:', monthsInSeason)
    seasonRepHours = selectRepresentativeHours(netDemandInMonths, hoursInMonthsNotSpecial,
                                               daysPerSeason)  # starting at 1
    otherSeasonHours = [hr for hr in hoursInMonthsNotSpecial if hr not in seasonRepHours]  # starting at 1
    return (seasonRepHours, otherSeasonHours)  # starting at 1


# Get 1d list of hours (starting at 1) in given list of months
# Input: 1d list of months
# Output: hours in given months (1d list, hours start at 1 in year)
def getHoursInMonths(months):
    daysPerMonth = [(1, 31), (2, 28), (3, 31), (4, 30), (5, 31), (6, 30), (7, 31), (8, 31), (9, 30), (10, 31), (11, 30),
                    (12, 31)]
    firstDayInMonthsAsDayInYear = getFirstDayInMonthsAsDayInYear(daysPerMonth)
    daysInMonths = []
    for month in months:
        firstDayInMonth = firstDayInMonthsAsDayInYear[month - 1]
        daysInMonth = [day for day in range(firstDayInMonth, firstDayInMonth + daysPerMonth[month - 1][1])]
        daysInMonths.extend(daysInMonth)
    hoursInMonths = []
    for day in daysInMonths: hoursInMonths.extend([val + 1 for val in range(day * 24, (day + 1) * 24)])
    return hoursInMonths  # starts @ 1


# Get 1d list of first day each month as day in year (starting at 0)
# Inputs: num days each month (list of tuples)
# Outputs: first day in each month (1d list) (days start at 0)
def getFirstDayInMonthsAsDayInYear(daysPerMonth):
    firstDayInMonthsAsDayInYear = []
    for idx in range(len(daysPerMonth)):
        if idx == 0:
            firstDayInMonthsAsDayInYear.append(0)
        else:
            firstDayInMonthsAsDayInYear.append(firstDayInMonthsAsDayInYear[idx - 1] + daysPriorMonth)
        (lastMonth, daysPriorMonth) = (daysPerMonth[idx][0], daysPerMonth[idx][1])
    return firstDayInMonthsAsDayInYear  # starts at 0


# Select representiative hours for months in a season
# Inputs: net demand in months (1d list), hours in months not already included
# as special hours (1d list, hours start @ 1), num rep days per season to select
# Outputs: rep hours for season (1d list) (hours start @ 1)
def selectRepresentativeHours(netDemandInMonths, hoursInMonthsNotSpecial, daysPerSeason):
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
    print('Lowest NRMSE:', lowestRmse / (max(netDemandInMonths) - min(netDemandInMonths)))
    return hoursLowestRmse


# Calculate RMSE b/wn 2 sets of data
def getRMSE(sampleData, originalData):
    sampleNLDC = sorted(sampleData)
    originalNLDC = sorted(originalData)
    squaredErrors = [(sampleNLDC[idx] - originalNLDC[idx]) ** 2 for idx in range(len(sampleData))]
    rmse = (sum(squaredErrors) / len(squaredErrors)) ** 0.5
    return rmse


##### ISOLATE DEMAND AND RE GENERATION FOR HOURS FOR CE
# Inputs: hours for CE model (1-8760 basis, 1d list), hourly demand, wind and solar
# for whole year for curr CE year (1d lists)
# Outputs: hourly demand, wind & solar for next CE run (1d lists)
# def isolateDemandAndREGenForCE(hoursForCE,demandScaled,hourlyWindGen,hourlySolarGen):
#     demandCE = [demandScaled[hr-1] for hr in hoursForCE] #-1 b/c hours in year start @ 1, not 0 like Python idx
#     hourlyWindGenCE = [hourlyWindGen[hr-1] for hr in hoursForCE] #-1 b/c hours in year start @ 1, not 0 like Python idx
#     hourlySolarGenCE = [hourlySolarGen[hr-1] for hr in hoursForCE] #-1 b/c hours in year start @ 1, not 0 like Python idx
#     return (demandCE,hourlyWindGenCE,hourlySolarGenCE) 


##### ISOLATE DEMAND AND RE GENERATION FOR HOURS FOR CE BY ZONE
# Inputs: hours for CE model (1-8760 basis, 1d list), dict of zone:hourly demand/wind/solar
# for whole year for curr CE year (1d lists)
# Outputs: dict of zones to hourly demand, wind or solar for next CE run
def isolateDemandAndREGenForCEZonal(hoursForCE, zonalDemandProfile, zonalHourlyWindGen, zonalHourlySolarGen):
    demandCEZonal, hourlyWindGenCEZonal, hourlySolarGenCEZonal = dict(), dict(), dict()
    for zone in zonalDemandProfile:
        demandCEZonal[zone] = [zonalDemandProfile[zone][hr - 1] for hr in
                               hoursForCE]  # -1 b/c hours in year start @ 1, not 0 like Python idx
        hourlyWindGenCEZonal[zone] = [zonalHourlyWindGen[zone][hr - 1] for hr in
                                      hoursForCE]  # -1 b/c hours in year start @ 1, not 0 like Python idx
        hourlySolarGenCEZonal[zone] = [zonalHourlySolarGen[zone][hr - 1] for hr in
                                       hoursForCE]  # -1 b/c hours in year start @ 1, not 0 like Python idx
    return (demandCEZonal, hourlyWindGenCEZonal, hourlySolarGenCEZonal)


########### CALCULATE SEASONAL WEIGHTS TO SCALE REP. DEMAND TO SEASON VALUE ####
# Inputs: hourly demand in curr CE year (1d list w/out headers), 1d list of
# representative hours per season (1-8760 basis), 1d list of regular (i.e. non-rep)
# hours per season (1-8760 basis)
# Outputs: map of season to weight to scale rep demand to full season demand (scalar)
def calculateSeasonalWeights(zonalDemandProfile, repHrsBySeason, regHrsBySeason):
    seasonDemandWeights, weightsList = dict(), [['Season', 'SeasonWeight']]
    demand = sumZonalData(zonalDemandProfile)
    for season in repHrsBySeason:
        (repHrs, regHrs) = (repHrsBySeason[season], regHrsBySeason[season])
        (repHourlyDemand, regHourlyDemand) = ([demand[hr - 1] for hr in repHrs],
                                              [demand[hr - 1] for hr in regHrs])
        seasonWeight = (sum(regHourlyDemand) + sum(repHourlyDemand)) / sum(repHourlyDemand)
        seasonDemandWeights[season] = seasonWeight
        weightsList.append([season, seasonWeight])
    return seasonDemandWeights, weightsList


########### GET PEAK DEMAND HOUR AND PLANNING RESERVE CAPACITY #################
# Inputs: dict of zonal to hourly demand, hours in CE, and margin as fraction
# Outputs: dict of zone: hour of peak demand in that zone
def getPeakDemandHourAndPlanningMarginCEZonal(demandCEZonal, planningMargin):
    peakDemandHourZonal, planningMarginZonal = dict(), dict()
    for zone in demandCEZonal:
        zoneDemand = demandCEZonal[zone]
        peakDemandHourZonal[zone] = zoneDemand.index(max(zoneDemand)) + 1  # +1 b/c idx of pos 0 = 0, want 1
        planningMarginZonal[zone] = max(zoneDemand) * (1 + planningMargin)
    return peakDemandHourZonal, planningMarginZonal


########### TEST FUNCTIONS #####################################################
def testGetRMSE():
    print('testing getRMSE')
    assert (almostEqual(getRMSE([1, 1, 1], [1, 1, 1]), 0))
    assert (almostEqual(getRMSE([1, 2, 10], [1, 2, 10]), 0))
    assert (almostEqual(getRMSE([1, 2, 0], [1, 2, 4]), ((1 + 1 + 4) / 3) ** 0.5))
    assert (almostEqual(getRMSE([1, 1, 10], [5, 5, 5]), ((16 + 16 + 25) / 3) ** 0.5))


def almostEqual(num1, num2):
    return abs(num2 - num1) < 0.005
