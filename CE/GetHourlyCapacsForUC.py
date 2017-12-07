#MIchael Craig, 17 may 2017
#Isolates hourly capacities for generators just to hours included in UC.
#Returns dict of gen:capacs

########### GET HOURLY CAPACITY VALUES FOR PARTICULAR UC DAY ###################
def getHourlyCapacitiesForDays(fleetUC,hourlyCapacsAllGens,hoursForUC):
    (hourlyCapacsUC,hourlyCapacsUCList) = (dict(),[])
    for row in fleetUC[1:]:
        genSymbol = createGenSymbol(row,fleetUC[0])
        hourlyCapacsDay = [hourlyCapacsAllGens[genSymbol][hr-1] for hr in hoursForUC]
        hourlyCapacsUC[genSymbol] = hourlyCapacsDay
        hourlyCapacsUCList.append([genSymbol] + hourlyCapacsDay)
    # write2dListToCSV(hourlyCapacsUCList,'hourlyCapacsUCDay' + str(hoursForUC[0]//24+1) + '.csv')
    return hourlyCapacsUC
