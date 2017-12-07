#Michael Craig
#April 13, 2017

#Assess hourly hydropower geneartion data from TVA

import numpy as np
import matplotlib.pyplot as plt
import os, csv
from AuxFuncs import *
from datetime import datetime

plt.style.use('ggplot')
rc = {'font.family':'Times New Roman','font.size':14,'text.color':'k',
    'axes.labelcolor':'k','xtick.color':'k','ytick.color':'k'}
plt.rcParams.update(**rc)

#Inputs
dataDir = 'C:\\Users\\mtcraig\\Desktop\\EPP Research\\Databases\\TVAHydroGenAndDemand\\FromTVA'
demandDataDir = 'C:\\Users\\mtcraig\\Desktop\\EPP Research\\Databases\\TVAHydroGenAndDemand'
monthList = [7,8,9] 
yearListInput = [2014,2015]
plotGenByHourOfDay = False
plotDailyGen = False
plotDailyGenAndDemand = False
plotWeeklyGenAndDemand = True

#Parameters
figCtr = 0
daysPerYear = 365
hoursPerDay = 24
daysInMonth = {5:31,6:30,7:31,8:31,9:30,10:31,11:30,12:31}


##################### BOXPLOTS OF GENERATION BY HOUR OF DAY #########################
if plotGenByHourOfDay == True:
    #Iterate over given years
    for year in yearListInput:
        figCtr += 1
        yearList = [year]  
        fileName = 'hydroDailyProfileSpring' + str(year)

        #Read data and trim first col (turbine flow is in cubic feet per second)
        hourlyHydroData1 = readCSVto2dList(os.path.join(dataDir,'hourlyTurbineFlow01To17Part1.csv'))
        hourlyHydroData2 = readCSVto2dList(os.path.join(dataDir,'hourlyTurbineFlow01To17Part2.csv'))

        #Merge datasets
        hourlyHydroData2 = [row[1:] for row in hourlyHydroData2]
        hourlyHydroData = [row1 + row2 for row1,row2 in zip(hourlyHydroData1,hourlyHydroData2)]
        dateCol = hourlyHydroData[0].index('Date')
        dam1Col = dateCol+1

        #Add columns w/ datetime info
        hourlyHydroData[0] += ['Year','Month','Day','Hour']
        yCol = hourlyHydroData[0].index('Year')
        mCol,dCol,hCol = yCol+1,yCol+2,yCol+3
        for row in hourlyHydroData[1:]:
            currDate = row[dateCol] #M/D/Y H:MM
            splitSlash = currDate.split('/')
            m,d,rest = int(splitSlash[0]),int(splitSlash[1]),splitSlash[2]
            splitYearAndTime = rest.split(' ')
            y,t = int(splitYearAndTime[0]),splitYearAndTime[1]
            h = int(t.split(':')[0])
            row += [y,m,d,h]

        #Isolate year and month data   
        years = [row[yCol] for row in hourlyHydroData[1:]]
        hydroSlim = [hourlyHydroData[0]] + [row for row in hourlyHydroData[1:] if row[yCol] in yearList]
        months = [row[mCol] for row in hydroSlim[1:]]
        hydroSlim = [hydroSlim[0]] + [row for row in hydroSlim[1:] if row[mCol] in monthList]
        fig = plt.figure(figCtr,figsize=(25,20))
        rows,cols = 5,6
        #Loop through each dam
        ctr = 0
        for col in range(dam1Col,yCol):
            #Put gen for given hour of day in new row of 2d list
            allGen = list()
            for h in range(hoursPerDay): #hour as 0-23
                hourRows = [i for i in range(len(hydroSlim)) if hydroSlim[i][hCol] == h]
                gen = [int(hydroSlim[idx][col]) for idx in hourRows if hydroSlim[idx][col] != '']
                allGen.append(gen)
            #Plot data
            ax = plt.subplot(rows,cols,col)
            bp = ax.boxplot(allGen,whis=[10,90])
            #Increase weight of median line
            for median in bp['medians']: median.set(linewidth=5)
            #Adjust axis labels and titles
            plt.title(hydroSlim[0][col])
            if ctr%cols == 0: plt.ylabel('Turb. Flow (cfs)')
            if col>(cols*(rows-1)): plt.xlabel('Hour of Day')
            else: ax.axes.get_xaxis().set_ticks([])
            ctr += 1

        fig.savefig(fileName+'.png',dpi=500,transparent=True, bbox_inches='tight', pad_inches=0.1)

        # plt.show()

##################### PLOT TOTAL GENERATION BY DAY #################################
if plotDailyGen == True:
    #Iterate over each year and month
    for year in yearListInput:
        for month in monthList:
            figCtr += 1  
            fileName = 'hydroDailyGenY' + str(year) + 'M' + str(month)

            #Read data and trim first col (turbine flow is in cubic feet per second)
            hourlyHydroData1 = readCSVto2dList(os.path.join(dataDir,'hourlyTurbineFlow01To17Part1.csv'))
            hourlyHydroData2 = readCSVto2dList(os.path.join(dataDir,'hourlyTurbineFlow01To17Part2.csv'))

            #Merge datasets
            hourlyHydroData2 = [row[1:] for row in hourlyHydroData2]
            hourlyHydroData = [row1 + row2 for row1,row2 in zip(hourlyHydroData1,hourlyHydroData2)]
            dateCol = hourlyHydroData[0].index('Date')
            dam1Col = dateCol+1

            #Add columns w/ datetime info
            hourlyHydroData[0] += ['Year','Month','Day','Hour']
            yCol = hourlyHydroData[0].index('Year')
            mCol,dCol,hCol = yCol+1,yCol+2,yCol+3
            for row in hourlyHydroData[1:]:
                currDate = row[dateCol] #M/D/Y H:MM
                splitSlash = currDate.split('/')
                m,d,rest = int(splitSlash[0]),int(splitSlash[1]),splitSlash[2]
                splitYearAndTime = rest.split(' ')
                y,t = int(splitYearAndTime[0]),splitYearAndTime[1]
                h = int(t.split(':')[0])
                row += [y,m,d,h]

            #Isolate year and month data   
            years = [row[yCol] for row in hourlyHydroData[1:]]
            hydroSlim = [hourlyHydroData[0]] + [row for row in hourlyHydroData[1:] if row[yCol]==year]
            months = [row[mCol] for row in hydroSlim[1:]]
            hydroSlim = [hydroSlim[0]] + [row for row in hydroSlim[1:] if row[mCol]==month]
            fig = plt.figure(figCtr,figsize=(25,20))
            rows,cols = 5,6
            #Loop through each dam
            ctr = 0
            for col in range(dam1Col,yCol):
                #Sum daily turbine flows
                dailyGen = list()
                days = daysInMonth[month]
                for day in range(1,days+1):
                    dayRows = [i for i in range(len(hydroSlim)) if hydroSlim[i][dCol] == day]
                    gen = sum([int(hydroSlim[idx][col]) for idx in dayRows if hydroSlim[idx][col] != ''])
                    dailyGen.append(gen)
                #Plot data
                ax = plt.subplot(rows,cols,col)
                ax.plot(dailyGen)
                #Adjust axis labels and titles
                plt.title(hydroSlim[0][col])
                if ctr%cols == 0: plt.ylabel('Daily Turb. Flow (cfs)')
                if col>(cols*(rows-1)): plt.xlabel('Day')
                else: ax.axes.get_xaxis().set_ticks([])
                ctr += 1

            fig.savefig(fileName+'.png',dpi=100,transparent=True, bbox_inches='tight', pad_inches=0.1)
        
    plt.show()

##################### PLOT TOTAL GENERATION VERSUS TOTAL DEMAND BY DAY #########
if plotDailyGenAndDemand == True:
    #Iterate over each year and month
    for year in yearListInput:
        for month in monthList:
            figCtr += 1  
            fileName = 'hydroDailyGenVersusDemandY' + str(year) + 'M' + str(month)

            #PROCESS DEMAND DATA
            #Import demand data
            hourlyDemandData = readCSVto2dList(os.path.join(demandDataDir,'tvademand.csv'))

            #Add columns w/ datetime info
            dateCol = hourlyDemandData[0].index('time')
            loadCol = hourlyDemandData[0].index('load')
            hourlyDemandData[0] += ['Year','Month','Day','Hour']
            yColDemand = hourlyDemandData[0].index('Year')
            mColDemand,dColDemand,hColDemand = yColDemand+1,yColDemand+2,yColDemand+3
            for row in hourlyDemandData[1:]:
                currDate = row[dateCol] #M/D/Y H:MM
                splitSlash = currDate.split('/')
                m,d,rest = int(splitSlash[0]),int(splitSlash[1]),splitSlash[2]
                splitYearAndTime = rest.split(' ')
                y,t = int(splitYearAndTime[0]),splitYearAndTime[1]
                h = int(t.split(':')[0])
                row += [y,m,d,h]

            #PROCESS HYDRO GENERATION DATA
            #Read data and trim first col (turbine flow is in cubic feet per second)
            hourlyHydroData1 = readCSVto2dList(os.path.join(dataDir,'hourlyTurbineFlow01To17Part1.csv'))
            hourlyHydroData2 = readCSVto2dList(os.path.join(dataDir,'hourlyTurbineFlow01To17Part2.csv'))

            #Merge datasets
            hourlyHydroData2 = [row[1:] for row in hourlyHydroData2]
            hourlyHydroData = [row1 + row2 for row1,row2 in zip(hourlyHydroData1,hourlyHydroData2)]
            dateCol = hourlyHydroData[0].index('Date')
            dam1Col = dateCol+1

            #Add columns w/ datetime info
            hourlyHydroData[0] += ['Year','Month','Day','Hour']
            yCol = hourlyHydroData[0].index('Year')
            mCol,dCol,hCol = yCol+1,yCol+2,yCol+3
            for row in hourlyHydroData[1:]:
                currDate = row[dateCol] #M/D/Y H:MM
                splitSlash = currDate.split('/')
                m,d,rest = int(splitSlash[0]),int(splitSlash[1]),splitSlash[2]
                splitYearAndTime = rest.split(' ')
                y,t = int(splitYearAndTime[0]),splitYearAndTime[1]
                h = int(t.split(':')[0])
                row += [y,m,d,h]

            #Isolate year and month data 
            years = [row[yCol] for row in hourlyHydroData[1:]]
            hydroSlim = [hourlyHydroData[0]] + [row for row in hourlyHydroData[1:] if row[yCol]==year]
            months = [row[mCol] for row in hydroSlim[1:]]
            hydroSlim = [hydroSlim[0]] + [row for row in hydroSlim[1:] if row[mCol]==month]
            
            years = [row[yColDemand] for row in hourlyDemandData[1:]]
            demandSlim = [hourlyDemandData[0]] + [row for row in hourlyDemandData[1:] 
                                                    if row[yColDemand]==year]
            months = [row[mColDemand] for row in demandSlim[1:]]
            demandSlim = [demandSlim[0]] + [row for row in demandSlim[1:] if row[mColDemand]==month]
            
            #Plot data
            fig = plt.figure(figCtr,figsize=(25,20))
            rows,cols = 5,6
            #Loop through each dam
            ctr = 0
            for col in range(dam1Col,yCol):
                #Sum daily turbine flows
                dailyGen,dailyDemand = list(),list()
                days = daysInMonth[month]
                for day in range(1,days+1):
                    dayRows = [i for i in range(len(hydroSlim)) if hydroSlim[i][dCol] == day]
                    gen = sum([int(hydroSlim[idx][col]) for idx in dayRows 
                                        if hydroSlim[idx][col] != ''])
                    dayRowsDemand = [i for i in range(len(demandSlim)) 
                                        if demandSlim[i][dColDemand] == day]
                    demand = sum([int(demandSlim[idx][loadCol]) for idx in dayRowsDemand 
                                        if demandSlim[idx][loadCol] != ''])
                    dailyGen.append(gen)
                    dailyDemand.append(demand)
                #Plot data
                ax = plt.subplot(rows,cols,col)
                ax.scatter(dailyDemand,dailyGen)
                #Adjust axis labels and titles
                plt.title(hydroSlim[0][col])
                if ctr%cols == 0: plt.ylabel('Daily Turb. Flow (cfs)')
                if col>(cols*(rows-1)): plt.xlabel('Daily Demand (MWh)')
                else: ax.axes.get_xaxis().set_ticks([])
                ctr += 1

            fig.savefig(fileName+'.png',dpi=100,transparent=True, bbox_inches='tight', pad_inches=0.1)
        
    plt.show()

##################### PLOT TOTAL GENERATION VERSUS TOTAL DEMAND BY WEEK #########
#Plot weekly demand versus generation for all months in year
if plotWeeklyGenAndDemand == True:
    #Iterate over each year
    for year in yearListInput:
        figCtr += 1  
        fileName = 'hydroWeeklyGenVersusDemandY' + str(year)

        #PROCESS DEMAND DATA
        #Import demand data
        hourlyDemandData = readCSVto2dList(os.path.join(demandDataDir,'tvademand.csv'))

        #Add columns w/ datetime info
        dateCol = hourlyDemandData[0].index('time')
        loadCol = hourlyDemandData[0].index('load')
        hourlyDemandData[0] += ['Year','Month','Day','Hour']
        yColDemand = hourlyDemandData[0].index('Year')
        mColDemand,dColDemand,hColDemand = yColDemand+1,yColDemand+2,yColDemand+3
        for row in hourlyDemandData[1:]:
            currDate = row[dateCol] #M/D/Y H:MM
            splitSlash = currDate.split('/')
            m,d,rest = int(splitSlash[0]),int(splitSlash[1]),splitSlash[2]
            splitYearAndTime = rest.split(' ')
            y,t = int(splitYearAndTime[0]),splitYearAndTime[1]
            h = int(t.split(':')[0])
            row += [y,m,d,h]

        #PROCESS HYDRO GENERATION DATA
        #Read data and trim first col (turbine flow is in cubic feet per second)
        hourlyHydroData1 = readCSVto2dList(os.path.join(dataDir,'hourlyTurbineFlow01To17Part1.csv'))
        hourlyHydroData2 = readCSVto2dList(os.path.join(dataDir,'hourlyTurbineFlow01To17Part2.csv'))

        #Merge datasets
        hourlyHydroData2 = [row[1:] for row in hourlyHydroData2]
        hourlyHydroData = [row1 + row2 for row1,row2 in zip(hourlyHydroData1,hourlyHydroData2)]
        dateCol = hourlyHydroData[0].index('Date')
        dam1Col = dateCol+1

        #Add columns w/ datetime info
        hourlyHydroData[0] += ['Year','Month','Day','Hour']
        yCol = hourlyHydroData[0].index('Year')
        mCol,dCol,hCol = yCol+1,yCol+2,yCol+3
        for row in hourlyHydroData[1:]:
            currDate = row[dateCol] #M/D/Y H:MM
            splitSlash = currDate.split('/')
            m,d,rest = int(splitSlash[0]),int(splitSlash[1]),splitSlash[2]
            splitYearAndTime = rest.split(' ')
            y,t = int(splitYearAndTime[0]),splitYearAndTime[1]
            h = int(t.split(':')[0])
            row += [y,m,d,h]

        #Isolate year and month data 
        years = [row[yCol] for row in hourlyHydroData[1:]]
        hydroSlim = [hourlyHydroData[0]] + [row for row in hourlyHydroData[1:] if row[yCol]==year]
        months = [row[mCol] for row in hydroSlim[1:]]
        hydroSlim = [hydroSlim[0]] + [row for row in hydroSlim[1:] if row[mCol] in monthList]
        
        years = [row[yColDemand] for row in hourlyDemandData[1:]]
        demandSlim = [hourlyDemandData[0]] + [row for row in hourlyDemandData[1:] 
                                                if row[yColDemand]==year]
        months = [row[mColDemand] for row in demandSlim[1:]]
        demandSlim = [demandSlim[0]] + [row for row in demandSlim[1:] if row[mColDemand] in monthList]
        
        #Set date of first Monday in July
        firstMondayDays = {2014:7,2015:6}
        firstMondayRow = [row[dColDemand] for row in demandSlim].index(firstMondayDays[year])

        #Plot data
        fig = plt.figure(figCtr,figsize=(25,20))
        rows,cols = 5,6
        #Loop through each dam
        ctr = 0
        for col in range(dam1Col,yCol):
            #Sum daily turbine flows and demand
            weeklyGen,weeklyDemand = list(),list()
            for weekStartIdx in range(firstMondayRow,len(demandSlim)-7*hoursPerDay,7*hoursPerDay):
                gen = sum([int(hydroSlim[idx][col]) for idx in range(weekStartIdx,weekStartIdx+7*hoursPerDay)
                                if hydroSlim[idx][col] != ''])
                demand = sum([int(demandSlim[idx][loadCol]) for idx in range(weekStartIdx,weekStartIdx+7*hoursPerDay)
                                if demandSlim[idx][loadCol] != ''])
                weeklyGen.append(gen)
                weeklyDemand.append(demand)
            #Plot data
            ax = plt.subplot(rows,cols,col)
            ax.scatter(weeklyDemand,weeklyGen)
            #Adjust axis labels and titles
            plt.title(hydroSlim[0][col])
            if ctr%cols == 0: plt.ylabel('Weekly Turb. Flow (cfs)')
            if col>(cols*(rows-1)): plt.xlabel('Weekly Demand (MWh)')
            else: ax.axes.get_xaxis().set_ticks([])
            ctr += 1

        fig.savefig(fileName+'.png',dpi=100,transparent=True, bbox_inches='tight', pad_inches=0.1)
        
    plt.show()
