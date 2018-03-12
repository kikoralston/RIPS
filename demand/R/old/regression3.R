# ************** REGRESSION 3: WITH DEW POINT *********************

library(car)

# component residual plot doesn't work with interactions so I will divide
# the data into seasons ans type of days to check Crplots

df.summer.weekday <- subset(df.final, season2 == "Summer" & 
                              type.day == "weekday")
# only for summer and weekdays
lin.reg.3 <- lm(demand ~ hour.of.day + tc.3 + tc.4 + tc.5 + 
                  tc.6 + dew.point, data=df.summer.weekday)

png(filename = "crPlot_dewpoint_summer.png", width = 2000, height = 2000,
    res = 300)
crPlots(lin.reg.3, terms = "dew.point", col = rgb(0, 0, 0, .5), 
        pch = 19)
dev.off()

df.winter.weekday <- subset(df.final, season2 == "Winter" & 
                              type.day == "weekday")
# only for winter and weekdays
lin.reg.4 <- lm(demand ~ hour.of.day + tc.3 + tc.4 + tc.5 + 
                  dew.point, data=df.winter.weekday)
png(filename = "crPlot_dewpoint_winter.png", width = 2000, height = 2000,
    res = 300)
crPlots(lin.reg.4, terms = "dew.point", col = rgb(0, 0, 0, .5), 
        pch = 19)
dev.off()

