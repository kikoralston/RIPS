
library(ggplot2)

output2file <- TRUE

# w.size <- 90*24 # size of window for moving average
# df.1$ma.3m <- stats::filter(df.1$demand, filter = (1/w.size)*rep(1,w.size), 
#                     sides = 1)
# w.size <- 365*24 # size of window for moving average
# df.1$ma.1y <- stats::filter(df.1$demand, filter = (1/w.size)*rep(1,w.size), 
#                     sides = 1)
# 
# df.summer <- data.frame(beg = seq(from = as.POSIXct("2006-06-22 0:00:00"), 
#                               length = 9, by = "year"),
#                         end = seq(from = as.POSIXct("2006-09-22 0:00:00"), 
#                                   length = 9, by = "year"))
# df.winter <- data.frame(beg = seq(from = as.POSIXct("2006-12-22 0:00:00"), 
#                                   length = 8, by = "year"),
#                         end = seq(from = as.POSIXct("2007-03-22 0:00:00"), 
#                                   length = 8, by = "year"))
# ticks.plot <- seq(from = as.POSIXct("2006-01-01 0:00:00"), 
#                   length = 9, by = "year")
# 
# g <- ggplot()
# g <- g + geom_rect(data=df.summer, 
#                    aes(xmin=beg, xmax=end, ymin=-Inf, ymax=+Inf), 
#                    fill=rgb(220, 220, 220, 130, maxColorValue = 255))
# g <- g + geom_rect(data=df.winter, 
#                    aes(xmin=beg, xmax=end, ymin=-Inf, ymax=+Inf), 
#                    fill=rgb(169, 169, 169, 130, maxColorValue = 255))
# g <- g + geom_line(data=df.1, aes(x=time, y=demand), colour = "darkgray")
# g <- g + geom_line(data=df.1, aes(x=time, y=ma.1y), colour = "black")
# g <- g + theme_bw() + scale_x_datetime(breaks=ticks.plot, date_labels = "%Y")
# 
# if (output2file) {
#   png(filename = "plot_TVA.png", width = 3000, height = 2000, res = 300)
# }
# print(g) # print plot
# if (output2file) {
#   dev.off()
# }
# top.values <- head(sort(df.1$demand, decreasing = TRUE, 
#                            index.return = TRUE)$x, n=50)
# top.pos <- head(sort(df.1$demand, decreasing = TRUE, 
#                            index.return = TRUE)$ix, n=50)
# b <- df.1$time[top.pos]


png(filename = "./out/test.png", width=1000, height=1500)
par(mfrow=c(3,1))
plot(x=df.final$time, y=df.final$load, type="l")
plot(x=df.final$time, y=df.final$temperature, type="l")
plot(x=df.final$time, y=df.final$dew.point, type="l")
dev.off()

# ------------------------------------------------------------------------
# Creates histograms and time-series plots for temperature
g <- ggplot(df.final, 
            aes(x=temperature, fill=season2)) +
  geom_density(position="identity", alpha=0.6, colour = NA) + theme_bw() +
  scale_fill_brewer(type = "qual", palette = 3) +
  guides(fill=guide_legend(title=NULL)) +
  theme(axis.title.y=element_blank(), legend.position=c(0.1,0.9))

png(filename = "hist_temperature.png", width = 2000, height = 2000, res = 300)
print(g)
dev.off()

df.summer <- data.frame(beg = seq(from = as.POSIXct("2006-06-22 0:00:00"), 
                                  length = 9, by = "year"),
                        end = seq(from = as.POSIXct("2006-09-22 0:00:00"), 
                                  length = 9, by = "year"))
df.winter <- data.frame(beg = seq(from = as.POSIXct("2006-12-22 0:00:00"), 
                                  length = 8, by = "year"),
                        end = seq(from = as.POSIXct("2007-03-22 0:00:00"), 
                                  length = 8, by = "year"))
ticks.plot <- seq(from = as.POSIXct("2006-01-01 0:00:00"), 
                  length = 9, by = "year")

df.temp <- df.final[, c("time", "temperature")]
df.temp <- df.temp[!is.na(df.temp$temperature), ]
w.size <- 90*24 # size of window for moving average
df.temp$ma.3m <- stats::filter(df.temp$temperature, filter = (1/w.size)*rep(1,w.size), 
                               sides = 1)
w.size <- 365*24 # size of window for moving average
df.temp$ma.1y <- stats::filter(df.temp$temperature, filter = (1/w.size)*rep(1,w.size), 
                               sides = 1)
g <- ggplot()
g <- g + geom_rect(data=df.summer, 
                   aes(xmin=beg, xmax=end, ymin=-Inf, ymax=+Inf), 
                   fill='yellow', alpha=0.2)
g <- g + geom_rect(data=df.winter, 
                   aes(xmin=beg, xmax=end, ymin=-Inf, ymax=+Inf), 
                   fill='blue', alpha=0.2)
g <- g + geom_line(data=df.temp, aes(x=time, y=temperature), 
                   colour = "darkgray")
g <- g + geom_line(data=df.temp, aes(x=time, y=ma.3m), colour = "red")
g <- g + geom_line(data=df.temp, aes(x=time, y=ma.1y), colour = "blue")
g <- g + theme_bw() + scale_x_datetime(breaks=ticks.plot, date_labels = "%Y")
g <- g + ylab(expression(paste("Temperature (", degree, "C)"))) +
  theme(axis.title.x=element_blank())

png(filename = "time_series_temperature.png", width = 3000, height = 2000, 
    res = 300)
print(g)
dev.off()

g <- ggplot() + 
  geom_point(data=df.final, aes(x=temperature, y=demand), 
             col = rgb(0, 0, 1, 0.3), size=0.4) + theme_bw() + 
  xlab(expression(paste("Temperature (", degree, "C)"))) +
  ylab("Demand (MW)")
png(filename = "temperature_demand.png", width = 2000, height = 2000,
    res = 300)
print(g)
dev.off()

png(filename = "dewpoint_demand.png", width = 2000, height = 2000, res = 300)
plot(x=df.final$dew.point, 
     y=df.final$demand, 
     pch = 19, cex = 0.4, 
     col = rgb(0, 0, 1, 0.3),
     xlab = expression(paste("Dew Point (", degree, "C)")),
     ylab = "Demand (MW)")
dev.off()

png(filename = "dewpoint_temperature.png", width = 2000, height = 2000,
    res = 300)
plot(x=df.final$dew.point, 
     y=df.final$temperature, 
     pch = 19, cex = 0.4, 
     col = rgb(0, 0, 1, 0.3),
     xlab = expression(paste("Dew Point (", degree, "C)")),
     ylab = expression(paste("Temperature (", degree, "C)")))
abline(a=0,b=1)
dev.off()

# ---------------------------------------------------------------
# Create data frames to plot weekdays, weekends, Seasons

mean.load.curve <- ddply(df.final, .(hour.of.day, type.day, season), 
                         summarise, 
                         mean = round(mean(demand), 2),
                         quant.5 = round(quantile(demand, 0.05),2),
                         quant.95 = round(quantile(demand, 0.95),2))

g <- ggplot(subset(mean.load.curve, 
                   type.day == "weekday" & season %in% c("Summer", "Winter")),
            aes(x=as.numeric(hour.of.day)-1, y=mean/1000, colour=season)) + 
  geom_ribbon(aes(ymin=quant.5/1000, ymax=quant.95/1000, fill=season), 
              colour=NA, alpha=0.2) +
  geom_line() + theme_bw() + guides(colour=guide_legend(title=NULL),
                                    fill=guide_legend(title=NULL)) +
  ylab("Demand (GW)") + xlab("hour of day") + 
  theme(legend.position=c(0.1,0.9))

png(filename = "daily_load_curve.png", width = 2000, height = 2000,
    res = 300)
print(g)
dev.off()


# -----------------------------------

annual.demand <- df.final %>% select(year, load) %>%
  group_by(year) %>% 
  summarize(annual.load = mean(load, na.rm = TRUE))
annual.demand <- as.data.frame(annual.demand)

pop.data <- df.final[ ,c(1,17:23)]
pop.data <- unique(pop.data)
nms.col <- names(pop.data)
nms.col <- gsub(".x","",nms.col)
names(pop.data) <- nms.col

gdp.data <- df.final[ ,c(1,24:30)]
gdp.data <- unique(gdp.data)
nms.col <- names(gdp.data)
nms.col <- gsub(".y","",nms.col)
names(gdp.data) <- nms.col

plotGdpDemand(gdp.data, annual.demand)
plot.population(pop.data, annual.demand)

annual.demand.log <- annual.demand
annual.demand.log$annual.load <- log(annual.demand$annual.load)

pop.data.log <- pop.data
pop.data.log[,c(2:8)] <- log(pop.data.log[,c(2:8)])

gdp.data.log <- gdp.data
gdp.data.log[,c(2:8)] <- log(gdp.data.log[,c(2:8)])

plotGdpDemand(gdp.data.log, annual.demand.log)
plot.population(pop.data.log, annual.demand.log)


