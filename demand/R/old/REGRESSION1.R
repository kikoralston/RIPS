# ***************** REGRESSION 1 ****************

lin.reg <- lm(demand ~ hour.of.day:season2:type.day + tc.1 + tc.2 + 
                tc.3 + tc.4 + tc.5 + tc.6, 
              data=df.final)

# ---------------------------------------------------------------
# Plotting results from Linear Regression

season.test <- "Winter"
type.of.day <- "weekday"
mean.temperature <- ddply(df.final, .(hour.of.day, season2), 
                          summarise, 
                          mean = mean(temperature, na.rm = TRUE))
t.c.mean <- createTempComponents(
  mean.temperature[mean.temperature$season2 == season.test, c("mean")])

idx.season <- grep(season.test, names(lin.reg$coefficients))
idx.day <- grep(type.of.day, names(lin.reg$coefficients))
idx.final <- base::intersect(idx.season, idx.day)

est.load.season <- lin.reg$coefficients[1] + 
  lin.reg$coefficients[idx.final] + 
  as.matrix(t.c.mean) %*% lin.reg$coefficients[2:7]
df.est <- data.frame(x = 0:23, y = est.load.season)

g <- ggplot() +
  geom_line(data=subset(mean.load.curve, 
                        type.day == type.of.day & season %in% c(season.test)),
            aes(x=as.numeric(hour.of.day)-1, y=mean/1000), colour = "red") + 
  geom_line(data = df.est, aes(x=x, y=y/1000), colour = "blue") + theme_bw() +
  ylab("Demand (GW)") + xlab("hour of day")

png(filename = "est_load_curve.png", width = 2000, height = 2000, 
    res = 300)
print(g)
dev.off()

lim.plot <- c(min(lin.reg$model$demand, lin.reg$fitted.values),
              max(lin.reg$model$demand, lin.reg$fitted.values))

png(filename = "fittedvalues.png", width = 2000, height = 2000, res = 300)
plot(x = lin.reg$model$demand, y = lin.reg$fitted.values, 
     pch = 19, cex = 0.3, 
     col = rgb(0, 0, 1, 0.3), 
     xlim = lim.plot, ylim = lim.plot,
     ylab = "Fitted Values", xlab = "Observed Values")
abline(a=0, b=1)
dev.off()

temp.plot <- -20:40
t.c <- createTempComponents(temp.plot)
pred.out <- as.matrix(t.c) %*% lin.reg$coefficients[2:7]
pred.out <- pred.out + lin.reg$coefficients[1]

png(filename = "temp_load_estcurve.png", width = 2000, height = 2000,
    res = 300)
plot(x=df.final$temperature, 
     y=df.final$demand, 
     pch = 19, cex = 0.4, 
     col = rgb(0.5, 0.5, 0.5, 0.15),
     xlab = expression(paste("Temperature (", degree, "C)")),
     ylab = "Demand (MW)")
lines(temp.plot, pred.out, type="l", col="black", lwd = 2, lty = "dashed")
dev.off()