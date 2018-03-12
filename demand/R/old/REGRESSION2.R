
# ********* REGRESSION 2: DIVIDE TEMPERATURE EFFECTS BY SEASON ********

library(reshape2)

lin.reg.2 <- lm(demand ~ hour.of.day:season2:type.day + tc.1:season2 + 
                  tc.2:season2 + tc.3:season2 + tc.4:season2 + tc.5:season2 + 
                  tc.6:season2, data=df.final)

# weekday summer
result.summer.weekday <- estimate.load.model.2(lin.reg.2, 
                                               df.final, 
                                               season.name = "Summer", 
                                               type.day = "weekday")
result.summer.weekday <- melt(result.summer.weekday,
                              id.vars = "hour")
result.summer.weekday$season <- "Summer"

# weekday winter
result.winter.weekday <- estimate.load.model.2(lin.reg.2, 
                                               df.final, 
                                               season.name = "Winter", 
                                               type.day = "weekday")
result.winter.weekday <- melt(result.winter.weekday,
                              id.vars = "hour")
result.winter.weekday$season <- "Winter"

result.total <- rbind(result.winter.weekday, result.summer.weekday)

g <- ggplot(result.total,
            aes(x=hour, y=value, colour = variable, linetype = season)) +
  geom_line() + theme_bw() + ylab("Demand (MW)") + xlab("hour of day") +
  guides(colour=guide_legend(title=NULL), linetype=guide_legend(title=NULL)) + 
  theme(legend.position=c(0.1,0.8), 
        legend.background = element_blank(),
        legend.box.just = "left")

png(filename = "compare_estimate_real_load.png", width = 3000, height = 2000,
    res = 300)
print(g)
dev.off()














# tc.coeff.2 <- matrix(data=0, nrow = 3, ncol = 6)
# for (i in 1:length(levels(df.final$season2))) {
#   seas.name <- levels(df.final$season2)[i]
#   tc.coeff.2[i, numbers.tc] <- get.temp.betas(lin.reg.2, seas.name)
# }
# 
# pred.out.2 <- as.matrix(t.c) %*% t(tc.coeff.2)
# pred.out.2 <- pred.out.2 + coef.lin.reg.2[1, 1]
# colnames(pred.out.2) <- levels(df.final$season2)
# 
# hour.effects <- coef.lin.reg.2[grep("hour.of.day", row.names(coef.lin.reg.2)), ]
# 
# mean.season.effect <- matrix(data=0, nrow= nrow(pred.out.2), 
#                              ncol = length(levels(df.final$season2)))
# for (i in 1:length(levels(df.final$season2))) {
#   i.r <- grep(levels(df.final$season2)[i], row.names(hour.effects))
#   season.hour.effects <- hour.effects[i.r, ]
#   mean.season.effect[ ,i] <- mean(season.hour.effects[ ,1])
# }
# 
# pred.out.2 <- pred.out.2 + mean.season.effect 
# 
# png(filename = "temp_load_estcurve2.png", width = 2000, height = 2000,
#     res = 300)
# plot(1,type="n")
# plot(x=df.final$temperature, 
#      y=df.final$demand, 
#      pch = 19, cex = 0.4, 
#      col = rgb(0.5, 0.5, 0.5, 0.15),
#      xlab = expression(paste("Temperature (", degree, "C)")),
#      ylab = "Demand (MW)")
# 
# col.line <- c("blue", "red", "black")
# for (i in 1:ncol(pred.out.2)) {
#   lines(temp.plot, pred.out.2[ ,i], type="l", col=col.line[i], lwd = 2, 
#         lty = "dashed")
# }
# legend("topleft", legend = colnames(pred.out.2), 
#        lty = "dashed", col = col.line)
# dev.off()