##### Creates plots for inferlab website about the project ######

library(grid)
library(gridExtra)

read.from.file <- TRUE

setwd("~/GoogleDrive/CMU/RIPS/R")

source("auxiliary_functions.R")

# prepares data frame for regression
source("prepare_data_frame.R")
if(read.from.file) {
  df.final <- readDataFromFile()
} else {
  df.final <- readAllData()
}

lm1 <- lm(load ~ I(temperature) + I(temperature^2) + I(temperature^3), 
          data = df.final)
lm1.sum <- summary(lm1)

x.temp.data <- cbind(1, c(-15:40), c(-15:40)^2 ,c(-15:40)^3)
est.load <- x.temp.data %*% coef(lm1)
est.load <- data.frame(x=c(-15:40), y = est.load)

g1 <- ggplot()
g1 <- g1 + geom_point(data=df.final, aes(x=temperature, y=load/1000), 
                    alpha=0.2, colour="blue")
g1 <- g1 + geom_line(data=est.load, aes(x=x, y=y/1e3), colour="red", size=1.2)
g1 <- g1 + theme_bw() + xlab(expression(paste("Temperature (", degree, "C)")))
g1 <- g1 + ylab("Hourly Load (GW)")
g1

source("regression4.R")

lm.1 <- regression.model.1(df.final, test.year = 2014)


#print(g)


png(file = "./plotout1.png", width = 2000, height = 1000, res = 300)
grid.arrange(g1, lm.1$gg.plot, ncol=2)
dev.off()

