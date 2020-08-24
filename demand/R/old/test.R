source('~/GoogleDrive/CMU/RIPS/R/downloadNOAA.R')

listStations(nameCity = "Logan")

pgh.stations <- listStations(nameCity = "Logan")
pgh.data <- readStationFileRemote(pgh.stations$USAF[5], 2013, 2016)
df.test <- pgh.data %>% group_by(date(time)) %>% 
  summarize(min=min(temp,na.rm=TRUE), 
            max=max(temp, na.rm=TRUE))
df.test.2 <- data.frame(year=year(df.test$`date(time)`), df.test)
df.test.2 <- data.frame(year=year(df.test$`date(time)`), df.test)
year(df.test.2$date.time.) <- min(year(df.test.2$date.time.), 
                                  na.rm = TRUE)
df.test.2 <- df.test.2[complete.cases(df.test.2), ]

library(plotly)

g <- ggplot(data=df.test.2)
g <- g + geom_ribbon(aes(x=date.time., ymin=min, ymax=max, 
                         fill=as.factor(year)), alpha=0.45)
g <- g + scale_x_date(date_labels = "%b", date_breaks = "1 months")
g <- g + theme_bw()
g <- g + guides(fill=guide_legend(title=NULL))

ggplotly(g)
