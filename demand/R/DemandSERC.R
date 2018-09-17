library('plyr')
library('dplyr')
library(lubridate)
library(tidyr)

library(ggplot2)
library(gridExtra)
library(grid)

df.companies <- read.csv(file = paste('/Users/kiko/Google Drive/CMU/RIPS/SERC companies/SERC.csv'),
                      stringsAsFactors = FALSE)
df.companies[df.companies$Name.ferc != "", c("Name")] <- df.companies[df.companies$Name.ferc != "", c("Name.ferc")]
companies <- unlist(df.companies %>% select(Name, FERC) %>% filter(FERC == 1) %>% select(Name))

source('readFERC.R')
source('rips_demand.R')
source('auxiliary_functions.R')

ferc.zip <- downloadFERCdata(ferc.zip = 'form714-database.zip')

df.results <- data.frame(name=as.character(),
                         year.ini=as.numeric(),
                         year.end=as.numeric(),
                         n.points=as.numeric(),
                         mean.load=as.numeric())

pb <- txtProgressBar(min = 0, max = length(companies), style = 3)
setTxtProgressBar(pb, 0)

data.ferc <- NULL

for (i in 1:length(companies)) {
  
  aux <- readFercNew(get.online = FALSE,
                     ferc.zip = ferc.zip,
                     name.util = companies[i])
  if (!is.null(aux)) {
    aux$name <- companies[i]
    zone <- df.companies$ipm.zone[df.companies$Name == companies[i]]
    aux$zone <- zone
    
    if(i==1){
      data.ferc <- aux
    } else{
      data.ferc <- rbind(data.ferc, aux)
    }
  }
  
  setTxtProgressBar(pb, i)
}

df.summary.1 <- data.ferc %>% group_by(name) %>% 
  summarize(year.ini=min(year(time)), year.end=max(year(time)),
            n.points=n(), mean.load=mean(load, na.rm = TRUE)) %>% 
  as.data.frame()

df.ferc.final <- data.ferc %>% filter(!is.na(zone)) %>% group_by(time, zone) %>% 
  summarize(load=sum(load)) %>% as.data.frame()

g <- ggplot() + geom_line(data=df.ferc.final, aes(x=time, y=load, colour=zone))

zones <- unique(df.ferc.final$zone)

path.demand <- '~/Documents/CE/Data/DemandData/'

df.zones <- read.csv(paste0(path.demand, 'map_zone_cell.csv'))

for (z in zones) {
  aux <- df.zones %>% filter(zone == z) %>% as.data.frame()
  
  for (i in 1:nrow(aux)) {
    fname <- paste0('RH_airT_19790101_20161231.', aux$X[i], '.csv')
    if (i == 1){
      df.temp <- read.csv(paste0(path.demand, '/trainning/', fname))
      df.temp$city <- aux$X[i]
    } else {
      df.temp.1 <- read.csv(paste0(path.demand, '/trainning/', fname))
      df.temp.1$city <- aux$X[i]
      
      df.temp <- rbind(df.temp, df.temp.1)
    }
  }
  names(df.temp)[1] <- 'time'
  df.temp$time <- as.POSIXct(strptime(as.character(df.temp$time), 
                                      "%Y-%m-%d %H:%M:%S"))
  df.temp <- df.temp %>% group_by(time) %>% summarise(rel_humid=mean(rel_humid),
                                                      air_T=mean(air_T))
  
  df.load <- df.ferc.final %>% filter(zone == z)
  
  df.final <- df.temp %>% left_join(df.load, by='time') %>% select(-zone)
  
  year.end <- max(year(df.ferc.final$time))
  year.ini <- min(year(df.ferc.final$time))
  
  df.final <- df.final %>% 
    transmute(time=time, temp=air_T, rh=rel_humid, 
              dp=convertRelHum2DewPoint(rh, temp), load=load) %>%
    as.data.frame()
  
  df.final <- df.final %>% filter(!is.na(time))
  
  df.final <- df.final %>% arrange(time)
  
  # there are some observations with demand = 0 which is not possible
  # set them to NA (not observed)
  df.final$load[df.final$load == 0] <- NA
  
  # there are some observations where temperature is NA but dew.point 
  # has a numerical value. This can cause error in the interaction term
  df.final$dp[is.na(df.final$temp)] <- NA
  df.final$temp[is.na(df.final$dp)] <- NA
  
  # create new columns for regression
  
  # calendar variables
  df.final <- compute.calendar.variables(df.final)
  df.final <- cbind(data.frame(year=year(df.final$time)), df.final)
  
  df.final <- df.final %>% filter(year >= year.ini & year <= year.end)
  
  df.final <- df.final %>% filter(year <= 2015)
  
  # fit regression model
  cat('\n-------- Fitting Regression model --------\n')
  
  temp.breaks <- c(-30, 0, 10, 20, 30)
  lm.model <- fit.lm.model(df=df.final, temp.breaks=temp.breaks, 
                           type.humidity = NULL)
  
  fname = paste0(path.demand, 'temp_coef_',z,'.csv')
  write.file.coef.temp(lm.hourly=lm.model, name.file=fname)
  
  fname = paste0(path.demand, 'fix_eff_hour_',z,'.csv')
  write.file.fix.eff.hour(lm.hourly=lm.model, name.file=fname)
  
  fname = paste0(path.demand, 'fix_eff_year_',z,'.csv')
  write.file.fix.eff.year(lm.hourly=lm.model, name.file=fname)
  
  fname = paste0(path.demand, 'intercept_',z,'.csv')
  write.file.intercept(lm.hourly=lm.model, name.file=fname)
}



