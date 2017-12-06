# Oficial script with functions for the 
# regression model that creates the projections of electricity demand
# RIPS project
# May 2017

# load other scripts and libraries ----
library(plyr)
library(dplyr)
library(tidyr)

source("auxiliary_functions.R")
source("ripsSpecialPlots.R")
source("prepare_data_frame.R")
source("regression4.R")
source("eiaData.R")
source("NBER.R")

source('Read_SSP_Data.R')
source('readUW.R')
source("CENSUS_DATA.R")
source('BEA_DATA.R')

read.data.frame <- function(read.from.file=TRUE, 
                            name.of.file = "./data/RipsDemand.rds",
                            clean.tva=TRUE) {
  # reads data frame with historical data used in the linear regression
  #
  #
  
  if(read.from.file) {
    df.final <- readDataFromFile(name.of.file = name.of.file)
  } else {
    df.final <- readAllData()
  }

  if (clean.tva){
    # clean data points in 2003-04-06 at 1AM 
    # (load drops to 4000 MW with no reason)
    df.final$load[df.final$load < 5000] <- NA
  }
  
  return(df.final)
}

fit.lm.model <- function(df, temp.breaks=c(0, 10, 20, 30),
                         type.humidity="dp") {
  # wrapper function to fit hourly regression model
  
  lm.hourly <- regression.model.1(data.set = df,
                                  temp.breaks = temp.breaks,
                                  type.humidity = type.humidity,
                                  idx.annual = NULL)
  
  return(lm.hourly)
}

get.ng.projections <- function(year.base=2015) {
  
  df.final <- read.data.frame()
  
  ng.price.base <- unique(unlist(df.final %>% dplyr::filter(year==year.base) %>% 
                                   select(ng.price.real)))
  
  # projections of natural gas price from EIA
  eia.ng.series <- "AEO.2016.REF2016.PRCE_REAL_RES_NA_NG_NA_NA_Y13DLRPMMBTU.A"
  api.call <- getApiCall(eia.ng.series)
  myfile <- fromJSON(api.call)
  proj.ng <- as.data.frame(myfile$series$data[[1]], stringsAsFactors = FALSE)
  names(proj.ng) <- c("YEAR", "VALUE")
  proj.ng <- transform(proj.ng, YEAR = as.numeric(YEAR),
                       VALUE = as.numeric(VALUE))
  proj.ng <- proj.ng[order(proj.ng$YEAR), ]
  proj.ng <- proj.ng[proj.ng$YEAR >= year.base, ]
  proj.ng$change <- proj.ng$VALUE/proj.ng$VALUE[1]
  proj.ng$ng.price <- ng.price.base * proj.ng$change
  
  return(proj.ng)
}


get.ssp.pop <- function(year.base=2015) {

  # ** Population projections
  
  x <- getMultipleStatesPop(states.list = c("US", "TN"))
  
  pop.tn.2015 <- unlist(x %>% dplyr::filter(year==year.base) %>% select(TN))
  
  z <- read.ssp()
  
  ssp.pop <- z[[1]]
  ssp.pop$SCENARIO <- as.factor(substr(ssp.pop$SCENARIO, 1, 4))
  ssp.pop <- ssp.pop[!is.na(ssp.pop$pop), ]
  ssp.pop <- ssp.pop %>% dplyr::select(YEAR, SCENARIO, pop)
  ssp.pop <- ssp.pop %>% dplyr::group_by(SCENARIO) %>% 
    dplyr::mutate(change=pop/pop[2],
                  pop.2=change*pop.tn.2015)
  ssp.pop <- as.data.frame(ssp.pop)
  
  return(ssp.pop)
}

get.ssp.gdp <- function(year.base=2015) {

  # ** GDP projections
  
  y <- getMultipleStatesGDP(states.list = c("US", "TN"))
  y <- y[complete.cases(y), ]
  gdp.tn.2015 <- unlist(y %>% dplyr::filter(year==year.base) %>% select(TN))
  z <- read.ssp()
  ssp.gdp <- z[[2]]
  ssp.gdp$SCENARIO <- as.factor(substr(ssp.gdp$SCENARIO, 1, 4))
  ssp.gdp <- ssp.gdp[!is.na(ssp.gdp$gdp), ]
  ssp.gdp <- ssp.gdp %>% dplyr::select(YEAR, SCENARIO, gdp)
  ssp.gdp <- ssp.gdp %>% dplyr::group_by(SCENARIO) %>% 
    dplyr::mutate(change=gdp/gdp[2],
                  gdp.2=change*gdp.tn.2015)
  ssp.gdp <- as.data.frame(ssp.gdp)
  return(ssp.gdp)
}


project.demand <- function(lm.hourly, socio.var=FALSE, 
                           year.base=2015, 
                           year.proj=2040,
                           temp.breaks=c(0, 10, 20, 30)) {
  
  uw.data <- read.uw(years.data=year.proj)
  uw.data <- compute.calendar.variables(uw.data)

  uw.data.1 <- dplyr::filter(uw.data, year==year.proj)
  uw.data.1 <- compute.calendar.variables(uw.data.1)

  if (socio.var) {
    # projections of demand with socioeconomic variables
    uw.data.1$year <- year.proj
  } else {
    # projection of demand without socioeconomic variables
    uw.data.1$year <- year.base
  }
  
  temp.components <- createTempComponents(uw.data.1$temp, temp.breaks)
  uw.data.1 <- cbind(uw.data.1, temp.components)
  rm(temp.components)

  if (socio.var) {
    # projections of demand with socioeconomic variables
    df.projection <- dplyr::filter(uw.data, year==year.proj)
    temp.components <- createTempComponents(df.projection$temp, temp.breaks)
    df.projection <- cbind(df.projection, temp.components)
    
    ssp.gdp <- get.ssp.gdp()
    ssp.pop <- get.ssp.pop()
    proj.ng <- get.ng.projections(df.final=df.final)
    
    gdp.proj <- unlist(ssp.gdp %>%
                         filter(YEAR == year.proj &
                                  SCENARIO == scenario.names[idx.ssp]) %>%
                         select(gdp.2))
    pop.proj <- unlist(ssp.pop %>%
                         filter(YEAR == year.proj &
                                  SCENARIO == scenario.names[idx.ssp]) %>%
                         select(pop.2))
    gdp.cap.tn <- gdp.proj*1e3/pop.proj
    
    ng.price.real <- unlist(proj.ng %>% filter(YEAR == year.proj) %>%
                              select(ng.price))
    
    df.projection$TN.GDP <- gdp.proj
    df.projection$TN.POP <- pop.proj
    df.projection$us.energyprod <- NA
    df.projection$ng.price.real <- ng.price.real
    df.projection$year <- 1993
    load.proj <- compute.performance(df.projection, lm.hourly$model,
                                     annual.model = lm.hourly$annual.model)$y.hat
  } else {
    # projection of demand without socioeconomic variables
    load.proj <- compute.performance(uw.data.1, lm.hourly$model, NULL)
  }
  load.proj$y.hat$year <- year.proj
#  load.curve.proj <- load.proj$y.hat %>% group_by(year, season, hour.of.day) %>% 
#    summarize(mean.value = mean(y.hat, na.rm=TRUE))
  load.proj <- load.proj$y.hat
  
  load.proj <- load.proj %>% select(time, hour.of.day,
                                    season, y.hat)
  names(load.proj)[4] <- "load"
  
  return(load.proj)
  
}

create.base.line <- function(name.file='data_memphis_1973.rds', 
                             period=c(2005, 2015)) {
  require(plyr)
  require(dplyr)
  require(lubridate)

  df <- readRDS(name.file)
  
  base.line.data <- df %>%
    filter(year(time) >= period[1] & year(time) <= period[2])
  
  return(base.line.data)
}

create.projection.weather <- function(url.uw="~/GoogleDrive/CMU/RIPS/UW/meteo_memphis",
                                      base.line.data, 
                                      type=c('original', 'paulina'),
                                      ref.period=c(2005, 2015), 
                                      future.period=c(2035, 2045)) {
  ## evaluate choices
  type <- match.arg(type)
  
  uw.data <- read.uw(url.uw, years.data=NULL)
  
  if (type == 'original') {
    uw.temp <- uw.data %>% transmute(time=time, value=temp)
    uw.dp <- uw.data %>% transmute(time=time, value=dp)
    base.line.temp <- base.line.data %>% transmute(time=time, value=temp)
    base.line.dp <- base.line.data %>% transmute(time=time, value=dp)
    
    uw.temp.2 <- data.transformation(data=uw.temp, base.line=base.line.temp, 
                                     ref.period=ref.period,
                                     future.period=future.period)
    uw.dp.2 <- data.transformation(data=uw.dp, base.line=base.line.dp,
                                   ref.period=ref.period,
                                   future.period=future.period)
    names(uw.temp.2)[2] <- 'temp'
    names(uw.dp.2)[2] <- 'dp'
    
    uw.data.1 <- data.frame(uw.temp.2, dp=uw.dp.2$dp)
  } else {
    uw.temp <- uw.data %>% transmute(time=time, value=temp)
    uw.dp <- uw.data %>% transmute(time=time, value=dp)
    base.line.temp <- base.line.data %>% transmute(time=time, value=temp)
    base.line.dp <- base.line.data %>% transmute(time=time, value=dp)
    
    uw.temp.2 <- data.transformation.paulina(data=uw.temp, 
                                             base.line=base.line.temp,
                                             ref.period=ref.period,
                                             future.period=future.period)
    uw.dp.2 <- data.transformation.paulina(data=uw.dp, base.line=base.line.dp,
                                           ref.period=ref.period,
                                           future.period=future.period)
    names(uw.temp.2)[2] <- 'temp'
    names(uw.dp.2)[2] <- 'dp'
    
    uw.data.1 <- data.frame(uw.temp.2, dp=uw.dp.2$dp)
  }
  
  uw.data.1 <- uw.data.1 %>% select(time, temp, dp)
  uw.data.1 <- compute.calendar.variables(uw.data.1)
  
  return(uw.data.1)
  
}

project.demand.2 <- function(lm.hourly, 
                             socio.var=FALSE, 
                             year.base=2015, 
                             year.proj=2040,
                             temp.breaks=c(0, 10, 20, 30),
                             weather.data) {
# Computes demand forecast using the transformation of the weather data
# suggested by aufhammer.
#
  
  if (socio.var) {
    # projections of demand with socioeconomic variables
    weather.data$year <- year.proj
  } else {
    # projection of demand without socioeconomic variables
    weather.data$year <- year.base
  }
  
  temp.components <- createTempComponents(weather.data$temp, temp.breaks)
  weather.data <- cbind(weather.data, temp.components)
  rm(temp.components)
  
  if (socio.var) {
    # projections of demand with socioeconomic variables
    df.projection <- dplyr::filter(weather.data, year==year.proj)
    temp.components <- createTempComponents(df.projection$temp, temp.breaks)
    df.projection <- cbind(df.projection, temp.components)
    
    ssp.gdp <- get.ssp.gdp()
    ssp.pop <- get.ssp.pop()
    proj.ng <- get.ng.projections(df.final=df.final)
    
    gdp.proj <- unlist(ssp.gdp %>%
                         filter(YEAR == year.proj &
                                  SCENARIO == scenario.names[idx.ssp]) %>%
                         select(gdp.2))
    pop.proj <- unlist(ssp.pop %>%
                         filter(YEAR == year.proj &
                                  SCENARIO == scenario.names[idx.ssp]) %>%
                         select(pop.2))
    gdp.cap.tn <- gdp.proj*1e3/pop.proj
    
    ng.price.real <- unlist(proj.ng %>% filter(YEAR == year.proj) %>%
                              select(ng.price))
    
    df.projection$TN.GDP <- gdp.proj
    df.projection$TN.POP <- pop.proj
    df.projection$us.energyprod <- NA
    df.projection$ng.price.real <- ng.price.real
    df.projection$year <- 1993
    load.proj <- compute.performance(df.projection, lm.hourly$model,
                                     annual.model = lm.hourly$annual.model)$y.hat
  } else {
    # projection of demand without socioeconomic variables
    load.proj <- compute.performance(weather.data, lm.hourly$model, NULL)
  }
  load.proj$y.hat$year <- year.proj
  #  load.curve.proj <- load.proj$y.hat %>% group_by(year, season, hour.of.day) %>% 
  #    summarize(mean.value = mean(y.hat, na.rm=TRUE))
  load.proj <- load.proj$y.hat
  
  load.proj <- load.proj %>% select(time, hour.of.day,
                                    season, y.hat)
  names(load.proj)[4] <- "load"
  
  # clean up NAs
  x <- load.proj$load
  x[is.na(x)] <- 0 
  cleaned.values <- stats::filter(x=x, filter=c(1/2, 0, 1/2), 
                                  method = 'convolution')
  load.proj$load[is.na(load.proj$load)] <- cleaned.values[is.na(load.proj$load)]
  
  # remove potential outliers due to strange values in temp and dp
  x <- load.proj$load
  m <- mean(x, na.rm=TRUE)
  idx.rows <- (x > 3*m)
  cleaned.values <- stats::filter(x=x, filter=c(1/2, 0, 1/2), 
                                  method = 'convolution')
  load.proj$load[idx.rows] <- cleaned.values[idx.rows]
  
  return(load.proj)
}



write.file.coef.temp <- function(lm.hourly, name.file='temp_coef.csv') {
  # Writes file with estimated coefficients for each of the temperature bins
  #
  # Args:
  #   lm.hourly: list created by function 'fit.lm.model' with the linear 
  #             regression object
  #   name.file: string with name of output file
  #
  # Returns:
  #   nothing
  
  names.coef <- names(lm.hourly$model$coefficients)
  
  # regular expression
  # matches 'tc.' followed by a digit and one optional digit
  # the '$' forces this to be in the end of the string
  idx.coef.temp <- grep('(tc\\.\\d(\\d)?)$', names.coef, perl=TRUE)
  
  coef.temp <- lm.hourly$model$coefficients[idx.coef.temp]
  
  names.bin <- lm.hourly$temp.bins
  
  # creates data frame and writes to file
  df.out <- data.frame(bin=names.bin,
                       value=coef.temp,
                       row.names = NULL)
  
  write.csv(df.out, file = name.file, row.names=FALSE)
}



write.file.coef.interaction <- function(lm.hourly, 
                                        name.file='interaction_coef.csv') {
  # Writes file with estimated coefficients for the interaction btw temperature
  # and air humidity
  #
  # Args:
  #   lm.hourly: list created by function 'fit.lm.model' with the linear 
  #             regression object
  #   name.file: string with name of output file
  #
  # Returns:
  #   nothing
  
  names.coef <- names(lm.hourly$model$coefficients)
  
  # regular expression
  # matches 'tc.' followed by a digit, one optional digit and a collon ':'
  idx.coef.temp <- grep('(tc\\.\\d(\\d)?:)', names.coef, perl=TRUE)
  
  coef.temp <- lm.hourly$model$coefficients[idx.coef.temp]
  
  names.bin <- lm.hourly$temp.bins
  
  # creates data frame and writes to file
  df.out <- data.frame(bin=names.bin,
                       value=coef.temp,
                       row.names = NULL)
  
  write.csv(df.out, file = name.file, row.names=FALSE)
}


write.file.fix.eff.hour <- function(lm.hourly, 
                                    name.file='fix_eff_hour.csv') {
  # Writes file with estimated hourly fixed effects
  #
  # Args:
  #   lm.hourly: list created by function 'fit.lm.model' with the linear 
  #             regression object
  #   name.file: string with name of output file
  #
  # Returns:
  #   nothing
  
  idx.hour.of.day <- grep("hour.of.day", names(lm.hourly$model$coefficients))
  
  coef.hour <- lm.hourly$model$coefficients[idx.hour.of.day]
  names.coef.hour <- names(lm.hourly$model$coefficients)[idx.hour.of.day]
  
  # extract hour of day from names
  # get final position of the hour of day value (use ':' as reference)
  pos2 <- regexpr(":", names.coef.hour) - 1
  pos1 <- rep(12, length(pos2))
  
  hour.of.day <- substr(names.coef.hour, pos1, pos2)
  
  #extract type of day
  pos1 <- regexpr("type.day", names.coef.hour)
  pos1 <- pos1 + attr(pos1, 'match.length')
  attributes(pos1) <- NULL
  names.aux <- substring(names.coef.hour, pos1)
  pos1 <- 1
  pos2 <- regexpr(":", names.aux) - 1
  
  type.day <- substr(names.aux, pos1, pos2)
  
  # extract season  
  pos1 <- regexpr("season", names.coef.hour)
  pos1 <- pos1 + attr(pos1, 'match.length')
  attributes(pos1) <- NULL
  names.aux <- substring(names.coef.hour, pos1)

  season <- names.aux

  # sets values with NA to zero
  coef.hour[is.na(coef.hour)] <- 0
  
  # creates data frame and writes to file
  df.out <- data.frame(season=season,
             type.day=type.day,
             hour.of.day=hour.of.day,
             value=coef.hour, row.names = NULL)
  
  write.csv(df.out, file = name.file, row.names=FALSE)
}

write.file.fix.eff.year <- function(lm.hourly, 
                                    name.file='fix_eff_year.csv') {
  # Writes file with estimated annual fixed effects
  #
  # Args:
  #   lm.hourly: list created by function 'fit.lm.model' with the linear 
  #             regression object
  #   name.file: string with name of output file
  #
  # Returns:
  #   nothing
  
  idx.year <- grep("year", names(lm.hourly$model$coefficients))
  
  fe.year <- lm.hourly$model$coefficients[idx.year]
  names.fe.year <- names(lm.hourly$model$coefficients)[idx.year]
  
  # extract year from names
  pos1 <- regexpr("\\)", names.fe.year) + 1

  year <- substring(names.fe.year, pos1)
  
  # creates data frame and writes to file
  df.out <- data.frame(year=year,
                       value=fe.year, 
                       row.names = NULL)
  
  write.csv(df.out, file = name.file, row.names=FALSE)
}

write.file.intercept <- function(lm.hourly, name.file='intercept.csv') {
  # Writes file with estimated intercept of the model
  #
  # Args:
  #   lm.hourly: list created by function 'fit.lm.model' with the linear 
  #             regression object
  #   name.file: string with name of output file
  #
  # Returns:
  #   nothing
  
  idx.intercept <- grep("Intercept", names(lm.hourly$model$coefficients))
  
  intercept <- lm.hourly$model$coefficients[idx.intercept]

  # creates data frame and writes to file
  df.out <- data.frame(value=intercept, 
                       row.names = NULL)
  
  write.csv(df.out, file = name.file, row.names=FALSE)
}


# temp.2015 <- df.data %>% filter(year == 2015) %>% 
#   select(season, hour.of.day, temp) %>% 
#   group_by(season, hour.of.day) %>% 
#   summarise(temp=mean(temp, na.rm=TRUE))
# 
# g <- ggplot()
# g <- g + geom_line(data=temp.2015, 
#                    aes(x=as.numeric(hour.of.day), y=temp)) + 
#   facet_wrap(~season)
# g
# 
# bb <- base.line.data %>% mutate(hour.of.day=hour(time),
#                                 season=getSeason(time)) %>%
#   select(season, hour.of.day, temp) %>% 
#   group_by(season, hour.of.day) %>% 
#   summarise(temp=mean(temp, na.rm=TRUE))
# 
# g <- g + geom_line(data=bb, aes(x=as.numeric(hour.of.day), y=temp),
#                    linetype='dashed') + facet_wrap(~season)
# g


# base line data averaging 10 years
# base.line.data.2 <- compute.calendar.variables(base.line.data)
# base.line.data.2 <- base.line.data.2 %>% 
#   filter(year(time)<=2015 & year(time) >= 2005) %>% 
#   select(time, temp, dp, hour.of.day, season) %>% 
#   mutate(month=month(time), day=day(time)) %>%
#   group_by(month, day, hour.of.day) %>%
#   summarise(temp=mean(temp, na.rm=TRUE), dp=mean(dp, na.rm=TRUE), 
#             season=season[1]) %>% 
#   mutate(time=ymd(paste(2012, month, day, sep='-'))) %>% as.data.frame()
# 
# data.1 <- uw.data %>% filter(year(time) >= 2005 & year(time) <= 2015)
# data.1$time <- data.1$time + years(30)
# 
# data.2 <- uw.data %>% filter(year(time) >= 2035 & year(time) <= 2045)
# 
# x <- left_join(data.1, data.2, by='time')
# delta.temp <- x$temp.y - x$temp.x
# delta.dp <- x$dp.y - x$dp.x
# 
# uw.data.1 <- data.frame(time=data.1$time,
#                         delta.temp=delta.temp,
#                         delta.dp=delta.dp)
# uw.data.1 <- uw.data.1 %>% 
#   mutate(month=month(time), day=day(time), hour.of.day=as.factor(hour(time))) %>%
#   left_join(base.line.data.2, by=c('month', 'day', 'hour.of.day')) %>%
#   transmute(time=time.x, hour.of.day=hour.of.day, 
#             temp=temp+delta.temp, dp=dp+delta.dp) %>% filter(!is.na(time))

