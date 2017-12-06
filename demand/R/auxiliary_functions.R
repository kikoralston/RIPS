# ************* Auxiliary functions ************************

library(ggplot2)
library(plyr)
library(dplyr)
library(lazyeval)

# --------------------------------------------------------------------
# Definitions of auxiliary functions

table.markdown <- function(df, name.file=''){
  # writes data frame in markdown table format to a txt file
  
  nc <- ncol(df)
  
  cat('| ', paste(names(df), collapse =" | "), ' |\n',
      '|', rep('----- |', nc), '\n', sep='', file = name.file)
  for (i in 1:nrow(df)) {
    cat('|', paste(df[i, ], collapse =" | "), '| \n', file = name.file)
  }
}

getSeason <- function(DATES) {
  # Find which season a particular date (or vector of dates) belongs to
  #
  # Args:
  #   DATES: vector of DATES ("POSIXlt" and "POSIXct" formats)
  #
  # Returns:
  #   vector with names of seasons (string)
  
  WS <- as.Date("2012-12-15", format = "%Y-%m-%d", tz='UTC') # Winter Solstice
  SE <- as.Date("2012-3-15",  format = "%Y-%m-%d", tz='UTC') # Spring Equinox
  SS <- as.Date("2012-6-15",  format = "%Y-%m-%d", tz='UTC') # Summer Solstice
  FE <- as.Date("2012-9-15",  format = "%Y-%m-%d", tz='UTC') # Fall Equinox
  
  # Convert dates from any year to 2012 dates
  d <- as.Date(strftime(DATES, format="2012-%m-%d", tz='UTC'))
  
  s <- ifelse (d >= WS | d < SE, "Winter", 
               ifelse (d >= SE & d < SS, "Spring",
                       ifelse (d >= SS & d < FE, "Summer", "Fall")))
  
  return(s)
}

createTempComponents <- function(temp,
                                 temp.breaks = c(-10, 0, 10, 20, 30)){
  # Computes temperature components for piecewise linear model
  #
  # Args:
  #   temp: vector of temperatures (degrees celsius)
  #   temp.breaks: temperature break points of piecewise linear function
  #
  # Returns:
  #   matrix with components in columns
  
  # replace NA with -9999
  temp[is.na(temp)] <- -9999
  
  nr <- length(temp)
  nb <- length(temp.breaks) # number of bounds
  nc <- nb + 1 # number of components = number bounds + 1
  
  temp.comp <- as.data.frame(matrix(data = 0, nrow = nr, ncol = nc))
  names(temp.comp) <- paste0("tc.", (1:nc))
  
  temp.comp[temp < temp.breaks[1], "tc.1"] <- temp[temp < temp.breaks[1]]
  temp.comp[temp >= temp.breaks[1], "tc.1"] <- temp.breaks[1]
  
  for(i in 2:nb) {
    idx.rows <- (temp < temp.breaks[i]) & (temp >= temp.breaks[i-1])
    temp.comp[idx.rows, i] <- temp[idx.rows] - temp.breaks[i-1]
    idx.rows <- (temp >= temp.breaks[i])
    temp.comp[idx.rows, i] <- temp.breaks[i] - temp.breaks[i-1]
  }
  idx.rows <- (temp > temp.breaks[nb])
  temp.comp[idx.rows, nc] <- temp[idx.rows] - temp.breaks[nb]
  
  return(temp.comp)
}

mean.value.by.tc <- function(df.data, name.var,
                             temp.breaks = c(-10, 0, 10, 20, 30)) {
  # Computes mean value of variable for each temperature bin
  # (necessary in order to plot when interaction terms are included)
  #
  # Args:
  #   df.data: data frame with data set
  #   name.var: name of variable
  #   temp.breaks: breakpoints of temperature intervals in stepwise linear model
  #
  # Returns:
  #   vector with mean dew point values for each temperature
  
  # include -Infinity and +Infinity in break values because of cut function
  temp.breaks2 <- c(-Inf, temp.breaks, Inf)
  
  # divides the range of temperature into intervals according to temp.breaks2
  temp.intervals <- cut(df.data$temp, temp.breaks2) 
  
  # creates data frame of ranges and dew.points
  df.1 <- data.frame(temp.intervals=temp.intervals, 
                    temp.name=df.data[, name.var])
  names(df.1)[2] <- name.var
  
  dots <- list(interp(~mean(var, na.rm=TRUE), var=as.name("dp")),
               interp(~quantile(var, 0.05, na.rm=TRUE), 
                      var=as.name("dp")),
               interp(~quantile(var, 0.95, na.rm=TRUE), 
                      var=as.name("dp")))
  
  # computes mean in each temperature range
  df.mean.value <- as.data.frame(df.1 %>% group_by(temp.intervals) %>%
    dplyr::summarize_(.dots = dots))
  names(df.mean.value)[2:4] <- c("mean", "q.5", "q.95")
  
  # drop possible NA group
  df.mean.value <- df.mean.value[complete.cases(df.mean.value), ]
  
  return(df.mean.value)
}

convertDewPoint2RelHum <- function(dew.point, temp) {
  # Converts dew point to relative humidity
  # uses inverse of formula (8) in:
  # http://journals.ametsoc.org/doi/pdf/10.1175/BAMS-86-2-225
  #
  # Args:
  #   dew.point: dew.point value in Celsius. Can be a vector or an atomic
  #   temp: air temperature in Celsius. Can be vector or atomic
  #
  #   OBS: if both arguments are vectors, they must be the same length 
  #
  # Returns:
  #   relative humidity in %

  if (length(dew.point) != length(temp)) {
    if (length(dew.point) > 1 && length(temp) > 1) {
      stop("\nArguments must be:\n",
               "   1) Both vectors of same length; or \n",
               "   2) One vector and one atomic value.\n\n",
               "Please check arguments!")
    }
  }
  
  # definition of constants
  A1 = 17.625 # dimensionless
  B1 = 243.04 # degress celsius
  C1 = 610.94 # Pascal
  
  exp.term <- A1 * dew.point / (B1 + dew.point) - A1 * temp / (B1 + temp)
  
  rh <- sapply(X = exp(exp.term), FUN = min, 1)

  return(100*rh)
}

convertRelHum2DewPoint <- function(rh, temp) {
  # Converts relative humidity to dew point
  # uses inverse of formula (8) in:
  # http://journals.ametsoc.org/doi/pdf/10.1175/BAMS-86-2-225
  #
  # Args:
  #   rht: relative humidity value in % (0 to 100). Can be a vector or an atomic
  #   temp: air temperature in Celsius. Can be vector or atomic
  #
  #   OBS: if both arguments are vectors, they must be the same length 
  #
  # Returns:
  #   dew point in Celsius
  
  if (length(rh) != length(temp)) {
    if (length(rh) > 1 && length(temp) > 1) {
      stop("\nArguments must be:\n",
           "   1) Both vectors of same length; or \n",
           "   2) One vector and one atomic value.\n\n",
           "Please check arguments!")
    }
  }
  
  # definition of constants
  A1 = 17.625 # dimensionless
  B1 = 243.04 # degress celsius
  C1 = 610.94 # Pascal
  
  dp <- B1*(log(rh/100)+A1*temp/(B1+temp)) / 
    (A1-log(rh/100)-A1*temp/(B1+temp)) 
  
  return(dp)
}

# NOT USED ANYMORE ----
# get.temp.betas <- function(lin.reg, season.name = "Winter") {
#   # Gets coefficients of temperature variables
#   #
#   # Args:
#   #   lin.reg: linear regression model
#   #   season.name: (string) name of season
#   #
#   # Returns:
#   #   vector with coefficients
#   
#   # get index of betas for temperature components
#   idx.season <- grep(season.name, names(lin.reg$coefficients))
#   idx.tc.total <- grep("tc.", names(lin.reg$coefficients))
#   idx.tc <- base::intersect(idx.season, idx.tc.total)
#   
#   return(lin.reg$coefficients[idx.tc])
# }
# 
# estimate.load.model.2 <- function(lin.reg, df.data, season.name = "Winter",
#                                   type.day.value = "weekday") {
#   
#   mean.temperature <- ddply(df.data, .(hour.of.day, season2), 
#                             summarise, 
#                             mean = mean(temperature, na.rm = TRUE))
#   t.c.mean <- createTempComponents(
#     mean.temperature[mean.temperature$season2 == season.name, c("mean")])
#   
#   # copy estimated coefficients (betas)
#   betas <- lin.reg$coefficients
#   
#   # set NAs to zero (some coefficients are automatically 
#   # dropped and set to NA to avoid collinearity)
#   betas[is.na(betas)] <- 0
#   
#   # get index of coefficients 
#   idx.season <- grep(season.name, names(betas))
#   idx.day <- grep(type.day.value, names(betas))
#   idx.hour.effect <- base::intersect(idx.season, idx.day)
#   
#   temp.betas <- get.temp.betas(lin.reg, season.name)
#   temp.betas[is.na(temp.betas)] <- 0
#   
#   idx.intercept <- which(names(betas) ==  "(Intercept)")
#   if (length(idx.intercept) == 0) {
#     beta0 <- 0
#   } else {
#     beta0 <- betas[idx.intercept]
#   }
#   
#   # compute estimate of load
#   est.load.season <- beta0 + betas[idx.hour.effect] + 
#     as.matrix(t.c.mean) %*% temp.betas
#   
#   mean.load <- df.data %>%
#     group_by(hour.of.day, type.day, season2) %>%
#     select(demand) %>%
#     summarise(
#       mean.value = mean(demand, na.rm = TRUE)) %>%
#     dplyr::filter(season2 == season.name & type.day == type.day.value)
#   
#   df.est <- data.frame(hour = 0:23, 
#                        estimate = est.load.season,
#                        real.mean = mean.load$mean.value)
#   
#   return(df.est)
# }