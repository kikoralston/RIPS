# ************************************************************************
#           functions to prepare and read data frame 
# ************************************************************************

readAllData <- function(save.to.file = TRUE, name.of.file = "",
                        remote.ferc = TRUE, remote.noaa = TRUE) {
  # downloads and compiles data base
  #
  # Args:
  #   save.to.file: TRUE if data base should be saved to file
  #   name.of.file: name of file (complete path)
  #   remote.ferc: (boolean) TRUE if FERC data should be read from 
  #                 remote repositories
  #   remote.noaa: (boolean) TRUE if NOAA data should be read from 
  #                 remote repositories
  #
  # Returns:
  #   df.final: data frame with data base

  # ---------------------------------------------------------------s
  # Sources other scripts
  source("auxiliary_functions.R")
  
  library(lubridate)
  
  # reads electricity demand data downloaded from FERC
  source("readFERC.R")
  df.ferc.1 <- readFercOld(flag.download = remote.ferc)
  df.ferc.2 <- readFercNew(get.online = remote.ferc)
  df.ferc.final <- rbind(df.ferc.1, df.ferc.2)
  rm(df.ferc.1, df.ferc.2)
  
  year.end <- max(year(df.ferc.final$time))
  year.ini <- min(year(df.ferc.final$time))
  
  # reads weather data from NOAA
  source("downloadNOAA.R")
  flag.noaa.local <- !(remote.noaa)
  #df.stations <- listStations(nameCity = "Memphis")
  #df.stations[ ,c(1:5, 10:11)]

  name.city <- "Memphis"
    
  if (flag.noaa.local) {
    # ls.files <- list.files(path = "./data", 
    #                        pattern = ".gz",
    #                        full.names = TRUE)  # list of gz files with weather data
    # df.temp <- readStationFileLocal(ls.files)
    df.weather <- read.csv("./data/weather_data.csv", stringsAsFactors = FALSE)
    df.weather$time <- as.POSIXct(df.weather$time, tz='UTC')
  } else {
    # df.temp <- readStationFileRemote(USAF = "723340", yearini = year.ini, 
    #                                  yearend = year.end)
    df.weather <- compileDataListStations(file.list = "./data/listStations.txt",
                                            year.ini = year.ini,
                                            year.end = year.end)
  }
  
  cols.weather <- grep(name.city, names(df.weather), ignore.case = TRUE,
                       perl = TRUE)
  
  df.weather.2 <- df.weather[, c(1, cols.weather)]
  names(df.weather.2) <- 
    gsub(paste0(".", name.city), "", names(df.weather.2), ignore.case = TRUE)
  
  # merges demand and temperature data
  df.final <- merge(df.ferc.final, df.weather.2, by="time", all = TRUE)
  rm(df.ferc.final)
  
  # ------------------------------------------------------------------------s
  # Prepare data frame for regression
  
  # there are some observations with demand = 0 which is not possible
  # set them to NA (not observed)
  df.final$load[df.final$load == 0] <- NA
  
  # there are some observations where temperature is NA but dew.point 
  # has a numerical value. This can cause error in the interaction term
  df.final$dp[is.na(df.final$temp)] <- NA
  
  # create new columns for regression
  
  # calendar variables
  df.final <- compute.calendar.variables(df.final)
  
  df.final <- cbind(data.frame(year=year(df.final$time)), df.final)
  
  if(FALSE){
    # socioeconomic data
    states.list <- c("US", "AL", "GA", "KY", "MS", "NC", "VA", "TN")
    source("CENSUS_DATA.R")
    pop <- getMultipleStatesPop(states.list)
    df.final <- merge(x=df.final, y=pop, by.x="year", by.y = "year", 
                      all.x = TRUE)
    rm(pop)
    
    source("BEA_DATA.R")
    gdp.us <- getMultipleStatesGDP(states.list, flag.plot = FALSE)
    df.final <- merge(x=df.final, y=gdp.us, by.x="year", by.y = "year", 
                      all.x = TRUE, suffixes=c(".POP", ".GDP"))
    rm(gdp.us)
    
    # read EIA data
    # read Natural gas prices
    source("eiaData.R")
    source("BLS_Data.R")
    ng.prices <- ReadHistoricalData(series.id = "NG.N3010US3.A")
    cpi.annual <- getCPI()
    ng.prices <- merge(x=ng.prices, y=cpi.annual, by.x = "V1", by.y = "year", 
                       all.x=TRUE)
    ng.prices$ng.price.real <- ng.prices$V2 * ng.prices$cpi[nrow(ng.prices)]/ng.prices$cpi
    ng.prices <- ng.prices[ ,c(1, 4)]
    names(ng.prices)[1] <- c("year")
    
    # read total energy from EIA
    series.name <- getEiaHistoricalSeries()[[1]]$series.id
    api.call <- getApiCall(seriesId = series.name)
    myfile <- fromJSON(api.call)
    tot.energy <- as.data.frame(myfile$series$data[[1]], stringsAsFactors = FALSE)
    tot.energy <- as.data.frame(sapply(tot.energy, as.numeric))
    names(tot.energy) <- c("year", "total.energy")
    
    df.final <- merge(x=df.final, y=ng.prices, by="year", all.x = TRUE)
    df.final <- merge(x=df.final, y=tot.energy, by="year", all.x = TRUE)
    df.final$us.energyprod <- df.final$US.GDP/df.final$total.energy  
  }
  
  # just to make sure that order is correct
  df.final <- df.final[order(df.final$time), ]
  
  if (save.to.file) {
    # save final data frame 
    if(name.of.file == "") {
      name.of.file <- "./data/RipsDemand.rds"
    }
    saveRDS(df.final, file = name.of.file) # first as binary
    write.csv(df.final, file = paste0(name.of.file, ".csv")) # then as csv
  }
  
  return(df.final)
}

readDataFromFile <- function(name.of.file = "./data/RipsDemand.rds") {
  # reads previously compiled database from file. 
  #
  # Args:
  #   name.of.file: name of rds file (complete path)
  #
  # Returns:
  #   df.final: data frame with data base
  
  df.final <- readRDS(file=name.of.file)
  
  return(df.final) 
}

compute.calendar.variables <- function(df) {
  
  year.end <- max(year(df$time))
  year.ini <- min(year(df$time))

  # day of week
  df$day.of.week <- weekdays(df$time)
  names.days <- unique(df$day.of.week)
  
  # set holidays
  library("timeDate")
  holi.days <- date(holidayNERC(year = c(year.ini:year.end)))
  df$day.of.week[date(df$time) %in% holi.days] <- "Holiday"
  
  # type of day
  df$type.day <- 
    as.factor(ifelse(df$day.of.week %in% c("Saturday","Sunday", "Holiday"), 
                     "weekend",
                     "weekday")) 
  
  # hour of day
  df$hour.of.day <- as.factor(hour(df$time))
  
  # season
  df$season <- getSeason(df$time)
  df$season <- as.factor(df$season)
  
  #df.final$season2[df.final$season %in% c("Spring", "Fall")] <- "Other"
  #df.final$season2 <- as.factor(df.final$season2)
  
  return(df)
}
