# ***************************************************************************
# Functions to download and read data from the NOAA repository
# ***************************************************************************

library(lubridate)
library(RCurl)

source("auxiliary_functions.R")


getFixColWidths <- function() {
  # returns widths of fixed columns in NOAA's ASCII file
  # (Reference: ftp://ftp.ncdc.noaa.gov/pub/data/noaa/ish-format-document.pdf)
  #
  # Args:
  #   none
  #
  # Returns:
  #   data frame with widths of each field
  
  widths <- data.frame(nv.char=4, usaf=6, wban=5, date=8, time=4, 
                       source=1, lat=6, lon=7, type=5, elev=5, fix.id=5, 
                       qual.name=4, wind.angle=3, wind.qual=1, 
                       w.obs.type=1, wind.speed=4, w.speed.qual=1, 
                       ceil=5, sky.qual=1, sky.code=1, cavok=1, hor.dist=6, 
                       dist.qual=1, dist.var=1, dist.var.qual=1,
                       temp=5, temp.qual=1, dew.point=5, dp.qual=1)
  return(widths)
}

getNoaaFtp <- function() {
  # FTP address of NOAA repository
  #
  # Args:
  #   none
  #
  # Returns:
  #   string with FTP address
  
  return("ftp://ftp.ncdc.noaa.gov/pub/data/noaa/");
}

getNoaaFileHistory <- function() {
  # name of file with list of names of weather stations
  #
  # Args:
  #   none
  #
  # Returns:
  #   string with name of file

  file.history <- 'isd-history.csv'
  ftp.stations.noaa <- getNoaaFtp()
  
  stations.noaa <- read.csv(file = paste0(ftp.stations.noaa, file.history),
                            colClasses = "character")
  return(stations.noaa);
}

listStations <- function(nameCity=NULL, state=NULL, nameCountry = "US") {
  # Connects to NOAA and finds stations that match name of city
  #
  # Args:
  #   nameCity: String with name of city (ignores letter case)
  #   state:    String with name of US state (using the two caracter standard code)
  #   nameCountry: String with code of Country. Defaults to "US". If NULL 
  #                is passed, all countries are included.
  #
  # Returns:
  #   data frame with list of stations
  
  cat("Connecting to NOAA...\n\n")
  cat("--------------------------------- \n\n")
  
  ftp.stations.noaa <- getNoaaFtp()
  stations.noaa <- getNoaaFileHistory()
  
  if (!is.null(nameCountry)) {
    stations.2 <-  stations.noaa[stations.noaa$CTRY == nameCountry, ]
    
    if(nrow(stations.2) == 0) {
      stop(paste0("Zero stations found in country ", nameCountry, "!"))
    }
  } else {
    stations.2 <-  stations.noaa
  }
  
  if (!is.null(nameCity)) {
    idx.stations.name <- grep(tolower(nameCity), 
                              tolower(stations.2$STATION.NAME))
  } else {
    idx.stations.name <- 1:nrow(stations.2)
  }
  
  if (!is.null(state)) {
    idx.stations.state <- 
      which(tolower(state) == tolower(stations.2$STATE))
  } else {
    idx.stations.state <- 1:nrow(stations.2)
  }

  idx.stations <- dplyr::intersect(idx.stations.name,
                                   idx.stations.state)
  
  df.stations <- stations.2[idx.stations,]
  
  return(df.stations)
}

checkStationFile <- function(USAF, WBAN, year) {
  # Connects to NOAA and checks if file for station exists
  #
  # Args:
  #   USAF: String with USAF code for station
  #   WBAN: String with WBAN code for station
  #   year: integer with year of file being searched
  #
  # Returns:
  #   TRUE if file is found. FALSE otherwise
  
  ftp.stations.noaa <- getNoaaFtp()
  
  strcom <- paste0("curl ", paste0(ftp.stations.noaa, year,"/",
                                   paste(USAF, WBAN, year, sep = "-"),
                                   ".gz"), " --ssl --head")
  status <- system(strcom, ignore.stderr = TRUE, ignore.stdout = TRUE)
  
  if (status !=0) {
    return(FALSE)
  } else {
    return(TRUE)
  }
}

downloadStationbyUSAF <- function(USAF, year, path.out="", 
                                  station.isd.data = NULL,
                                  save.to.temp=TRUE) {
  # Connects to NOAA and downloads .gz file with weather station data
  #
  # Args:
  #   USAF: (String)  USAF code for station
  #   year: (integer) year of data
  #   path.out: (String) path where downloaded files should be saved (defaults 
  #             to empty string).
  #   station.isd.data: (String data frame) data frame (same columns 
  #                     as file isd-history.csv) with registry data for  
  #                     station being searched. If this is different from NULL 
  #                     (default value), then USAF argument is not used.
  #   save.to.temp: (boolean) TRUE if downloaded file should be saved to 
  #                 a temporary file
  #
  # Returns:
  #   (string) complete path name of file downloaded. If there is an error
  #   in the download, returns NULL.
  
  ftp.stations.noaa <- getNoaaFtp();
  
  if(is.null(station.isd.data)) {
    cat("Connecting to NOAA...\n\n")
    cat("--------------------------------- \n\n")
    
    # read file with station data
    stations.noaa <- getNoaaFileHistory()
    
    # filter row with data for station being searched
    station.isd.data <- stations.noaa[stations.noaa$USAF == USAF, ]
  }
    
  name.file <- paste(station.isd.data$USAF[1], station.isd.data$WBAN[1], 
                     year, sep = "-")
  
  cat("\n--------------------------------- \n")
  cat("Downloading station ", station.isd.data$STATION.NAME, " for year ",
      year, "\n")
  
  file.out <- NULL
  if (save.to.temp) {
    file.out <- tempfile()
  } else {
    file.out <- paste0(name.file, ".gz")
  }
  if (checkStationFile(station.isd.data$USAF, station.isd.data$WBAN, year)) {
    download.file(url = paste0(ftp.stations.noaa, year,"/", 
                               name.file, ".gz"),
                  destfile = file.out)
  } else {
    cat("File not found in NOAA repository!!\n")
    if (save.to.temp) unlink(file.out)
    file.out <- NULL
  }
  return(file.out)
}

readNoaaFile <- function(name.file, path.file = "") {
  # Reads NOAA gz file and creates data frame with data
  #
  # Args:
  #   name.file: (String)  name of file (local)
  #   path.file: (String) path of file (local)
  #
  # Returns:
  #   data frame with data (time, temperature, dew.point)
  
  col.widths <- getFixColWidths()
  
  b <- read.fwf(file = gzfile(description = paste0(path.file, name.file)), 
                widths = unlist(col.widths), colClasses="character", 
                col.names = names(col.widths))
  hour <- sapply(b$time, substr, start=1, stop=2)
  yy <- as.POSIXct(paste(b$date, paste(hour, "00", sep = ":")),
                   format = "%Y%m%d %H:%M", tz='UTC')
  
  # read temperature
  temp.1 <- as.numeric(b$temp)
  temp.1[temp.1 == 9999] <- NA
  temp.1 <- temp.1/10
  
  # read dew point
  dew.point.1 <- as.numeric(b$dew.point)
  dew.point.1[dew.point.1 == 9999] <- NA
  dew.point.1 <- dew.point.1/10
  
  df.out <- data.frame(time=yy, temp=temp.1, dp=dew.point.1) %>% 
    group_by(time) %>% 
    summarize(temp = mean(temp, na.rm=TRUE),
              dp = mean(dp, na.rm=TRUE))
  
  df.out$rh <- convertDewPoint2RelHum(dew.point = df.out$dp, 
                                      temp = df.out$temp)
  df.out <- as.data.frame(df.out)
  
  # df.out <-  data.frame(time = as.POSIXct(names(temp.hour)),
  #                       temp = as.numeric(temp.hour),
  #                       dew.point = as.numeric(dewpoint.hour),
  #                       rh = rh)
  row.names(df.out) <- NULL
  return(df.out)
}

readStationFileRemote <- function(USAF, yearini, yearend,
                                  station.isd.data=NULL,
                                  save.to.temp=TRUE) {
  # Downloads and Reads data of a NOAA weather station from NOAA's remote 
  # repository
  #
  # Args:
  #   USAF: (String)  USAF code for station
  #   yearini: (integer) initial year of data
  #   yearend: (integer) final year of data
  #   station.isd.data: (String data frame) data frame (same columns as file 
  #                     isd-history.csv) with registry data for station being 
  #                     searched. If this is different from NULL (default 
  #                     value), then USAF argument is not used.
  #   save.to.temp: (boolean) TRUE if downloaded file should be saved to 
  #                 a temporary file
  #
  # Returns:
  #   data frame with data
  
  if(is.null(station.isd.data)) {
    cat("Connecting to NOAA...\n\n")
    cat("--------------------------------- \n\n")
    
    ftp.stations.noaa <- getNoaaFtp();
    # read file with station data
    stations.noaa <- getNoaaFileHistory()
    
    # filter row with data for station being searched
    station.isd.data <- stations.noaa[stations.noaa$USAF == USAF, ]
  }  

  # allocate data frame
  df.out <- data.frame(time = as.POSIXct(character()), 
                       temperature = double())
  
  for (year in yearini:yearend){
    name.file <- downloadStationbyUSAF(USAF = USAF, year = year, 
                                       station.isd.data = station.isd.data,
                                       save.to.temp = save.to.temp)
    #name.file <- paste0(paste(station.isd.data$USAF[1], 
    #                          station.isd.data$WBAN[1], year, sep = "-"),
    #                    ".gz")
    if (!is.null(name.file)) {
      df.aux <- readNoaaFile(name.file = name.file)
      df.out <- rbind(df.out, df.aux)
    }
    if (save.to.temp && !is.null(name.file)) {
      unlink(name.file)
    }
  }
  
  return(df.out)
}

readStationFileLocal <- function(ls.files) {
  # Reads data of a NOAA weather station from previously downloaded gz files
  #
  # Args:
  #   ls.files: (vector of String)  names of gz files
  #
  # Returns:
  #   data frame with data
  
  # allocate data frame
  df.out <- data.frame(time = as.POSIXct(character()), 
                       temperature = double())
  
  for (i in 1:length(ls.files)){
    name.file <- ls.files[i]
    df.aux <- readNoaaFile(name.file = name.file)
    df.out <- rbind(df.out, df.aux)
  }
  
  return(df.out)
}

checkStationFileOld <- function(USAF, WBAN, year) {
  # Connects to NOAA and checks if file for station exists
  # (SLOWER VERSION)
  #
  # Args:
  #   USAF: String with USAF code for station
  #   WBAN: String with WBAN code for station
  #   year: integer with year of file being searched
  #
  # Returns:
  #   TRUE if file is found. FALSE otherwise
  
  cat("Connecting to NOAA...\n\n")
  cat("--------------------------------- \n\n")
  
  ftp.stations.noaa <- getNoaaFtp();
  
  url.content <- getURL(paste0(ftp.stations.noaa, year,"/"), 
                        ftp.use.epsv = FALSE, dirlistonly = TRUE)
  df.files <- as.data.frame(strsplit(url.content, split = "\n"), 
                            stringsAsFactors = FALSE)
  names(df.files) <- c("ftp.file")
  
  # get column with file names
  df.files$ftp.file <- (substring(df.files$ftp.file,
                                  nchar(df.files$ftp.file)-19,
                                  nchar(df.files$ftp.file)))
  list.stations <- paste0(paste(USAF, WBAN, year, sep = "-"), ".gz")
  flag.exists <- (list.stations %in% df.files$ftp.file)
  return(flag.exists)
}

compileDataListStations <- function(file.list, year.ini, year.end, 
                                    save.csv=FALSE) {
  # Compiles a data frame with weather data for a list of stations
  #
  # Args:
  #   file.list:  (string) complete path of csv file with station info.
  #               (USAF, WBAN, station name, country, state)
  #   yearini: (integer) initial year of data
  #   yearend: (integer) final year of data
  #
  # Returns:
  #   data frame with weather data for each station

  listOfStations <- read.csv(file=file.list, colClasses = "character",
                             comment.char = "#")
  ls.data <- list()
  for (i in 1:nrow(listOfStations)) {
    ls.data[[i]] <- 
      readStationFileRemote(USAF=NULL, yearini = year.ini, yearend = year.end, 
                            station.isd.data = listOfStations[i, ])
    full.name <- listOfStations$STATION.NAME[i]
    name.station <- substr(full.name, 1, (regexpr(" ", full.name, 
                                                  perl = TRUE)-1))
    
    names(ls.data[[i]]) <- c("time", paste(names(ls.data[[i]])[-1],
                                           name.station,
                                           sep = "."))
  }
  
  if(nrow(listOfStations) == 1) {
    final.df <- ls.data[[1]]
    final.df <- final.df %>% dplyr::filter(!(is.na(time)))
  } else {
    for (i in 2:length(ls.data)) {
      if (i == 2) {
        final.df <- dplyr::full_join(ls.data[[1]], ls.data[[2]], 
                          by="time")
      } else {
        final.df <- dplyr::full_join(final.df, ls.data[[i]], 
                          by="time")
      }
      final.df <- final.df %>% dplyr::filter(!(is.na(time)))
    }
  }
  
  
  weather.data <- final.df
  idx.rh <- grep("rh.", names(weather.data))
  weather.data[ ,idx.rh] <- round(weather.data[ ,idx.rh], 2)
  
  save(weather.data, 
       file = "./data/weather_data.Rdata") # first as binary
  if (save.csv) {
    write.csv(weather.data, file = "./data/weather_data.csv",
              row.names = FALSE) # then as csv
  }
  
  return(weather.data)
}