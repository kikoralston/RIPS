# BLS_Data

# -------------------- UTILITY FUNCTIONS --------------------

getBlsApiKey <- function() {
  # Returns BLS API key code
  #
  # Args:
  #   none
  #
  # returns:
  #   string with API key code
  #
  return("00f798b763b94cef8a96cfb42af90d25")
}

getBlsApiCall <- function(seriesId) {
  
  bls.api.key <- getBlsApiKey()
  
  api.call <- paste0("http://api.bls.gov/publicAPI/v2/timeseries/data/?api_key=",
                     bls.api.key,"&","series_id=", seriesId)
  
  return(api.call)
}

getCPI <- function() {
  
  library(lubridate)
  library(RCurl)
  library(XML)
  library(plyr)
  library(dplyr)
  
  url.bls <- "http://download.bls.gov/pub/time.series/cu/cu.data.1.AllItems"
  
  a <- read.table(file=url.bls, sep = "\t", header = TRUE, 
                  stringsAsFactors = FALSE)
  a$series_id <- trimws(a$series_id)
  
  a <- a[a$series_id == a$series_id[1], ]
  
  cpi.annual <- a %>% group_by(year) %>% summarise(cpi = mean(value, 
                                                              na.rm=FALSE))
  cpi.annual <- as.data.frame(cpi.annual)
  names(cpi.annual) <- c("year", "cpi")
  return(cpi.annual)
}

if(FALSE){
  a <- ReadHistoricalData(series.id = "NG.N3010US3.A")
  cpi.anual <- getCPI()
  b <- merge(x=a, y=cpi.annual, by.x = "V1", by.y = "year", all.x=TRUE)
  b$ng.price.real <- b$V2 * b$cpi[nrow(b)]/b$cpi
  plot(x=b$V1, y=b$ng.price.real, type="l")
}


