# ************************************************************************
#           functions to download relevant BEA data
# ************************************************************************


#library(data.table)
library(RCurl)
library(gridExtra)
library(jsonlite)
library(plyr)
library(dplyr)
library(tidyr)
library(ggplot2)

getBEAFipsCodes <- function() {
  url.census <- paste0("http://www2.census.gov/geo/docs/reference/codes/",
                       "files/national_county.txt")
  df.1 <- read.csv(file=url.census, colClasses = "character")
  df.1 <- unique(df.1[ ,c(1:2)])
  names(df.1) <- c("state", "geofips")
  df.1$geofips <- paste0(df.1$geofips, "000")
  row.names(df.1) <- NULL
  return(df.1)
}

getBEAKey <- function() {
  # Returns API key to access BEA Data
  # (email from 2016-09-17)
  #
  # Args:
  #
  # Returns:
  #   (string) API key
  
  api.key <- "6E8E3542-05A0-48A2-A606-003BEF0CC9D5"
  
  return(api.key)
}

getIndustryId <- function() {
  # NAIC
  industryId <- list(agric=3, mine=6, util=10, constr=11, manuf=12, whole=34, 
                     retail=35, transp=36, infor=45, finan=50, prof=59, edu=68, 
                     art=74, other=81, gov=82)
  return(industryId)
}

getIndustrySector <- function() {
  industrySector <- list(agric=1, mine=1, util=2, constr=2, manuf=2, whole=3, 
                         retail=3, transp=3, infor=3, finan=3, other=3, gov=3)
  return(industrySector)
}

getIndustryIdSIC <- function() {
  # SIC
  industryId <- list(agric=3, mine=6, util=46, constr=11, manuf=12, whole=47, 
                     retail=48, transp=37, infor=45, finan=49, other=57, 
                     gov=72)
  return(industryId)
}


mapNaicSic <- function(naicCode) {
  # maps to industry code from NAIC methodology to SIC methodology
  #
  # Arg:
  
  map.out <- data.frame(naic = c(1, 3, 6, 10, 11, 12, 34, 35, 36, 45, 
                                 50, 59, 68, 74, 81, 82),
                        sic =  c(1, 3, 6, 46, 11, 12, 47, 48, 37, 45, 
                                 49, 57, 57, 57, 57, 72))
  return(map.out)
}

getNationalData <- function(tableId=3, frequency="A", years = "ALL") {
  # Get national data from BEA
  # http://www.bea.gov/API/bea_web_service_api_user_guide.htm (Appendix B)
  #
  # Args:
  # tableId:  (integer) BEA code for table to be imported
  #             Default: 3 (Table 1.1.3. Real GDP)
  # frequency: (string) frequency of data. "A" (annual), "Q" (quarterly),
  #             "M" (monthly). Default: "A"
  # years: (string or numeric vector) years of data. Default: ALL
  
  if (years == "ALL") {
    years <- "X"
  } else {
    # Concatenate year vector 
    years <- paste(years, collapse = ",")
  }
  
  key <- getBEAKey()
  dataset <- "NIPA"
  
  # Concatenate URL
  this_request <- paste0(
    "http://www.bea.gov/api/data?&UserID=", key,
    "&method=GetData&DataSetName=", dataset,
    "&TableID=", tableId,
    "&Frequency=", frequency,
    "&Year=", years,
    "&ResultFormat=JSON"
  )
  
  # Pull data, converting from JSON to an R list object
  data_json <- fromJSON(this_request)
  
  # Pick out dataframe
  data_frame <- data_json$BEAAPI$Results$Data
  types.cols <- data_json$BEAAPI$Results$Dimensions
  
  # convert to numeric
  name.value <- types.cols$Name[types.cols$IsValue==1]
  data_frame[,name.value] <- gsub(",", "", data_frame[,name.value])
  data_frame[,name.value] <- as.numeric(data_frame[,name.value])
  
  return(data_frame)
}

getRegionalData <- function(component="RGDP_SAN", industryId="1", 
                            geoFips="00000", years = "ALL") {
  # Get regional data from BEA
  # http://johnricco.github.io/2016/07/21/bea-api-r/
  #
  # Args:
  # component:  (string) BEA code for type of data to be fetched. 
  #             Default: RGDP_SAN (Real GDP by State)
  # industryId: (string) NAIC numeric code for Industry Id. 
  #             Default: 1 (all industries) 
  # geoFips:  (String) numeric code for GeoFips code (location). 
  #           Default: 00000 (US)
  # years: (string or numeric vector) years of data. Default: ALL
  
  # Concatenate year vector 
  if (years != "ALL") {
    years <- paste(years, collapse = ",")
  }
  
  key <- getBEAKey()
  dataset="RegionalProduct"
  
  # Concatenate URL
  this_request <- paste0(
    "http://www.bea.gov/api/data?&UserID=", key,
    "&method=GetData&DataSetName=", dataset,
    "&Component=", component,
    "&IndustryId=", industryId,
    "&GeoFips=", geoFips,
    "&Year=", years,
    "&ResultFormat=JSON"
  )
  
  # Pull data, converting from JSON to an R list object
  data_json <- fromJSON(this_request)
  
  # Pick out dataframe
  data_frame <- data_json$BEAAPI$Results$Data
  types.cols <- data_json$BEAAPI$Results$Dimensions
  
  # convert to numeric
  name.value <- types.cols$Name[types.cols$IsValue==1]
  data_frame[,name.value] <- as.numeric(data_frame[,name.value])
  
  return(data_frame)
}

getStateGDP <- function(state.name = "US", allInd=TRUE) {
  
  if (allInd) {
    indId <- list(all=1)
  } else {
    indId <- getIndustryId()
  }
  
  if (state.name == "US") {
    code.state <- "00000"
  } else { 
    code.state <- getBEAFipsCodes() %>% dplyr::filter(state==state.name) %>% 
      dplyr::select(geofips)
    code.state <- as.character(code.state)
  }
  
  #NAICS
  for(nameInd in names(indId)) {
    ll <- getRegionalData(component="RGDP_SAN",
                          industryId = indId[[nameInd]], 
                          geoFips = code.state)
    
    if(nameInd == names(indId)[1]) {
      naics <- ll[ ,c("TimePeriod", "DataValue")]
      names(naics) <- c("year", nameInd)
    }else{
      naics <- data.frame(naics, ll$DataValue)
      names(naics)[ncol(naics)] <- nameInd
    }
  }
  
  if (!allInd) {
    # aggregate columns in other
    cols.other <- c("prof", "edu", "art", "other")
    other <- rowSums(naics[ ,cols.other])
    naics <- naics[ , !(names(naics) %in% cols.other)]
    naics <- data.frame(naics, other)

    indId.2 <- getIndustryIdSIC()
  } else {
    indId.2 <- list(all=1)
  }
  
  #SIC
  for(nameInd in names(indId.2)) {
    ll <- getRegionalData(component="RGDP_SAS",
                          industryId = indId.2[[nameInd]], 
                          geoFips = code.state)
    
    if(nameInd == names(indId.2)[1]) {
      sic <- ll[ ,c("TimePeriod", "DataValue")]
      names(sic) <- c("year", nameInd)
    }else{
      sic <- data.frame(sic, ll$DataValue)
      names(sic)[ncol(sic)] <- nameInd
    }
  }
  
  #make sure order of columns is the same
  if (!allInd) naics <- naics[ ,names(sic)]
  nr <- nrow(sic)
  sic.base <- sic[nr, ]
  for (i in 2:ncol(sic)) {
    sic[ ,i] <- sic[ ,i]/unlist(sic.base[i])
  }
  for (i in 2:ncol(sic)) {
    sic[ ,i] <- sic[ ,i]*unlist(naics[1, i])
  }
  sic <- sic[-nr, ]
  
  gdp <- rbind(sic,naics)
  
  return(gdp)
}

getMultipleStatesGDP <- function(states.list, flag.plot = FALSE) {
  
  cat("\n------------------------------\n")
  cat("Reading GDP data from BEA ... \n")
  
  for (st in states.list) {
    if(st == states.list[1]) {
      df <- getStateGDP(state.name = st, allInd = TRUE)
      names(df)[names(df) == "all"] <- st
    } else {
      df1 <- getStateGDP(state.name = st, allInd = TRUE)
      names(df1)[names(df1) == "all"] <- st
      df <- merge(x=df, y=df1, all=TRUE)
    }
  }

  cat("Done! \n")
  cat("------------------------------\n")
  
  if (flag.plot) {
    df.long <- gather_(df, key_col="state", value_col="gdp", 
                       gather_cols = states.list, factor_key = TRUE)
    df.long$year <- as.numeric(df.long$year)
    g <- ggplot(df.long)
    g <- g + geom_line(aes(x=year, y=gdp, colour=state))
    g <- g + geom_point(aes(x=year, y=gdp, shape=state))
    g <- g + xlab("Year")
    g <- g + ylab("GDP")
    print(g)
  }
  
  return(df)
}


plotStateGDPSector <- function(state) {
  state.data <- getStateGDP(state.name = state, 
                            allInd=FALSE)
  sectors <- getIndustrySector()
  
  primary <- rowSums(state.data[ ,names(sectors[sectors == 1])])
  secondary <- rowSums(state.data[ ,names(sectors[sectors == 2])])
  tertiary <- rowSums(state.data[ ,names(sectors[sectors == 3])])
  
  total <- primary + secondary + tertiary
  
  primary.perc <- primary/total
  secondary.perc <- secondary/total
  tertiary.perc <- tertiary/total
  
  gdp.sector <- data.frame(year=as.numeric(state.data$year),
                           primary, secondary, tertiary)
  gdp.sector.2 <- data.frame(year=as.numeric(state.data$year),
                             primary=primary.perc, 
                             secondary=secondary.perc, 
                             tertiary=tertiary.perc)
  
  gdp.sector <- gdp.sector[complete.cases(gdp.sector), ]
  gdp.sector.2 <- gdp.sector.2[complete.cases(gdp.sector.2), ]
  
  gdp.sector <- gather(gdp.sector, key="sector", 
                       value="value", primary:tertiary,
                       factor_key = TRUE)
  gdp.sector.2 <- gather(gdp.sector.2, key="sector", 
                       value="value", primary:tertiary,
                       factor_key = TRUE)
  
  g <- ggplot(data=gdp.sector)
  g <- g + geom_line(aes(x=year, y=value, colour=sector)) + 
    geom_point(aes(x=year, y=value, shape=sector))
  print(g)
  g <- ggplot(data=gdp.sector.2)
  g <- g + geom_line(aes(x=year, y=value, colour=sector)) + 
    geom_point(aes(x=year, y=value, shape=sector))
  print(g)
}

readFileBEANationalData <- function() {
  # Reads csv file from BEA with US GDP data
  #
  # Args:
  #
  # Returns:
  #   vector of strings where each element is a line of the csv file
  
  address.bea.1 <-
    "http://www.bea.gov//national/nipaweb/SS_Data/Section1All_csv.csv"
  
  conn <- url(address.bea.1)
  
  orig.data <- readLines(conn)
  close(conn)
  
  return(orig.data)
}

readNationalGDP <- function() {
  # Reads BEA data with US GDP data and filters Total Real GDP data
  # (also converts strings to numbers)
  #
  # Args:
  #
  # Returns:
  #   data frame with YEAR and GDP.US values
  
  # reads CSV file form BEA and allocates it to array
  orig.data <- readFileBEANationalData()
  
  # get line with National Real GDP (code A191RX1)
  irow.gdp.real <- grep("A191RX1", orig.data)[1]
  
  year.row <- orig.data[irow.gdp.real-1]
  data.row <- orig.data[irow.gdp.real]
  
  table.1 <- read.table(textConnection(year.row), sep=",",
                        stringsAsFactors = FALSE)
  table.2 <- read.table(textConnection(data.row), sep=",",
                        stringsAsFactors = FALSE)
  table.1[1] <- NA
  
  table.final <- rbind(table.1, table.2)
  col.data <- !(is.na(table.1))
  table.final <- table.final[ ,col.data]
  table.final <- as.data.frame(t(table.final))
  
  table.final[ ,1] <- as.numeric(as.character(table.final[ ,1]))
  table.final[ ,2] <- as.character(table.final[ ,2])
  
  # erase thousand separator and convert to numeric
  table.final[ ,2] <- as.numeric(gsub(",", "", table.final[ ,2]))
  
  # convert from Billion to Million (State GDP is in million)
  table.final[ ,2] <- table.final[ ,2] * 1000
  
  names(table.final) <- c("year", "GDP.US")
  row.names(table.final) <- NULL
  
  return(table.final)
}