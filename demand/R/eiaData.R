library(jsonlite)
library(reshape2)
library(ggplot2)
library(grid)
library(gridExtra)

library(plyr)
library(dplyr)

# -------------------- UTILITY FUNCTIONS --------------------

getEiaApiKey <- function() {
  # Returns EIA API key code
  #
  # Args:
  #   none
  #
  # returns:
  #   string with API key code
  #
  return("E08FB58BE2C96C433494026A09B6F5A4")
}

getApiCall <- function(seriesId) {
  
  eia.api.key <- getEiaApiKey()
  
  api.call <- paste0("http://api.eia.gov/series/?api_key=",
                     eia.api.key,"&","series_id=", seriesId)
  
  return(api.call)
}

# -------------------- FUNCTIONS TO READ HISTORICAL DATA --------------------

getEiaHistoricalSeries <- function() {
  # Creates and Returns list with relevant EIA historical series
  #
  # Args:
  #   none
  #
  # returns:
  #   list with EIA series
  #
  
  # creates empty series
  seriesList <- list()
  
  # populate list
  # each element of the list is another list with two string elements:
  #   - series.id:  a string with the id of the series that will be used in the
  #                 api call
  #   - descrip:    a string with description of the series
  
  # Total Primary Energy Consumption, Annual
  seriesList[[1]] <- list(series.id="TOTAL.TETCBUS.A", 
                          descrip="Total Primary Energy Consumption, Annual")
  # U.S. Price of Natural Gas Delivered to Residential Consumers, Annual
  seriesList[[2]] <- list(series.id="NG.N3010US3.A", 
                          descrip=paste0("U.S. Price of Natural Gas Delivered",
                                         " to Residential Consumers, Annual"))
  
  return(seriesList)
}

showEiaSeries <- function() {
  # Prints in console list of EIA series defined in function 
  # "getEiaHistoricalSeries"
  #
  # Args:
  #   none
  #
  # returns:
  #   nothing
  #
  
  print((getEiaHistoricalSeries()))
}

ReadHistoricalData <- function(series.id) {
  # Reads historical data from EIA
  #
  # Args:
  #   series.name: name of series
  #
  # returns:
  #   data frame
  #
  
  api.call <- getApiCall(seriesId = series.id)
  myfile <- fromJSON(api.call)
  df.out <- as.data.frame(myfile$series$data[[1]], stringsAsFactors = FALSE)
  df.out <- as.data.frame(sapply(df.out, as.numeric))
  
  return(df.out)
}

ReadTotalEnergyData <- function(print.plot = TRUE, 
                                name.file = NULL,
                                percent = FALSE) {
  # Reads data from TOTAL ENERGY and plots grid by sector
  #
  # Args:
  #   print.plot: (boolean) TRUE if plot should be created
  #   name.file: (string) name of file to save (full path). Default value is
  #               NULL in which case plot is printed in default output
  #   TRUE: (boolean) TRUE if annual share should be plotted.
  #
  # returns:
  #   list with data frames with annual energy consumption by source
  #
  source("ggplotExtra.R")
  
  col.pallete <- c('#8c510a','#bf812d','#dfc27d','#f6e8c3',
                   '#f5f5f5','#c7eae5','#80cdc1','#35978f','#01665e')

  sources <- c("Coal", "Natural Gas", "Petroleum", "Hydroelectric", 
               "Geothermal", "Solar/PV", "Wind", "Biomass",
               "Electricity Retail")
  
  # convert to factors
  sources.f <- factor(sources, levels = sources)
  
  eia.base.url <- "http://www.eia.gov/totalenergy/data/browser/csv.cfm?tbl="
  
  df.tables <- data.frame(residential = "T02.02",
                          comercial = "T02.03",
                          industrial = "T02.04",
                          transportation = "T02.05",
                          stringsAsFactors = FALSE)
  
  list.out <- list()
  
  for (j in 1:4) {
    a <- read.csv(file = paste0(eia.base.url, df.tables[1, j]))
    
    idx.annual.values <- grep("\\d\\d\\d\\d13",a$YYYYMM, perl = TRUE)
    a <- a[idx.annual.values, ]
    
    idx.total.values <- grep("Total",a$Description, perl = TRUE)
    a <- a[-idx.total.values, ]
    
    idx.losses.values <- grep("Losses",a$Description, perl = TRUE)
    a <- a[-idx.losses.values, ]
    
    a$YYYYMM <- as.numeric(substr(as.character(a$YYYYMM), 1, 4))
    names(a)[names(a) == "YYYYMM"] <- "YEAR"
    
    a$Value <- as.numeric(as.character(a$Value))
    a$Unit <- as.character(a$Unit)
    a$Description <- as.character(a$Description)
    
    a$source <- NA
    
    for (i in 1:length(sources)) {
      a$source[grep(sources[i], a$Description, perl = TRUE)] <- sources[i]
    }

    a$source <- factor(a$source, levels = sources)

    # create new data frame that includes all 9 types of sources
    b <- rep(seq(range(a$YEAR)[1], range(a$YEAR)[2]), 
             rep(length(sources), length(unique(a$YEAR))))
    c <- data.frame(YEAR = b, 
                    source = rep(factor(sources, levels = sources), 
                                 length(unique(a$YEAR))))
    d <- merge(a, c, by=c("YEAR", "source"), all.y = TRUE)    
    
    d$Value[is.na(d$Value)] <- 0
    d <- d[order(d$YEAR, d$source), ]
    
    d <- d  %>% group_by(YEAR) %>% 
      select(YEAR, Value, source) %>% 
      mutate(percent = Value/sum(Value, na.rm=TRUE))
    
    list.out[[names(df.tables)[j]]] <- d
  }
  
  list.ggplot <- list()
  
  # sum.eia.data <- final.df  %>% group_by(YEAR, source) %>% 
  #   select(YEAR, Value, source) %>% 
  #   summarize(consumption = mean(Value, na.rm=TRUE))
  # 
  # sum.eia.data <- sum.eia.data %>% group_by(YEAR) %>% 
  #   mutate(percent = consumption/sum(consumption, na.rm = TRUE))
  
  if (print.plot) {
    
    ifelse(percent, plot.variable <- "percent", plot.variable <- "Value")
    
    for (i in 1:4) {
      g <- ggplot(data=list.out[[i]])
      g <- g + geom_area(aes_string(x="YEAR", y=plot.variable, 
                                    fill="source"))
      g <- g + scale_fill_manual(values = col.pallete,
                                 breaks=rev(sources.f))
      g <- g + theme_bw() + guides(fill=guide_legend(title=NULL))
      g <-  g + annotate("text", x=mean(range(list.out[[i]]$YEAR)), 
                         y=Inf, label=names(df.tables)[i], 
                         vjust=1.5, size=4, fontface="bold")
      list.ggplot[[i]] <- g
    }
    
    if (!is.null(name.file)) {
      png(filename = name.file, width=3000, height=3000, res=300)
    }
    #grid.arrange(grobs = list.ggplot)
    grid_arrange_shared_legend(grobs = list.ggplot, 
                               ncol = 2, nrow = 2, position = "right")
    #print(g)
    if (!is.null(name.file)) {
      dev.off()
    }
  }
  
  return(list.out)
}

# -------------------- FUNCTIONS TO READ NEMS PROJECTIONS --------------------

# code to natural gas price
# AEO.2016.REF2016.PRCE_REAL_RES_NA_NG_NA_NA_Y13DLRPMMBTU.A

EiaScenarioNames <- function() {
  # Defines a list with the names and codes of EIA's AEO scenarios to be used
  # with EIA API
  #
  # Args:
  #   none
  #
  # returns:
  #   scenarios: a named list with scenario codes
  #
  
  scenarios <- list()
  
  scenarios$ref.2015 <- "REF2015"
  scenarios$high.econ <- "HIGHMACRO"
  scenarios$low.econ <- "LOWMACRO"
  scenarios$high.oil <- "HIGHPRICE"
  scenarios$low.oil <- "LOWPRICE"
  scenarios$high.oil.gas <- "HIGHRESOURCE"
  scenarios$ref.2014 <- "AEO2014FULL"
  
  return(scenarios)
}


getEiaScenarioCode <- function(name.scenario = "ref.2015") {
  
  scenarios <- EiaScenarioNames()
  
  out <- scenarios[[name.scenario]]
  
  if (is.null(out)) {
    cat("\n ******** Name of Scenario informed not found! **********")
    cat("\n\n Valid names: ", paste(names(scenarios), collapse = ", "), "\n\n")
    cat("\n ********     Returning reference scenario     **********\n")
    out <- scenarios$ref.2015
  }
  return(out)
}

readEiaSercCentralDemand <- function(scenario = "ref.2015",
                                     print.result = FALSE) {
  # Reads SERC Central Electricity Demand Projection from EIA repository
  # (name of series: "Electricity Demand : Total Sales")
  #
  # Args:
  # scenario: string with name of scenario to be imported
  # print.result: boolean. TRUE if result should be printed in screen
  #
  # returns:
  # data.eia: data frame with years and demand
  #
  
  list.scenarios <- EiaScenarioNames()
  code.scenario <- getEiaScenarioCode(scenario) 
  
  if (code.scenario == list.scenarios[[1]] && scenario != "ref.2015") {
    scenario <- "ref.2015"
  }
  
  series.name <- paste0("AEO.2015.", code.scenario, 
                        ".CNSM_NA_ELEP_NA_ELC_NA_SERCCNT_BLNKWH.A")
  
  api.call <- getApiCall(series.name)
  myfile <- fromJSON(api.call)
  
  data.eia <- as.data.frame(myfile$series$data[[1]], stringsAsFactors = FALSE)
  names(data.eia) <- c("YEAR", "VALUE")
  data.eia <- transform(data.eia, YEAR = as.numeric(YEAR),
                        VALUE = as.numeric(VALUE))
  data.eia <- data.eia[order(data.eia$YEAR), ]
  names(data.eia)[2] <- scenario
  row.names(data.eia) <- NULL
  
  if (print.result) {
    cat("\n--------------------------------------------------------\n")
    cat(myfile$series$name, "\n")
    cat(myfile$series$units, "\n")
    print(data.eia)
    cat("--------------------------------------------------------\n")
  }
  
  return(data.eia)
}

plotSercDemand <- function(save.file = FALSE, 
                           name.file = "./out/EiaDemand.png") {
  # Reads and plots Demand scenarios from EIA for SERC Central 
  #
  # Args:
  #   save.file
  #   name.file
  #
  # returns:
  #   nothing
  #
  
  list.scenarios <- EiaScenarioNames()
  n.scenarios <- length(scenarios)
  names.scenarios <- names(list.scenarios)
  
  for (i in 1:n.scenarios) {
    if (i == 1) {
      data.eia <- readEiaSercCentralDemand(scenario = names.scenarios[i])
    } else {
      df.temp <- readEiaSercCentralDemand(scenario = names.scenarios[i])
      data.eia <- cbind(data.eia, df.temp[ ,2])
      names(data.eia)[ncol(data.eia)] <- names(df.temp)[2]
    }
  }
  # convert to long format
  data.eia.long <- melt(a, id.vars = "YEAR", variable.name = "Scenario")
  
  # convert from billion kWh to GWavg
  data.eia.long$value <- data.eia.long$value * 1e9 / (8760 * 1e6)
  
  g <- ggplot(data = data.eia.long)
  g <- g + geom_line(aes(x=YEAR, y=value, colour = Scenario))
  g <- g + geom_point(aes(x=YEAR, y=value, colour = Scenario, 
                          shape = Scenario))
  g <- g + theme_bw() + xlab("Year") + ylab("Demand (GWavg)")
  g <- g + guides(colour=guide_legend(title=NULL), 
                  shape=guide_legend(title=NULL))
  g <- g + theme(legend.position=c(0,1), legend.justification=c(0,1))
  print(g)
  
  if(save.file){
    png(filename = name.file)
    print(g)
    dev.off()
  }
}

# -------------------- FUNCTIONS TO READ 923 FORM --------------------

# jgc <- function(){
#   #http://stackoverflow.com/questions/21937640/handling-java-lang-outofmemoryerror-when-writing-to-excel-from-r
#   .jcall("java/lang/System", method = "gc")
# }
#options(java.parameters = "- Xmx1024m")
#library(xlsx)

download.data.eia <- function (url.file, keep.file = FALSE){
  
  if (keep.file) {
    if(!dir.exists("./data/EIA/")) {
      dir.create(path = "./data/EIA/", recursive = TRUE)
    }
    name.file <- basename(as.character(url.file))
    ff <- paste0("./data/EIA/", name.file)
  } else {
    ff <- tempfile()
  }
  
  # first try to download using method "auto" (default)
  cat("Downloading file ", name.file, "...\n")
  code <- 1
  tryCatch({code <- download.file(url = url.file, destfile = ff)},
           error={function(e){}})
  
  if(code != 0) {
    # error downloading in the first attempt
    # try again with method "curl"
    cat("Error downloading file! Trying with method curl...")
    code <- download.file(url = url.file, destfile = ff, method = "curl", 
                          extra = "-L -#", quiet=TRUE)
    if (code == 0) {
      cat("Download OK!\n\n")
    } else {
      cat("Error downloading file with curl.\n\n")
      stop()
    }
  }
  
  # return path to downloaded file
  return(ff)
}

filter.utility.data <- function(df.data, name.cols = "GEN", 
                                name.utility = "tennessee valley") {
  # Filters data from EIA historical files
  #
  # Args:
  # df.data: (data frame) original EIA data
  # name.cols: (string) name (or part of) of columns to be extracted
  #             (NET ?GEN or GEN)
  # name.utility: (string) name (or part of) of rows to be extracted
  
  # "operator name" "utilname"
  
  cols.gen <- grep(name.cols, names(df.data), perl = TRUE, ignore.case = TRUE)
  if (length(cols.gen) > 12) {
    cols.gen <- head(cols.gen, -1) # drop last column
  }
  
  for(i in cols.gen) {
    df.data[ ,i] <- as.numeric(df.data[ ,i])
  }
  # find column with utility names
  cols.utility <- which(tolower(names(df.data)) == "operator name" | 
                          tolower(names(df.data)) == "utilname")
  
  gen.data <- df.data[ ,c(cols.utility, cols.gen)]
  
  filtered.data <- gen.data[grep(name.utility, gen.data[ ,1], 
                                 perl = TRUE, ignore.case = TRUE), ]
  
  total.monthly.gen <- colSums(filtered.data[ ,-1], na.rm=TRUE)
  
  return(total.monthly.gen)
}

read.form.923 <- function() {
  
  source("excel.R")
  
  library(lubridate)
  library(RCurl)
  library(XML)
  library(ggplot2)
  
  keep.downloaded.file <- TRUE
  
  setwd("~/GoogleDrive/CMU/RIPS/R")
  
  source("NBER.R")
  list.dirs(recursive = FALSE)
  
  name.utility <- "tennessee valley"
  
  # READ FILES 1970 - 2000 ----
  
  years <- 1970:2000
  
  url.eia <- "http://www.eia.gov/electricity/data/eia923/xls/utility/"
  
  base.name <- "f759"
  
  #list.files(pattern = ".xls")
  
  dates.obs <- seq.Date(from=as.Date(paste0(years[1],"-01-01")),
                        to=as.Date(paste0(years[length(years)],"-12-01")),
                        by = "month")
  df.out <- data.frame(date=dates.obs, gen=NA)
  
  for(y in years) {
    #jgc()
    gc()
    if(y < 1996) {
      name.file <- paste0(base.name, y, "u.xls")
    } else {
      name.file <- paste0(base.name, y, "mu.xls")
    }
    
    if (y < 1990) {
      # prior to 1990 the data is in kWh
      # http://www.eia.gov/electricity/data/eia923/xls/F759layout_um.txt
      fct.adj <- 1e-3
    } else {
      fct.adj <- 1
    }
    
    f.csv <- paste0(tempdir(),"/out.csv")
    f.xls <- download.data.eia(url.file=paste0(url.eia, name.file), 
                               keep.file = keep.downloaded.file)
    convert.xls2csv(xls.file = f.xls, csv.file = f.csv)
    if (!keep.downloaded.file) {
      unlink(f.xls)
    }
    
    # get converted csv files
    list.csv <- list.files(path=tempdir(), pattern = "out.csv", 
                           full.names = TRUE)
    
    a <- read.csv(file=list.csv[1])
    unlink(list.csv)
    gc()
    
    total.monthly.gen <- filter.utility.data(df.data = a, name.cols = "GEN",
                                             name.utility = name.utility)
    #adjust from kWh to MWh if year < 1990
    total.monthly.gen <- total.monthly.gen * fct.adj
    
    df.out$gen[year(df.out$date) == y] <- total.monthly.gen
  }
  
  # READ FILES 2001+ ----
  
  # read whole URL to string
  url.eia <- "http://www.eia.gov/electricity/data/eia923/index.html"
  xx <- getURL(url.eia)
  # parse string to XML
  zz <- htmlParse(xx)
  # get XML nodes marked as tbody (tables) and convert to R list
  ww <- getNodeSet(zz, path = "//tbody")
  xxxx <- xmlToList(ww[[1]], simplify = TRUE)
  
  # get only table with links to Excel Files (had to go through list to find)
  xxxx <- xxxx[[3]]$table
  
  # get elements of table that are table rows (marked as "tr")
  irows <- which(names(xxxx) == "tr")
  
  # read table with links to ZIP files
  df.links.zip <- data.frame(year=numeric(), links=as.character())
  for (i in irows) {
    if(!is.null(xxxx[[i]])) {
      if (length(xxxx[[i]]) > 2) {
        yy <- as.numeric(substring(xxxx[[i]][[1]]$text, 1, 4))
        html.link <- xxxx[[i]][[3]]$a$.attrs
        if(yy < 2016) {
          df.links.zip <- rbind(df.links.zip,
                                data.frame(year=yy, 
                                           links=html.link))
        }
      }
    }
  }
  df.links.zip <- df.links.zip[order(df.links.zip$year), ]
  years <- df.links.zip$year
  
  dates.obs <- seq.Date(from=as.Date(paste0(years[1],"-01-01")),
                        to=as.Date(paste0(years[length(years)],"-12-01")),
                        by = "month")
  
  df.out.2 <- data.frame(date=dates.obs, gen=NA)
  url.eia.2 <- "http://www.eia.gov/electricity/data/eia923/"
  
  for (i in 1:nrow(df.links.zip)) {
    gc()
    # download zip file
    f.zip <- download.data.eia(url.file=paste0(url.eia.2, df.links.zip$links[i]), 
                               keep.file = keep.downloaded.file)
    # get name of excel file in zip file (larger file)
    a <- unzip(zipfile = f.zip, list = TRUE)
    name.xls <- a$Name[which(a$Length == max(a$Length))]
    
    #extract excel file and read to R
    f.xls <- unzip(zipfile = f.zip, files = name.xls, exdir = tempdir())
    
    if (!keep.downloaded.file) {
      unlink(f.zip)
    }
    gc()
    
    f.csv <- paste0(tempdir(),"/out.csv")
    convert.xls2csv(xls.file = f.xls, csv.file = f.csv) 
    unlink(f.xls)
    gc()
    
    # get converted csv files
    list.csv <- list.files(path=tempdir(), pattern = "out.csv", 
                           full.names = TRUE)
    
    a.1 <- read.csv(file = list.csv[1], header = FALSE, stringsAsFactors = FALSE)
    unlink(list.csv)
    #a.1 <- read.xlsx2(file=qq, sheetIndex = 1, stringsAsFactors=FALSE)
    
    # find row with headers
    i.header <- grep("plant id", a.1[ ,1], perl = TRUE, ignore.case = TRUE)
    
    # get header
    header <- a.1[i.header, ]
    names(header) <- NULL
    
    #delete rows before header
    a.1 <- a.1[-c(1:(i.header+1)), ]
    names(a.1) <- header
    
    total.monthly.gen <- filter.utility.data(df.data = a.1, 
                                             name.cols = "NET ?GEN",
                                             name.utility = name.utility)
    
    df.out.2$gen[year(df.out.2$date) == df.links.zip$year[i]] <- total.monthly.gen
  }
  
  df.out <- rbind(df.out, df.out.2)
  hours.month <- as.numeric(((df.out$date %m+% months(1)) - df.out$date)*24)
  df.out$gen <- df.out$gen/hours.month
  rm(df.out.2)
  df.out$ma <- stats::filter(df.out$gen, filter=rep(1/12, 12), 
                             method = "convolution", sides = 1)
  
  df.rec <- recessions()
  df.rec <- df.rec[year(df.rec$Peak.month) >= min(year(df.out$date)) ,]
  
  g <- ggplot()
  g <- g + geom_line(data=df.out, aes(x=date, y=gen))
  g <- g + geom_line(data=df.out, aes(x=date, y=ma), colour="red")
  
  g <- g + theme_bw()
  
  g <- g + geom_rect(data=df.rec, aes(xmin=Peak.month, xmax=Trough.month), 
                     ymin=-Inf, ymax=Inf, alpha=.2, fill="darkgray")
  
  if(FALSE) {
    # download weather data from 1970s ----
    source("downloadNOAA.R")
    
    # weather.1 <- readStationFileRemote(USAF=723340, yearini=1970, yearend=2015, 
    #                                    station.isd.data=NULL, 
    #                                    save.to.temp = FALSE)
    list.noaa.files <- list.files(path = "./data/NOAA/", pattern = ".gz", 
                                  full.names = TRUE)
    
    weather.1 <- readStationFileLocal(ls.files = list.noaa.files)
    
    weather.monthly <- weather.1 %>% group_by(year(time), month(time)) %>% 
      summarize(temp=mean(temp, na.rm=TRUE),
                dp=mean(dp, na.rm=TRUE),
                rh=mean(rh, na.rm=TRUE))
    weather.monthly <- as.data.frame(weather.monthly)
    weather.monthly$time <- as.Date(paste(weather.monthly$`year(time)`, 
                                          weather.monthly$`month(time)`, 
                                          1, sep = "-"), 
                                    format = "%Y-%m-%d")
    weather.monthly <- 
      weather.monthly[-which(is.na(weather.monthly$`year(time)`)), ]
    weather.monthly <- weather.monthly[ ,c("time", "temp", "dp", "rh")]
    
    df.out.2 <- merge(x=df.out, y=weather.monthly, 
                      by.x="date", by.y="time", all = TRUE)
    
    plot(df.out.2$temp, df.out.2$gen)
    df.out.2$year <- year(df.out.2$date)
    
    g <- ggplot(data=df.out.2)
    g <- g + geom_point(aes(x=temp,y=gen,colour=year)) + theme_bw()
    g
    
    #states.list <- c("US", "AL", "GA", "KY", "MS", "NC", "VA", "TN")
    source("CENSUS_DATA.R")
    
    source("BEA_DATA.R")
    gdp.us <- readNationalGDP()
  }
}

