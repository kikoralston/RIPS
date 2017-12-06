# ************************************************************************
#           functions to download CENSUS data
# ************************************************************************

# List with all datasets in Census API
# http://api.census.gov/data.html

library(RCurl)
library(XML)
library(reshape2)
library(ggplot2)
library(jsonlite)
library(plyr)
library(dplyr)
library(lubridate)

get.links.county.by.state <- function(){
  # download html
  html <- getURL('https://www.census.gov/geo/reference/codes/cou.html')
  
  # parse string to XML
  doc <- htmlParse(html)
  
  # get XML nodes
  states <- xpathApply(doc, 
                       path = "//optgroup[@label=\"States\"]/option", 
                       xmlValue)
  link.states <- xpathApply(doc, 
                            path = "//optgroup[@label=\"States\"]/option", 
                            xmlGetAttr, "value")
  
  df.list.states.links <- data.frame(state=unlist(states),
                                     link=unlist(link.states))
  
  return(df.list.states.links)
}

get.df.counties.state <- function(name.state){
  
  df <- get.links.county.by.state()
  
  if(!(tolower(name.state) %in% tolower(df$state))) {
    stop('State \'', name.state, "\' was not found in CENSUS. Check name.")
  }
  
  url.state <- df %>% filter(tolower(state) == tolower(name.state)) %>%
    select(link) %>% unlist() %>% as.character()
  
  ff <- tempfile()
  download.file(url = url.state, destfile = ff, quiet = TRUE)
  t <- read.table(file=ff, sep = ',', stringsAsFactors = FALSE)
  unlink(ff)
  t <- t %>% select(V1, V4)
  names(t) <- c('region', 'subregion')
  t$subregion <- tolower(t$subregion)
  t$subregion <- trimws(gsub('county', '', t$subregion, ignore.case = TRUE))
  t$subregion <- trimws(gsub('\\.', '', t$subregion, ignore.case = TRUE))
  return(t)
}


getStateList <- function() {
  url.census <- 'https://www.census.gov/geo/reference/ansi_statetables.html'
  a <- getURL(url.census)
  b <- htmlParse(a)
  # get XML nodes marked as tbody (tables) and convert to R list
  ww <- getNodeSet(b, path = "//tbody")
  
  # get first table
  www <- readHTMLTable(doc = ww[[1]], 
                       header = FALSE,
                       stringsAsFactors = FALSE)
  
  names(www) <- c('Name', 'FIPS', 'USPS')
  return(www)
}

getStateFIPS <- function(usps.state) {
  df.1 <- getStateList()
  state.fips <- df.1$FIPS[tolower(df.1$USPS) %in% tolower(usps.state)]
  return(state.fips)
}

getStateName <- function(fips.state) {
  df.1 <- getStateList()
  state.name <- df.1$Name[tolower(df.1$FIPS) %in% tolower(fips.state)]
  return(state.name)
}

getCensusKey <- function() {
  # Returns API key to access Census Data
  # (email from 2016-09-18)
  #
  # Args:
  #
  # Returns:
  #   (string) API key
  
  api.key <- "91716cc7f14d8f9b1603ccc7d392471e8f40ef8f"
  
  return(api.key)
}

getApiPrefixState <- function() {
  
  # https://www.census.gov/data/developers/data-sets/popest-popproj/popest.2000-2010_Intercensals.html
  
  list.apis <- list()
  
  #list of apis with annual intercensus data for each decade
  list.apis$annual.pop.1990 <- 
    paste0("https://api.census.gov/data/1990/pep/int_charagegroups?", 
           "get=AGEGRP,RACE_SEX,YEAR,HISP,POP&for=county:*&in=state:")
  list.apis$annual.pop.2000 <- 
    "https://api.census.gov/data/2000/pep/int_population?get=DATE,POP&for=state:"
  list.apis$annual.pop.2010 <- 
    "https://api.census.gov/data/2015/pep/population?get=DATE,POP&for=state:"
  
  return(list.apis)
}

getApiPrefixUS <- function() {
  
  list.apis <- list()
  
  # APIS for national population
  list.apis$annual.us.pop.2010 <- 
    "https://api.census.gov/data/2015/pep/population?get=DATE,POP&for=us:1"
  list.apis$annual.us.pop.2000 <- 
    "https://api.census.gov/data/2000/pep/int_population?get=DATE,POP&for=us:1"
  list.apis$annual.us.pop.1990 <- 
    "https://api.census.gov/data/1990/pep/int_natrespop?get=YEAR,TOT_POP"
  
  return(list.apis)
}

getApiPrefix <- function() {
  list.apis <- c(getApiPrefixUS(), getApiPrefixState())
  return(list.apis)
}

filterPopulationData <- function(df, type.data) {
  
  # intercensus data has some adjustments to the base data
  # 1) census data are reported in April 1st (first element)
  # 2) then they adjust it to July (second element) because all other annual 
  # values are referecend to July (this is the value we will use)
  
  # the second element of the postcensus (current decade) has data
  # a adjusted data value for April 1st. This will also be igored
  
  # For each decade I will only consider form year 0 to year 9
  # (10 values)
  
  # the state data from 1990 decade already includes only the values we need
  # (years 0 to 9).

  nc <- nchar(type.data)
  type.year <- substr(type.data, nc-3, nc)
  
  if(type.year == "2000") {
    # intercensus 2000
    df <- df[c(-1), ]
    df <- df[c(1:10), ]
    names(df)[names(df) == "DATE"] <- "year"
    df$year <- c(2000:2009)
  } else { 
    if (type.year == "2010") {
      # post census 2010
      names(df)[names(df) == "DATE"] <- "year"
      df <- df[c(-1, -2), ]
      df$year <- c(2010:2015)
    } else { 
      # intercensus 1990
      names(df)[names(df) == "YEAR"] <- "year"
      df <- df[(df$year %in% c(90:99) | df$year %in% c(1990:1999)), ]
      df$year <- c(1990:1999)
    }
  }
  
  #valid column headers
  val.col <- c("year", "POP", "TOT_POP")
  df <- df[ ,names(df) %in% val.col]
  
  #change TOT_POP to POP
  if ("TOT_POP" %in% names(df)) {
    names(df)[names(df) == "TOT_POP"] <- "POP"
  }
  
  return(df)
}

getAnnualPopData <- function(state="47", type.data="annual.pop.2000"){
  
  if (!(type.data %in% names(getApiPrefix()))) {
    cat("\nArgument \'type.data\' not valid!\n\nMust be one of the following:\n")
    cat(paste(names(getApiPrefix()), collapse = "\n"))
    stop()
  }
  
  key <- getCensusKey()
  
  if(state == "00") {
    prefix <- getApiPrefixUS()[[type.data]]
    this_request <- paste0(prefix, "&key=", key)
  } else {
    prefix <- getApiPrefixState()[[type.data]]
    this_request <- paste0(prefix, state, "&key=", key)
  }
  
  # read data from Census
  a <- getURL(this_request)
  
  if(substr(a, 1, 6) == "html") {
    cat("\nError reading census data.")
    stop()
  }
  
  if(a == "") {
    # nothing found. return NULL
    return(NULL)
  }
  
  # format of data is a python's list style
  # convert to data frame
  # delete "[" and "]"
  b <- gsub("[\\[\\]]", "", a, perl = TRUE)
  conn <- textConnection(b)
  dd <- read.csv(file=conn)
  close.connection(conn)
  
  nc <- nchar(type.data)
  if(state != "00" && substr(type.data, nc-3, nc) == "1990") {
    # the 1990 data set only exists divided by county, agre group, race etc.
    # here we sum it all up to get the total by state
    ee <- dd %>% group_by(YEAR, state) %>% 
      select(YEAR, state, POP) %>% summarize(POP = sum(POP, na.rm=TRUE))
    dd <- as.data.frame(ee)
  }
  
  dd <- filterPopulationData(dd, type.data)

  # return data frame
  return(dd)
}

combineAllCensusData <- function(state="01", flag.plot = FALSE) {
  
  if(state == "00") {
    data.sets <- names(getApiPrefixUS())
  } else {
    data.sets <- names(getApiPrefixState())
  }
  
  df.out <- NULL
  cat("\n-------------------------------------\n")
  for (i in 1:length(data.sets)) {
    cat("Reading ", data.sets[[i]], "for state: ", state,"... ")
    df.1 <- getAnnualPopData(type.data = data.sets[[i]], state = state)
    if (is.null(df.1)) {
      cat("Data not found!\n")
    } else {
      if(i == 1) {
        df.out <- df.1
      } else {
        df.out <- rbind(df.out, df.1)
      }
      cat("Done!\n")
    }
  }
  cat("-------------------------------------\n")
  #make sure order of years is correct
  
  if(!is.null(df.out)) {
    df.out <- df.out[order(df.out$year), ]
    if (flag.plot) {
      g <- ggplot(df.out)
      g <- g + geom_line(aes(x=year, y=POP/1e3))
      g <- g + ggtitle(getStateName(state))
      g <- g + xlab("Year")
      g <- g + ylab("Population (thousands)")
      print(g)
    }
  }
  return(df.out)
}

getMultipleStatesPop <- function(states.list, flag.plot = FALSE) {
  
  list.codes <- getStateFIPS(states.list)
  names(list.codes) <- getStateName(list.codes)
  
  df <- NULL
  for (st in list.codes) {
    df1 <- combineAllCensusData(state = st, flag.plot = FALSE)
    if (!is.null(df1)) {
      df1 <- df1[,c("year", "POP")]
      names(df1)[names(df1) == "POP"] <- names(list.codes)[list.codes==st]
      if(st == list.codes[1]) {
        df <- df1
      } else {
        df <- merge(x=df, y=df1, all=TRUE)
      }
    }
  }
  
  if (flag.plot) {
    df.long <- tidyr::gather_(df, key_col="state", value_col="pop", 
                              gather_cols = states.list, factor_key = TRUE)
    df.long$DATE <- as.numeric(df.long$DATE)
    g <- ggplot(df.long)
    g <- g + geom_line(aes(x=DATE, y=pop, colour=state))
    g <- g + geom_point(aes(x=DATE, y=pop, shape=state))
    g <- g + xlab("Year")
    g <- g + ylab("GDP")
    print(g)
  }
  
  return(df)
}

readNationalPop <- function() {
  library(lubridate)
  url.census <- "https://www.census.gov/popest/data/national/totals/pre-1980/tables/popclockest.txt"

  a <- readLines(url.census)
  first.col <- sapply(a, FUN = substr, 1, 13, USE.NAMES = FALSE)
  first.col <- trimws(first.col)
  first.col <- as.Date(first.col, format = "%B %d, %Y")
  
  ir.data <- which(!is.na(first.col))
  first.col <- first.col[ir.data] 
  
  a <- a[ir.data]
  
  pop <- sapply(a, 
                FUN = function(x) as.numeric(gsub(",", "", 
                                                  trimws(substr(x, 17, 32)))),
                USE.NAMES = FALSE)
  
  df.pop <- data.frame(year = year(first.col), pop=pop)
  df.pop <- df.pop[order(df.pop$year), ]
  
  df.pop.2 <- combineAllCensusData(state = "00")
  
  return(df.pop)
}

readHistoricalStatePop <- function() {
  url.census <- paste0("https://www.census.gov/popest/data/state/asrh/1980s/", 
                       "80s_st_totals.html")
  xx <- getURL(url.census)
  # parse string to XML
  zz <- htmlParse(xx)
  # get XML nodes marked as tbody (tables) and convert to R list
  ww <- getNodeSet(zz, path = "//tbody")
  xxxx <- xmlToList(ww[[1]], simplify = TRUE)
  
  links.census <- data.frame(year.ini=numeric(), year.end=numeric(),
                             url=character())
  # get links to files
  for (i in seq(1, length(xxxx), by=2)){
    # data of table is only in odd indexes
    years.desc <- unlist(strsplit(xxxx[[i]]$a$text, split = "-"))
    year.ini <- as.numeric(years.desc[1])
    year.end <- as.numeric(paste0(substr(years.desc[1], 1, 2), years.desc[2]))
    links.census <- rbind(links.census,
                          data.frame(year.ini=year.ini, 
                                     year.end=year.end,
                                     url=xxxx[[i]]$a$.attrs[2]))
  }
  
  list.data <- list()
  for (j in 1:nrow(links.census)) {
    url.census <- paste0("https://www.census.gov", links.census$url[j])
    conn.url <- url(description=url.census)
    a <- readLines(con=conn.url)
    close(conn.url)
    names.rows <- sapply(a, FUN=substring, 1, 2, USE.NAMES = FALSE)
    names.states <- getStateList()[ ,1]
    codes.states.num <- as.numeric(getStateList()[ ,2])
    
    # check if columns FIPS exist in header and return row of header
    fips.header <- grep("\\bFIPS?\\b", a, perl = TRUE, ignore.case = TRUE)
    if(length(fips.header) > 0) {
      idx.header <- unlist(strsplit(a[fips.header], split = " ")) != ""
      header <- unlist(strsplit(a[fips.header], split = " "))[idx.header]
      pos.fips <- grep("\\bFIPS?\\b", header, perl = TRUE, ignore.case = TRUE)
    }
    
    # data for some decades is in thousands
    # look for the word "thousands" in description and create conversion factor
    thousands.flag <- length(grep("thousands", a, perl = TRUE, 
                                  ignore.case = TRUE)) > 0
    if (thousands.flag) {
      conv.factor <- 1e3
    } else {
      conv.factor <- 1
    }
    
    firstletter <- paste0(unique(sapply(names.states, substring, 1, 1, 
                                        USE.NAMES = FALSE)), 
                          collapse = "")
    secondletter <- paste0(unique(sapply(names.states, substring, 2, 2, 
                                         USE.NAMES = FALSE)), 
                           collapse = "")
    reg.class <- paste0("\\b[",firstletter,"]","[",secondletter,"]\\b")
    
    irows <- grep(reg.class, a, perl = TRUE)
    
    a <- a[irows]
    b <- sapply(a, 
                FUN = function(x){
                  unlist(strsplit(x, " "))[unlist(strsplit(x, " ")) !=""]}, 
                USE.NAMES = FALSE)
    
    # delete from with US data
    us.row <- grep("\\bUS\\b", b, perl = TRUE, ignore.case = TRUE)
    
    if(length(us.row) > 0) {
      b <- b[-us.row]
    }
    
    # delete thousands separator (comma)
    b <- sapply(b, FUN=function(y){gsub(",", "", y)}, 
                simplify = FALSE)
    
    n.list <- length(b)
    
    for (i in 1:(n.list/2)) {
      tt.2 <- as.numeric(b[[i+(n.list/2)]])
      if(length(fips.header) > 0) {
        # if there is a column with FIPS
        tt.2 <- tt.2[!(is.na(tt.2) | tt.2 %in% codes.states.num)]
      } else {
        tt.2 <- tt.2[!(is.na(tt.2))]
      }
      
      b[[i]] <- c(b[[i]], tt.2)
    }
    b <- b[1:(n.list/2)]
    c <- as.data.frame(matrix(data=unlist(b), 
                              nrow=length(b), 
                              ncol=length(b[[1]]), 
                              byrow = TRUE), stringsAsFactors = FALSE)
    if(length(fips.header) > 0) {
      # if there is a column with FIPS, delete it
      c <- c[ ,-1]
    }
    
    for (icol in 2:ncol(c)) {
      # convert population values to numeric
      c[ ,icol] <- as.numeric(c[ ,icol])
    }
    
    # multiply conversion factor
    c[ ,2:ncol(c)] <- c[ ,2:ncol(c)] * conv.factor
    
    # order states alphabeticaly
    c <- c[order(c[ ,1]), ]
    
    # filter extra and repeated years
    if (links.census$year.end[j] %% 10 == 0) {
     # if last year is XXX0, delete last year
     c <- c[ ,-ncol(c)]
    }
    
    if (ncol(c) > 11) {
      # if there are still more than 10 years, first year is duplicated
      # delete first point (second column)
      c <- c[ ,-2]
    }

    # initial.decade <- as.numeric(paste0("19", decades[j]))
    names(c) <- c("state", seq(from=links.census$year.ini[j], 
                               to=(links.census$year.ini[j]+ncol(c)-2)))
    
    #put in long format
    c <- tidyr::gather(c, "year", "population", 2:ncol(c))
    
    list.data[[j]] <- c
  }
  
  for (i in 1:length(list.data)) {
    if (i == 1) {
      complete.pop <- list.data[[1]]
    } else {
      complete.pop <- rbind(complete.pop, list.data[[i]])
    }
    #complete.pop <- cbind(complete.pop, list.data[[i]][ ,-1])
  }
  
  complete.pop <- complete.pop[order(complete.pop$year, complete.pop$state), ]
  complete.pop$state <- as.factor(complete.pop$state)
  complete.pop$year <- as.numeric(complete.pop$year)
  
  # read data from 1990+
  names.states <- names.states[-which(names.states == "US")]
  states.data.1990 <- getMultipleStatesPop(names.states)
  
  states.data.1990 <- tidyr::gather(states.data.1990, "state", "population", 
                                    2:ncol(states.data.1990))
  states.data.1990 <- states.data.1990[ ,c("state", "year", "population")]
  states.data.1990$state <- as.factor(states.data.1990$state)
  
  complete.pop <- rbind(complete.pop, states.data.1990)
  
  library(ggplot2)
  library(plotly)
    
  g <- ggplot(data=complete.pop)
  g <- g + geom_line(aes(x=year, y=population, colour=state))
  g <- g + geom_point(aes(x=year, y=population, bg=state), shape=21)
  g <- g + scale_colour_discrete()
  
  ggplotly(g)
}


