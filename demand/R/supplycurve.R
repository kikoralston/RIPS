source("excel.R")

library(plyr)
library(dplyr)
library(tidyr)
library(XML)
library(RCurl)


fill.up.costs.thermal <- function(df) {
  # this function fills up missing data on heat rate and fuel costs 
  # for thermal plants in df
  #
  # rule:
  # heat rates will be filled with national avarages from EIA (Table 8.1)
  #
  # fuel will be filled with state average for selected fuel
  
  average.data <- read.eia.fuel.costs.states()
  national.average <- average.data %>% filter(State == 'U.S. Total')
  national.average$State <- NULL
  
  # use default sets from R with states info
  df.states <- data.frame(state.name, state.abb, stringsAsFactors = FALSE)

  # join data frames
  average.data <- average.data %>% left_join(df.states, by=c('State' = 'state.name')) 
  average.data$State <- NULL
  
  df$PSTATABB <- as.character(df$PSTATABB)
  df$PLFUELCT <- as.character(df$PLFUELCT)
  
  df <- df %>% left_join(average.data, by=c(c('PSTATABB' = 'state.abb'), 
                                            c('PLFUELCT' = 'fuel')))
  
  df$avg.fuel.cost <- ifelse(is.na(df$avg.fuel.cost), df$cost, df$avg.fuel.cost)
  df$PLHTRT <- ifelse(is.na(df$PLHTRT), df$heat.rate, df$PLHTRT)
  
  df$cost <- NULL
  df$heat.rate <- NULL
  
  # if there is still NA, use national average
  df <- df %>% left_join(national.average, by=c('PLFUELCT' = 'fuel'))
  df$avg.fuel.cost <- ifelse(is.na(df$avg.fuel.cost), df$cost, df$avg.fuel.cost)
  df$PLHTRT <- ifelse(is.na(df$PLHTRT), df$heat.rate, df$PLHTRT)
  df$cost <- NULL
  df$heat.rate <- NULL
  
  return(df)
}


get.utility.data <- function(utility.id=18642, bacode=NULL,
                             e.grid.data, needs.data, eia.data) {
  
  source('CENSUS_DATA.R')

  # Read E grid data ----
  # PLHTRT (heat rate) is in Btu/kwh
  if (!is.null(bacode)) {
    # if balancing authority code is informed, use this
    utility.data <- e.grid.data %>% filter(BACODE == bacode) %>% 
      select(PSTATABB, ORISPL, PNAME, CNTYNAME, LAT, LON, PLFUELCT, CAPFAC, 
             NAMEPCAP, NBFACTOR, PLNGENAN, PLHTRT)
    
  } else {
    # otherwise use utility.id (default to TVA)
    utility.data <- e.grid.data %>% filter(UTLSRVID == utility.id) %>% 
      select(PSTATABB, ORISPL, PNAME, CNTYNAME, LAT, LON, PLFUELCT, CAPFAC, 
             NAMEPCAP, NBFACTOR, PLNGENAN, PLHTRT)
    
  }
  
  utility.data <- utility.data %>% filter(!(PLFUELCT %in% c('OTHF', 'OFSL')))
  
  # added this line to take into consideration the energy in pumped hydro plant at TVA
  # originally the capacity factor was zero
  utility.data$CAPFAC <- round(abs(utility.data$PLNGENAN)/(utility.data$NAMEPCAP*8760), 3)
  
  summary.by.source <- utility.data %>% dplyr::group_by(PLFUELCT) %>%
    dplyr::summarize(capacity = sum(NAMEPCAP, na.rm=TRUE))
  summary.by.source <- as.data.frame(summary.by.source)
  print(summary.by.source)
  
  # Read NEEDS data ----
  utility.data.needs <- needs.data %>% 
    filter(ORIS.Plant.Code %in% utility.data$ORISPL) %>%
    select(Plant.Name, ORIS.Plant.Code, PlantType, State.Name, State.Code,
           FIPS5, Capacity..MW., Heat.Rate..Btu.kWh.) %>% 
    group_by(ORIS.Plant.Code) %>%
    dplyr::summarise(Plant.Name = Plant.Name[1],
                     PlantType = PlantType[1],
                     State.Name = State.Name[1],
                     State.Code = State.Code[1],
                     FIPS5 = FIPS5[1], 
                     Cap.MW = sum(Capacity..MW.), 
                     HR.Btu.kWh = mean(Heat.Rate..Btu.kWh.))
  
  # merge capacity data from NEEDS with eGRID data
  needs.merge <- utility.data.needs %>% select(ORIS.Plant.Code, Cap.MW)
  
  utility.data <- utility.data %>% left_join(needs.merge, 
                                             by = c("ORISPL" = "ORIS.Plant.Code"))
  
  utility.data <- utility.data %>% mutate(AVAILFAC = pmin(Cap.MW/NAMEPCAP, 1, 
                                                          na.rm = TRUE))
  
  # create the capacity data that will be used by the LP
  # WIND, BIOMASS: CAPACITY FACTOR
  # THERMAL: AVAILABLE POWER
  # HYDRO: AVAILABLE POWER (but there will be an energy constraint for these plants)
  utility.data$FINAL.CAP <- rep(0, nrow(utility.data))
  
  #idx1 <- utility.data$PLFUELCT %in% c('HYDRO', 'WIND', 'BIOMASS')
  
  idx1 <- utility.data$PLFUELCT %in% c('WIND', 'BIOMASS')
  utility.data$FINAL.CAP[idx1] <- (utility.data$CAPFAC * utility.data$NAMEPCAP)[idx1]
  utility.data$FINAL.CAP[!idx1] <- utility.data$Cap.MW[!idx1]
  
  utility.data <- utility.data[ ,c("PSTATABB", "ORISPL", "PNAME", "CNTYNAME", "LAT",
                                   "LON", "PLFUELCT", "CAPFAC", "NAMEPCAP", "Cap.MW",
                                   "AVAILFAC", "FINAL.CAP", "NBFACTOR", 
                                   "PLNGENAN", "PLHTRT")]
  
  utility.data <- utility.data %>% filter(!is.na(FINAL.CAP))
  
  # Read Eia 923 data ----
  
  utility.fuel <- eia.data %>% filter(Plant.Id %in% utility.data$ORISPL)

  utility.fuel.2 <- utility.fuel %>% group_by(Plant.Id) %>%
    dplyr::summarise(avg.fuel.cost = {(sum(FUEL_COST*Total.Heat,na.rm=TRUE)/
                                         sum(Total.Heat, na.rm=TRUE))/100})
  
  utility.data <- utility.data %>% left_join(utility.fuel.2, 
                                             by = c('ORISPL' = 'Plant.Id'))
  
  # fill up heat rate and fuel costs of thermal plants with NA or zero
  utility.data <- fill.up.costs.thermal(utility.data)

  # fill up data for nuclear (source: NEEDS and EPA IPM)
  hr <- get.national.heat.rates()
  utility.data$PLHTRT[utility.data$PLFUELCT == "NUCLEAR"] <- hr$heat.rate[hr$type == 'NUCLEAR']
  utility.data$avg.fuel.cost[utility.data$PLFUELCT == "NUCLEAR"] <- 0.751
  
  utility.data$gen.cost.fuel <- utility.data$PLHTRT*(1e-3)*utility.data$avg.fuel.cost
  utility.data$gen.cost.fuel[utility.data$PLFUELCT %in% c('HYDRO', 'WIND', 'SOLAR')] <- 0
  
  utility.data$total.cost <- utility.data$gen.cost.fuel

  # compute summary and print on screen
  cat('\n--------------------------------------------')
  cat('\nSummary of power plant portfolio:\n\n')
  print(utility.data %>% group_by(PLFUELCT) %>% 
          dplyr::summarise(capacity=sum(FINAL.CAP, na.rm = TRUE), 
                           cost = weighted.mean(x=total.cost, 
                                                w=FINAL.CAP, na.rm=TRUE)))
  print(sum(utility.data$FINAL.CAP, na.rm = TRUE))
  cat('\n\n--------------------------------------------\n')
  
  return(utility.data)
  
}

plot.supply.curve <- function(utility.data) {
  library(ggplot2)
  library(colorspace)

  # sort plants by generation cost
  if ('name' %in% names(utility.data)) {
    # if there is a column name (with the company name) then group rows by company    
    utility.data.2 <- utility.data[order(utility.data$name,
                                         utility.data$total.cost, 
                                         utility.data$PNAME), ]
    utility.data.2 <- utility.data.2 %>% group_by(name) %>% 
      mutate(sum.cap=cumsum(FINAL.CAP)) %>% as.data.frame()
  } else {
    utility.data.2 <- utility.data[order(utility.data$total.cost, 
                                         utility.data$PNAME), ]
    utility.data.2$sum.cap <- cumsum(utility.data.2$FINAL.CAP)
  }

  
  shape.pallete <- c(1, 0, 2, 3, 4, 6, 5, 7, 8)
  colour.pallete <- c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99',
                      '#e31a1c','#fdbf6f','#ff7f00','#cab2d6')
  
  g <- ggplot()
  g <- g + geom_point(data=utility.data.2,
                      aes(x=sum.cap/1e3, y=total.cost,
                          shape=PLFUELCT,
                          colour=PLFUELCT),
                      stroke = 1.1)
  g <- g + theme_bw() 
  g <- g + scale_color_manual(values = colour.pallete)
  g <- g + scale_shape_manual(values = shape.pallete)
  g <- g + theme(legend.position=c(0.05,0.95), legend.justification=c(0, 1))
  g <- g + guides(shape=guide_legend(title=NULL), 
                  colour=guide_legend(title=NULL))
  g <- g + xlab('Cumulative Capacity (GW)')
  g <- g + ylab('Short-run Marginal Cost ($/MWh)')
  
  return(g)
}

plot.supply.curve.facet <- function(utility.data) {
  
  g <- plot.supply.curve(utility.data = utility.data)
  g <- g + facet_wrap(~ name , ncol = 4)
  
  return(g)
  
}

read.list.utilities.serc <- function() {
  # the data in this csv file was organized in the google spreadsheet
  # https://goo.gl/r8eqVa
  df <- read.csv(file='/Users/kiko/GoogleDrive/CMU/RIPS/SERC companies/serc_companies_codes.csv')
  return(df)
}


read.needs.data <- function() {
  cat('\n\nReading NEEDS data to get estimated availability factor ...\n')
  
  needs.url <- 'https://www.epa.gov/sites/production/files/2015-08/needs_v515.xlsx'
  
  needs.file <- tempfile()
  download.file(url = needs.url, destfile = needs.file)
  
  f.csv <- tempfile(fileext = ".csv")
  f.xls <- needs.file
  
  convert.xls2csv(xls.file = f.xls, csv.file = f.csv)
  
  unlink(needs.file)
  unlink(f.xls)
  gc()
  
  needs.data <- read.csv(file=paste0(f.csv,'.0')) # plant data (skip first row)
  
  return(needs.data)
}

read.e.grid.data <- function(egrid.path='/Users/kiko/GoogleDrive/CMU/RIPS/egrid/all_egrid2014v2_files/',
                             egrid.file='eGRID2014_Data_v2.xlsx'){ 
  cat('\n\nReading E-GRID data ...\n')
  
  f.csv <- paste0(tempdir(),"/out.csv")
  f.xls <- paste0(egrid.path, egrid.file)
  
  convert.xls2csv(xls.file = f.xls, csv.file = f.csv)
  
  list.csv <- list.files(path=tempdir(), pattern = "out.csv", 
                         full.names = TRUE)
  
  plant.data <- read.csv(file=list.csv[4], skip = 1) # plant data (skip first row)
  unlink(list.csv)
  gc()
  
  return(plant.data)
}

read.eia.923.data <- function(eia1.path = '/Users/kiko/GoogleDrive/CMU/RIPS/egrid/',
                              eia1.file = 'f923_2014.zip'){
  # Fields in data frame
  #
  # Quantity: Quantity of fuel received in tons, barrel, or Mcf. Numeric
  # Average_Heat Content: Heat content of the fuel in millions of Btus per physical unit to the nearest 0.01 percent.
  # Fuel-Cost: All costs incurred in the purchase and delivery of the fuel to the plant in cents per million Btu(MMBtu) to the nearest 0.1 cent. Numeric.
  
  
  cat('\n\nReading EIA 923 data to get estimated fuel costs ...\n')
  
  f.zip <- paste0(eia1.path, eia1.file)
  
  # get name of excel file in zip file (larger file)
  a <- unzip(zipfile = f.zip, list = TRUE)
  name.xls <- a$Name[which(a$Length == max(a$Length))]
  
  #extract excel file and read to R
  f.xls <- unzip(zipfile = f.zip, files = name.xls, exdir = tempdir())
  
  f.csv <- paste0(tempdir(),"/out.csv")
  convert.xls2csv(xls.file = f.xls, csv.file = f.csv) 
  unlink(f.xls)
  gc()
  
  list.csv <- list.files(path=tempdir(), pattern = "out.csv", 
                         full.names = TRUE)
  eia.data <- read.csv(file=list.csv[5], skip = 4) # plant data (skip first row)
  eia.data <- eia.data %>% select(Plant.Id, Plant.Name, Plant.State, FUEL_GROUP, 
                                    QUANTITY, Average.Heat.Content, FUEL_COST,
                                    Operator.Name, Operator.Id)
  eia.data$Total.Heat <- eia.data$QUANTITY * eia.data$Average.Heat.Content
  eia.data$FUEL_COST <- as.numeric(as.character(eia.data$FUEL_COST))

  state.data <- eia.data %>% group_by(Plant.State, FUEL_GROUP) %>%
    mutate(total.cost=FUEL_COST*Total.Heat) %>%
    summarise(avg.state.cost=sum(total.cost, na.rm=TRUE)/sum(Total.Heat, na.rm=TRUE))
  
  national.data <- eia.data %>% group_by(FUEL_GROUP) %>%
    mutate(total.cost=FUEL_COST*Total.Heat) %>%
    summarise(avg.national.cost=sum(total.cost, na.rm=TRUE)/sum(Total.Heat, na.rm=TRUE))
  
  eia.data.2 <- eia.data %>% left_join(state.data, by= c('Plant.State', 'FUEL_GROUP'))
  
  eia.data.3 <- eia.data.2 %>% left_join(national.data, by= c('FUEL_GROUP'))
  
  eia.data.3$FUEL_COST <- ifelse(is.na(eia.data.3$FUEL_COST), 
                                 eia.data.3$avg.state.cost, 
                                 eia.data.3$FUEL_COST)
  
  eia.data.3$FUEL_COST <- ifelse(eia.data.3$FUEL_COST == 0, 
                                 eia.data.3$avg.national.cost, 
                                 eia.data.3$FUEL_COST)
  
  unlink(list.csv)
  rm(eia.data)
  rm(eia.data.2)
  
  eia.data.3$avg.state.cost <- NULL
  eia.data.3$avg.national.cost <- NULL
  gc()
  
  return(eia.data.3)
}

read.eia.web.report <- function(){
  
  cat('\n\nReading O&M costs from EIA website ...\n')
  eia.url.oem <- 'https://www.eia.gov/electricity/annual/html/epa_08_04.html'
  
  a <- getURL(eia.url.oem)
  b <- htmlParse(a)
  # get XML nodes marked as tbody (tables) and convert to R list
  ww <- getNodeSet(b, path = "//tbody")
  
  # get first table (operation and maintenance)
  # see names of columns in url
  www <- readHTMLTable(doc = ww[[1]], 
                       header = FALSE,
                       stringsAsFactors = FALSE)
  
  costs.oem <- as.numeric(www[www[ ,1] == 2014, ])[-1]
  
  df.oem <- data.frame(fuel=c('NUCLEAR', 'COAL',	'HYDRO', 'OTHER'),
                       op=costs.oem[1:4],
                       maint=costs.oem[5:8])
  
  return(df.oem)
  
  # utility.data <- utility.data %>% left_join(df.oem, 
  #                                            by = c('PLFUELCT' = 'fuel'))
  # 
  # utility.data$op[is.na(utility.data$op)] <- df.oem$op[df.oem$fuel == 'OTHER']
  # utility.data$maint[is.na(utility.data$maint)] <- df.oem$maint[df.oem$fuel == 'OTHER']
  # 
  # # total generation cost
  # utility.data$total.cost <- rowSums(cbind(utility.data$gen.cost.fuel, 
  #                                          utility.data$op,
  #                                          utility.data$maint), 
  #                                    na.rm=TRUE) 
}


download.data.utilities <- function(){
  
  needs.data <- read.needs.data()
  e.grid.data <- read.e.grid.data()
  eia.data <- read.eia.923.data()
  
  list.serc <- read.list.utilities.serc()
  
  x <- lapply(X=list.serc$eia.code, FUN=get.utility.data, 
              e.grid.data = e.grid.data, needs.data = needs.data,
              eia.data = eia.data)
  names(x) <- list.serc$name
  
  x <- ldply(x, data.frame, .id = 'name')
  
  return(x)
  
}


compute.fuel.cost.state <- function(eia.data) {
  # this function use data in form EIA 923 to compute average costs of fuel
  # for each source and each state
  
  fuel.data <- eia.data %>% group_by(Plant.State, FUEL_GROUP) %>%
    mutate(total.cost=FUEL_COST*Total.Heat) %>%
    summarise(avg.cost=sum(total.cost, na.rm=TRUE)/sum(Total.Heat, na.rm=TRUE))
  
  national.data <- eia.data %>% group_by(FUEL_GROUP) %>%
    mutate(total.cost=FUEL_COST*Total.Heat) %>%
    summarise(avg.cost=sum(total.cost, na.rm=TRUE)/sum(Total.Heat, na.rm=TRUE))
}

read.eia.fuel.costs.states <- function() {
  
  states.names <- getStateList()$Name
  states.names <- c(states.names, 'U.S. Total')
  
  # table 7.17 (average cost of coal by state)
  a <- getURL('https://www.eia.gov/electricity/annual/html/epa_07_17.html')
  b <- htmlParse(a)
  # get XML nodes marked as tbody (tables) and convert to R list
  ww <- getNodeSet(b, path = "//tbody")
  # get first table
  www <- readHTMLTable(doc = ww[[1]], 
                       header = FALSE,
                       stringsAsFactors = FALSE)
  www <- www[ ,c(1,2)]
  names(www) <- c('State', 'cost')
  coal <- www %>% filter(tolower(State) %in% tolower(states.names))
  coal$fuel <- 'COAL'
  
  # table 7.18 (average cost of petroleum by state)
  a <- getURL('https://www.eia.gov/electricity/annual/html/epa_07_18.html')
  b <- htmlParse(a)
  # get XML nodes marked as tbody (tables) and convert to R list
  ww <- getNodeSet(b, path = "//tbody")
  # get first table
  www <- readHTMLTable(doc = ww[[1]], 
                       header = FALSE,
                       stringsAsFactors = FALSE)
  www <- www[ ,c(1,2)]
  names(www) <- c('State', 'cost')
  oil <- www %>% filter(tolower(State) %in% tolower(states.names))
  oil$fuel <- 'OIL'

  # table 7.20 (average cost of ng by state)
  a <- getURL('https://www.eia.gov/electricity/annual/html/epa_07_20.html')
  b <- htmlParse(a)
  # get XML nodes marked as tbody (tables) and convert to R list
  ww <- getNodeSet(b, path = "//tbody")
  # get first table
  www <- readHTMLTable(doc = ww[[1]], 
                       header = FALSE,
                       stringsAsFactors = FALSE)
  www <- www[ ,c(1,2)]
  names(www) <- c('State', 'cost')
  ng <- www %>% filter(tolower(State) %in% tolower(states.names))
  ng$fuel <- 'GAS'

  df.out <- rbind(coal, oil, ng)
  df.out$cost <- as.numeric(df.out$cost)
  
  hr <- get.national.heat.rates()
  
  df.out <- df.out %>% left_join(hr, by=c('fuel' = 'type'))
  
  return(df.out)
  
}

get.national.heat.rates <- function() {
  
  # Table 8.1. Average Operating Heat Rate for Selected Energy Sources,
  # 2005 through 2015 (Btu per Kilowatthour)
  #
  
  theurl <- getURL("https://www.eia.gov/electricity/annual/html/epa_08_01.html")
  tables <- readHTMLTable(theurl)
  table.heat.rates <- tables[[2]]
  table.heat.rates <- as.data.frame(lapply(table.heat.rates, 
                                           FUN={function(x){as.numeric(as.character(x))}}))
  table.heat.rates <- table.heat.rates %>% gather(key=type, value=heat.rate,
                                                  Coal:Nuclear) %>% 
    group_by(type) %>% summarise(heat.rate=mean(heat.rate, na.rm=TRUE)) %>%
    as.data.frame()
  
  table.heat.rates$type[table.heat.rates$type == 'Coal'] <- 'COAL'
  table.heat.rates$type[table.heat.rates$type == 'Natural.Gas'] <- 'GAS'
  table.heat.rates$type[table.heat.rates$type == 'Nuclear'] <- 'NUCLEAR'
  table.heat.rates$type[table.heat.rates$type == 'Petroleum'] <- 'OIL'
  return(table.heat.rates)
  
}


