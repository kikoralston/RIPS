# This script runs the complete pipeline of analysis

# load other source files

# get data frame with data from TVA power plants
source("supplycurve.R")

# function to create demand projection
source("rips_demand.R")

# power plant dispatch LP
source("dispatch_LP.R")

divide.load.daily.curves <- function(range.dates, load.data, n.clusters=NA) {
  # divides data frame with hourly load data into list where each element
  # is one 24 hour day of data
  # performs this using parallel computing
  
  if(is.na(n.clusters)) n.clusters <- parallel::detectCores()-1
  
  parallelCluster <- parallel::makeCluster(n.clusters)
  
  tryCatch({
    list.out <- parallel::parLapply(cl=parallelCluster, 
                                    X=range.dates, 
                                    fun={function(d, y){
                                      dplyr::filter(y, 
                                                    lubridate::date(time) == d)}},
                                    y=load.data)}, 
    error = function(e) print(e))
  
  # Shutdown cluster neatly
  if(!is.null(parallelCluster)) {
    parallel::stopCluster(parallelCluster)
    parallelCluster <- c()
  }
  
  return(list.out)
}

dispatch.future <- function(lm.model, gcm.data.path, base.line.data,
                            utility.data, n.clusters=NA,
                            type=c('original', 'paulina'),
                            ref.period=c(2005, 2015), 
                            future.period=c(2035, 2045)) {
  # This function runs a complete pipeline of analysis for a given 
  # GCM future scenario and given weather and load historical data
  
  ## evaluate choices
  type <- match.arg(type)
  
  # get future weather
  future.weather <- create.projection.weather(url.uw=gcm.data.path,
                                              base.line.data = base.line.data,
                                              type = type,
                                              ref.period=ref.period,
                                              future.period=future.period)
  
  # Future demand
  cat('\n-------- Computing future demand --------\n')
  load.future <- project.demand.2(lm.hourly=lm.model, 
                                  weather.data = future.weather)
  
  range.dates <- unique(date(load.future$time))
  list.load <- divide.load.daily.curves(range.dates = range.dates, 
                                        load.data = load.future,
                                        n.clusters = n.clusters)
  # compute dispatch for future
  cat('\n-------- Running dispatch for future --------\n')
  df.results.future <- run.multiple.LP(list.load = list.load, 
                                       utility.data = utility.data,
                                       n.clusters = n.clusters)
  return(df.results.future)
}

dispatch.base.line <- function(lm.model, utility.data, base.line.data,
                               n.clusters=NA) {
  
  # estimate typical load for this baseline weather data
  cat('\n-------- Computing baseline demand --------\n')
  base.line.data <- compute.calendar.variables(base.line.data)
  base.line.load <- project.demand.2(lm.hourly = lm.model,
                                     weather.data = base.line.data)
  
  # divide total load in daily curves
  range.dates <- unique(date(base.line.load$time))
  list.load <- divide.load.daily.curves(range.dates = range.dates, 
                                        load.data = base.line.load,
                                        n.clusters = n.clusters)
  
  
  # run dispatch for the base line weather
  cat('\n-------- Running dispatch for baseline --------\n')
  df.results <- run.multiple.LP(list.load = list.load, 
                                utility.data = utility.data, 
                                n.clusters = n.clusters)
  
  return(df.results)
  
}

plot.mean.dispatch <- function(df.results) {
  # line plot with mean dispatch of each type in each season
  
  # remove rows with demand data
  df.result <- df.result %>% filter(!(plant == 'Load'))
  
  xx <- df.results %>% group_by(year, type, season, hour) %>% 
    dplyr::summarise(value=sum(gen, na.rm=TRUE))
  
  shape_pallete <- c(21, 22, 23, 24, 25)
  g <- ggplot()
  g <- g + geom_line(data=xx, aes(x=hour, y=value, colour=type)) + 
    geom_point(data=xx, aes(x=hour, y=value, shape=type, colour=type), 
               fill='white')
  g <- g + scale_shape_manual(values = shape_pallete)
  g <- g + facet_grid(year ~ season)
  g <- g + theme_bw()
  return(g)
}

plot.mean.cost <- function(df.results) {
  # bar plot with dispatch cost in each season mean dispatch of each type in each season
  
  # remove rows with demand data
  df.result <- df.result %>% filter(!(plant == 'Load'))
  
  yy <- df.results %>% dplyr::mutate(total.cost=gen*cost) %>% 
    group_by(year, type, season) %>% 
    dplyr::summarise(value=sum(total.cost, na.rm=TRUE))
  yy <- as.data.frame(yy)
  
  yyy <- yy %>% group_by(year, season) %>% summarise(value=sum(value)) 
  
  yyyy <- ((yyy %>% filter(year == 2040))$value/(yyy %>% filter(year == 2015))$value) -1
  names(yyyy) <- levels(yyy$season)
  print(yyyy)
  
  g <- ggplot()
  g <- g + geom_bar(data=yy, aes(x=as.factor(year), y=value/1e6, fill=type), stat='identity', 
                    width = 0.3)
  g <- g + theme_bw()
  g <- g + facet_wrap( ~ season, ncol = 2)
  g <- g + ylab('Cost (MM$)') + xlab('Year')
  return(g)
}


plot.capacity.factors <- function(df.result, utility.data) {
  
  # remove rows with demand data
  df.result <- df.result %>% filter(!(plant == 'Load'))
  
  df.1 <- left_join(df.results, select(utility.data, PNAME, Cap.MW), 
                    by=c("plant" = "PNAME"))
  
  df.1 <- df.1 %>% group_by(year, season, type) %>% 
    summarize(tot.gen=sum(gen, na.rm=TRUE), tot.cap=sum(Cap.MW, na.rm=TRUE)) %>%
    mutate(cap.factor = round(tot.gen/tot.cap,3)) %>% as.data.frame()
  
  g <- ggplot(data=df.1)
  g <- g + geom_bar(aes(x=as.factor(year), y=cap.factor*100, fill=type), 
                    position = "dodge", stat='identity', width = 0.65)
  g <- g + facet_wrap( ~ season, ncol = 2) + theme_bw()
  g <- g + xlab("year") + ylab("Capacity Factor (%)")
  return(g)
}


run.single.complete.case <- function(n.clusters = NA) {
  
  path.data <- './data/'
  gcm.data.path <- "../UW/meteo_memphis.txt"
  
  e.grid.data <- readRDS(file=paste0(path.data, 'egrid_data.rds'))
  needs.data <- readRDS(file=paste0(path.data, 'needs_data.rds'))
  eia.data <- readRDS(file=paste0(path.data, 'eia_data.rds'))
  
  utility.data <- get.utility.data(e.grid.data = e.grid.data, 
                                   needs.data = needs.data,
                                   eia.data = eia.data)
  
  df.historical <- read.data.frame(name.of.file=paste0(path.data, 
                                                       'RipsDemand.rds'))
  # fit regression model
  cat('\n-------- Fitting Regression model --------\n')
  lm.model <- fit.lm.model(df=df.historical)
  
  # get baseline weather data
  cat('\n-------- Reading Base line data --------\n')
  base.line.data <- create.base.line(name.file=paste0(path.data, 
                                                      'data_memphis_1973.rds'))
  
  # compute base line (present day) case
  df.results.present <- dispatch.base.line(lm.model, utility.data, 
                                           base.line.data, 
                                           n.clusters = n.clusters)
  
  # compute future case
  df.results.future <- dispatch.future(lm.model, gcm.data.path, base.line.data,
                                       utility.data,
                                       n.clusters = n.clusters)
  
  # combine results
  df.results.present <- data.frame(year=2015, df.results.present)
  df.results.future <- data.frame(year=2040, df.results.future)
  
  # combine results in a single data frame
  df.results <- rbind(df.results.present, df.results.future)
  
  return(df.results)
}

run.batch.cases <- function(json.file, n.clusters=NA) {
  # this function will run a complete batch of cases for a single 
  # region/balance authority. Information of this batch must be specified 
  # in the json.file. The json.file has field name 'cases' which contains one
  # array where each element is a list with 6 named fields:
  #
  # - name: string with name of case (typically of balancing area being simulated)
  #         which will be used to name the rds a file with the results
  # - data.frame: string with path to a rds file with  
  #               the data frame to be used to fit the regression model
  # - utility.data: string with complete path to rds file with the data 
  #                 frame containing utility power plant data 
  # - base.line: string with paths to a rds file with the  
  #              data frame containing historical weather data that will be  
  #              used to compute the baseline case and future cases
  # - ref.period: string of the form 'yyyy;yyyy' with initial and final year 
  #               of the period to be used as baseline case
  # - future.period: string of the form 'yyyy;yyyy' with initial and final year 
  #               of the period to be considered in the future projections
  # - type    : when creating the future projections use original or paulina
  #             method
  # - gcm.data: array of strings with paths to txt files with the GCM outputs
  #             for this specific region/balance authority
  #
  # Args:
  #   json.file: string with complete path to json file with information of batch cases
  #   n.clusters: integer with number of clusters to use
  
  require(jsonlite)
  
  # list.cases$cases will be a list of lists
  # each element of the list will be a list with 6 fields (specified above)
  list.cases <- read_json(path=json.file, simplifyVector = FALSE)
  
  n.cases <- length(list.cases$cases)
  
  for (j in 1:n.cases) {
    case <- list.cases$cases[[j]]
    
    cat('\n-------- ', case$name,' --------\n', sep='')
    
    ref.period <- as.numeric(unlist(strsplit(case$ref.period, split = ';')))
    future.period <- as.numeric(unlist(strsplit(case$future.period, split = ';')))
    
    utility.data <- readRDS(file = case$utility.data)
    
    df.historical <- readRDS(file=case$data.frame)
    
    # fit regression model
    cat('\n-------- Fitting Regression model --------\n')
    lm.model <- fit.lm.model(df=df.historical)
    
    # get baseline weather data
    cat('\n-------- Reading Base line data --------\n')
    base.line.data <- create.base.line(name.file=case$base.line,
                                       period=ref.period)
    
    list.results <- list()
    
    # compute base line (present day) case
    list.results[['present']] <- dispatch.base.line(lm.model, utility.data, 
                                                    base.line.data, 
                                                    n.clusters = n.clusters)
    
    # compute batch of future gcms
    n.gcm <- length(case$gcm.data)
    for (i in 1:n.gcm) {
      cat('\n---------------------------------------------\n')
      cat('                   GCM ', i, sep = '')
      cat('\n---------------------------------------------\n')
      df.results.future <- dispatch.future(lm.model, case$gcm.data[[i]], 
                                           base.line.data, utility.data,
                                           n.clusters = n.clusters,
                                           type = case$type,
                                           ref.period=ref.period, 
                                           future.period=future.period)
      name.gcm <- paste0('gcm.', i)
      list.results[[name.gcm]] <- df.results.future
    }
    
    saveRDS(list.results, file = paste0(case$name, '.rds'))
  }
  
  return(TRUE)
}


simulate.weather <- function(json.file, n.clusters=NA) {
  # creates simulation of weather data
  #
  # Args:
  #   json.file: string with complete path to json file with information of batch cases
  #   n.clusters: integer with number of clusters to use
  
  require(jsonlite)
  
  # list.cases$cases will be a list of lists
  # each element of the list will be a list with 6 fields (specified above)
  list.cases <- read_json(path=json.file, simplifyVector = FALSE)
  
  n.cases <- length(list.cases$cases)
  
  for (j in 1:n.cases) {
    case <- list.cases$cases[[j]]
    
    cat('\n-------- ', case$name,' --------\n', sep='')
    
    list.results <- list()

    # get baseline weather data
    cat('\n-------- Reading Base line data --------\n')
    
    # compute base line (present day) case
    base.line.data <- create.base.line(name.file=case$base.line)
    list.results[['present']] <- base.line.data
    
    # compute batch of future gcms
    n.gcm <- length(case$gcm.data)
    for (i in 1:n.gcm) {
      cat('\n---------------------------------------------\n')
      cat('                   GCM ', i, sep = '')
      cat('\n---------------------------------------------\n')
      
      name.gcm <- paste0('gcm.', i)
      list.results[[name.gcm]] <- create.projection.weather(url.uw=case$gcm.data[[i]],
                                                            base.line.data = base.line.data,
                                                            type = case$type)
    }
    
    saveRDS(list.results, file = paste0(case$name, '.rds'))
  }
  
  return(TRUE)
}

