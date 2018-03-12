# July 2017
#
# this file contains all functions to execute the MILP that calculates
# the optimal dispatch of a set of power plants given a projected demand
# and analyze the results
#

run.LP <- function(hourly.demand,
                   utility.data, 
                   turndow.constraints=data.frame(type= c("OIL", "NUCLEAR", 
                                                          "COAL", "GAS", 
                                                          "HYDRO", "BIOMASS", 
                                                          "WIND", "AUX"),
                                                  min.value = c(0, 0.5, 0.5, 0,
                                                                0, 0, 0, 0)),
                   print.screen=TRUE, use.binary=TRUE) {
  
  require(lpSolve, quietly=TRUE, warn.conflicts=FALSE)
  require(plyr, quietly=TRUE, warn.conflicts=FALSE)
  require(dplyr, quietly=TRUE, warn.conflicts=FALSE)
  
  # create optimization problem ----
  #
  # Args:
  #   hourly.demand: data frame with 2 columns (time and load). The size of this
  #                  data frame will define the horizon of the simulation
  #   utility.data: data frame with data of power plants owned by an utility
  #                 (check format of columns)
  #   turndow.constraints: data frame with minimum generation value for each type.
  #                 (see default values above)
  #   print.screen: if TRUE it prints info of the problem on screen
  #   use.binary: if TRUE the problem includes the binary variables for
  #               the turndown constraints
  #
  # Returns:
  #   sol: list with solution of LP
  #
  # Decision variables:
  # x_{p,h} : dispatched power by plant p in hour h
  # P: number of power plants
  # H: number of hours in simulation horizon
  # p \in \{1, \ldots , P \}
  # h \in \{1, \ldots , H \}
  
  if (print.screen) cat('\n--------------------------------------------\n')
  ptm <- proc.time()
  if (print.screen) timestamp()
  if (print.screen) cat('Initializing LP...\n')
  
# general parameters ----
  M = 1e10 # big M
  P = nrow(utility.data)
  H = nrow(hourly.demand)
  
  if (print.screen) cat('P =', P,'\n')
  if (print.screen) cat('H =', H,'\n')

  # auxiliary vectors
  cap.plants <- utility.data$Cap.MW
  cost.plants <- utility.data$total.cost
  type.plants <- utility.data$PLFUELCT
  capfac.plants <- utility.data$CAPFAC
  
  load <- hourly.demand$load
  
  # Objective function ----
  
  # cost vector
  
  # __marginal cost of each plant ----
  cost.obj <- rep(cost.plants, H)
  
  # get indexes of hydro plants
  idx.hyd <- which(type.plants == "HYDRO")
  num_hyd <- length(idx.hyd)
  
  # __add binary variables for plants with minimum dispatch constraints ----
  # get fuel types that have turndown constraints
  tdc.fueltypes <- turndow.constraints$type[turndow.constraints$min.value > 0]
  flag.tdc <- type.plants %in% tdc.fueltypes
  idx.tdc <- which(flag.tdc)
  n.plants.tdc <- sum(flag.tdc)
  
  # __write cost vector ----
  
  if (use.binary) cost.obj <- c(cost.obj, rep(M, n.plants.tdc))
  
  # __number of decision variables (including auxiliary binary) ----
  n.real <- P*H
  n.bin <- ifelse(use.binary, n.plants.tdc, 0)
  n.var <- n.real + n.bin
  
  idx.real <- 1:n.real
  if(use.binary) {
    idx.bin <- c((n.real+1):n.var)
  } else {
    idx.bin <- NULL
  }
  
  if (print.screen) cat('# of decision variables =', n.real,'\n')
  if (print.screen) cat('# of auxiliary binary variables =', n.bin,'\n')
  
  # __number of constraints
  if(use.binary){
    n.con <- H + P*H + H*length(idx.tdc) + H*length(idx.tdc) + num_hyd
  } else{
    n.con <- H + P*H + H*length(idx.tdc) + num_hyd
  }
  
  if (print.screen) cat('# of contraints =', n.con,'\n')
  
  # Constraints ----
  if (print.screen) cat('Creating Constraint Matrix...\n')
  
  # create empty matrix for constraints
  A.con <- NULL
  b.rhs <- NULL
  const.dir <- NULL
  
  # __for each hour of the day, \sum_p x_{p,h} >= dem_h ----
  # (H contraints)
  if (print.screen) cat('generation > demand \n')
  
  # auxiliary matrix
  A.con.aux <- matrix(0, nrow = H,  n.var)
  
  for (h in 1:H) {
    # row with constraint for hour h
    idx1 <- (h-1)*P + 1
    idx2 <- h*P
    A.con.aux[h, idx1:idx2] <- rep(1, P)
  }
  A.con <- rbind(A.con, A.con.aux)
  b.rhs <- c(b.rhs, load)
  const.dir <- c(const.dir, rep(">=", H))
  
  # __for each hour of the day, x_{p,h} <= CAP_p ----
  # (H*P contraints)
  if (print.screen) cat('generation < capacity \n')
  
  if(use.binary){
    A.con <- rbind(A.con, cbind(diag(nrow = P*H), 
                                matrix(0, nrow=P*H, ncol=n.plants.tdc)))
  } else{
    A.con <- rbind(A.con, diag(nrow = P*H))
  }
  
  b.rhs <- c(b.rhs, rep(cap.plants, H))
  const.dir <- c(const.dir, rep("<=", P*H))
  
  # __minimum dispatch constraints ----
  # if plant is on it must generate >= minimum value 
  # (H*len(idx.tdc) contraints)
  if (print.screen) cat('minimum dispatch constraints \n')
  
  for (i in idx.tdc) {
    # mimimum value of generation of this type of power plant
    min.value <- (turndow.constraints %>% 
                    filter(as.character(type) == as.character(type.plants[i])) %>%
                    select(min.value))$min.value
    
    # for this plant allocate auxiliary matrix with zeros (size H by n.var)
    m.aux <- matrix(0, nrow = H, ncol = n.var)
    # define indexes of plant variables in aux matrix
    idx.aux <- cbind(c(1:H), seq(from=i, by=P, length.out = H))
    # set selected elements of matrix to 1
    m.aux[idx.aux] <- 1
    
    if(use.binary) {
      i.bin <- H*P + which(idx.tdc == i)
      m.aux[ , i.bin] <- -1*min.value*cap.plants[i]
      rhs.value <- 0
    } else {
      rhs.value <- min.value*cap.plants[i]
    }
    
    # add to constraint matrix
    A.con <- rbind(A.con, m.aux)
    b.rhs <- c(b.rhs, rep(rhs.value, H))
    const.dir <- c(const.dir, rep(">=", H))
  }
  
  # __either on or off for all h constraints ----
  # for each plant p with minimum dispatch 
  # x_{p, h} <= M*y_p  \forall h <=> x_{p, h} - M*y_p <= 0  \forall h
  # (H*len(idx.tdc) contraints)
  if (use.binary) {
    for (i in idx.tdc) {
      # allocate auxiliary matrix with zeros
      m.aux <- matrix(0, nrow = H, ncol = n.var)
      # define indexes of plant variables in aux matrix
      idx.aux <- cbind(c(1:H), seq(from=i, by=P, length.out = H))
      # set selected elements of matrix to 1
      m.aux[idx.aux] <- 1
      
      i.bin <- H*P + which(idx.tdc == i)
      
      m.aux[ , i.bin] <- -1*M
      
      # add to constraint matrix
      A.con <- rbind(A.con, m.aux)
      b.rhs <- c(b.rhs, rep(0, H))
      const.dir <- c(const.dir, rep("<=", H))
    }
  }

  # __energy constraint for Hydro plants ----
  # for each plant p with type HYDRO
  # \sum_h x_{p, h} <= CAPFAC_p * CAP_p * H   \forall p \in HYDRO
  # (num_hyd contraints)
  if (print.screen) cat('energy constraint for hydro \n')
  
  m.aux <- matrix(0, nrow = num_hyd, ncol = n.var)
  
  # define indexes in m.aux that will be set to 1
  a <- as.vector(sapply(X=idx.hyd, FUN=seq, by=P, length.out=H))
  idx.aux <- cbind(base::rep(seq(from=1, to=num_hyd), rep(H, num_hyd)), a) 
  
  m.aux[idx.aux] <- 1
    
  # add to constraint matrix
  A.con <- rbind(A.con, m.aux)
  b.rhs <- c(b.rhs, (capfac.plants * cap.plants * H)[idx.hyd])
  const.dir <- c(const.dir, rep("<=", num_hyd))
  
  if (print.screen) cat('Constraint Matrix complete!\n\n')
  
  if (print.screen) cat('Running solver... \n\n')
  # Run optimization problem ----
  sol <- lp (direction = "min", 
             objective.in=cost.obj, 
             const.mat=A.con, 
             const.dir=const.dir, 
             const.rhs=b.rhs,
             binary.vec=idx.bin) 
  if (print.screen) cat('Done!\n')
  elapsed.time <- proc.time() - ptm
  if (print.screen) cat('Elapsed time: ', elapsed.time[3], 'sec\n', sep = '')
  if (print.screen) cat('--------------------------------------------\n')
  return(sol)
}

parse.solution <- function(x.sol, plants.data, hourly.demand, 
                           normalize=FALSE) {
  # returns data frame with solution of the LP
  # in a shape easier to visualize
  #
  # xsol: vector of length P*H with solution of decision variables
  # plants.data: data frame with data of power plants
  # hourly.demand: vector with values of load for each hour
  # normalize: if TRUE dispatch results are fraction of capacity
  #
  # returns:
  # matrix.out: 
  #     representation of xsol as a P x H matrix
  
  P <- nrow(plants.data)
  H <- length(hourly.demand)
  matrix.out <- matrix(x.sol[1:(P*H)], nrow = P)
  matrix.out <- round(matrix.out, 2)
  if (normalize) {
    cap.matrix <- matrix(rep(plants.data$FINAL.CAP, H), nrow = P)
    matrix.out <- matrix.out / cap.matrix
  }
  
  df.out <- data.frame(plant=c('Load', as.character(plants.data$PNAME)), 
                       type=c(NA, plants.data$PLFUELCT),
                       cost=c(NA, plants.data$total.cost),
                       rbind(matrix(hourly.demand, nrow=1), matrix.out),
                       stringsAsFactors=FALSE)
  
  names(df.out)[4:ncol(df.out)] <- formatC(1:H, width=2, 
                                           format='d', flag = '0')
  
  return(df.out)
}

plot.dispatch.by.type <- function(df.solution) {
  require(plyr)
  require(dplyr)
  require(tidyr)
  require(ggplot2)
  
  df <- as.data.frame(df.solution %>% group_by(type) %>% 
                        summarise_at(vars(num_range(prefix='', range=1:24, width = 2)), sum) %>%
                        gather(key=hour, value=value, `01`:`24`))
  df$hour <- as.numeric(df$hour)
  
  df$type <- factor(df$type, levels = rev(c("HYDRO", "NUCLEAR", 
                                            "COAL", "GAS", "OIL")))
  
  g <- ggplot(data = df) + 
    geom_area(aes(x=hour, y=value, fill=type)) +
    theme_bw()
  g <- g + xlab("Hour of the day") + ylab("Total Generation (MWh)")
  return(g)
}

plot.compare.total.gen.by.type <- function(df1, df2, case.names) {
  require(plyr)
  require(dplyr)
  require(tidyr)
  require(ggplot2)
  
  df1$case <- case.names[1]
  df2$case <- case.names[2]
  
  df.total <- rbind(df1, df2)
  
  df <- df.total %>% group_by(type, case) %>% 
    summarise_at(vars(num_range(prefix='', range=1:24, width = 2)), sum) %>%
    gather(key=hour, value=value, `01`:`24`) %>% group_by(case, type) %>%
    summarise(tot.gen = sum(value)/1e3)
  df <- as.data.frame(df)
  
  df$type <- factor(df$type, levels = rev(c("HYDRO", "NUCLEAR", 
                                            "COAL", "GAS", "OIL")))
  df$case <- as.factor(df$case)
  g <- ggplot(data = df) + 
    geom_bar(aes(x=case, y=tot.gen, fill=type), stat='identity', width = 0.6) +
    theme_bw()
  g <- g + xlab("Case") + ylab("Total Generation (GWh)")
  g
  return(g)
}


compute.solution <- function(x.sol, plants.data, H=24) {
  # the LP contains artificial binary variables to emulate the 
  # minimum dispatch constraints. This function recomputes the objective value
  # without these binary variables.
  
  P <- nrow(plants.data)
  # get solutions of "real" generation variables
  x.real <- matrix(x.sol[1:(P*H)], nrow = P) 
  
  # cost vector of the generation variables
  cost.obj <- rep(plants.data$total.cost, H)
  
  # compute final objective value cost
  obj.value.final <- sum(x.real * cost.obj)
  
  return(obj.value.final)
}


compute.cost.by.type <- function(df.sol) {
  # computes data frame with generation cost for each type
  #
  # Args:
  # df.sol: data frame with solution of LP. The format of df.sol must be
  #         the same as the output of the function 'parse.solution'
  #
  # Returns:
  # df.cost:
  #
  
  require(plyr)
  require(dplyr)
  require(tidyr)

  df.cost <- df.sol %>% 
    mutate_at(vars(num_range(prefix='', 1:H, width = 2)),
              {function(x){x*df.sol$cost}}) %>%
    select(-cost) %>% gather(key=hour, value=gen.cost, `01`:`24`) %>%
    group_by(type) %>% summarise(value=sum(gen.cost)/1e6)
  df.cost <- as.data.frame(df.cost)
  
  return(df.cost)
}

wrapper.runLP <- function(hourly.demand,
                          utility.data, 
                          turndow.constraints=data.frame(type= c("OIL", "NUCLEAR", 
                                                                 "COAL", "GAS", 
                                                                 "HYDRO", "BIOMASS", 
                                                                 "WIND"),
                                                         min.value = c(0, 0.5, 
                                                                       0.5, 0,
                                                                       0, 0, 0)),
                          print.screen=FALSE, use.binary=TRUE) {
  # wrapper function to run LP and format solution
  
  ptm <- proc.time()
  day.value <- unique(lubridate::date(hourly.demand$time))
  season <- unique(as.character(hourly.demand$season))
  
  H = nrow(hourly.demand)
  result.lp <- run.LP(hourly.demand,
                      utility.data, 
                      turndow.constraints,
                      print.screen = print.screen,
                      use.binary = use.binary)
  obj.val <- compute.solution(result.lp$solution, 
                              plants.data = utility.data,
                              H=H)
  df.solution <- parse.solution(result.lp$solution, 
                                plants.data = utility.data,
                                hourly.demand=hourly.demand$load)
  list.out <- list(day=day.value,
                   season=season,
                   status=result.lp$status,
                   obj.value=obj.val,
                   df.sol=df.solution)
  
  elapsed.time <- proc.time() - ptm
  cat('---- Day ', format(day.value, '%Y-%m-%d'), ' complete! ', 
      'Elapsed time: ', elapsed.time[3], 'sec-----\n', sep = '')
  return(list.out)
}


consolidate.daily.results <- function(list.results){
  # reads list with results of daily dispatches and consolidates into
  # a single data frame in long format
  
  require(plyr)
  require(dplyr)
  require(tidyr)
  
  # for each element of the list, construct a data frame
  # with columns for the day and the season
  a <- lapply(X=list.results, FUN={function(x){
    return(data.frame(status=x$status, day=x$day, 
                      season=x$season, x$df.sol))}})
  # merge all elements of the list in a single data frame
  b <- ldply(a, data.frame)
  
  # convert to long format
  c <- b %>% gather(key="hour", value="gen", X01:X24)
  
  # convert hour column to numeric
  c$hour <- as.numeric(gsub("X", "", c$hour))
  
  return(c)
}

run.multiple.LP <- function(list.load, utility.data, run.parallell=TRUE, 
                            n.clusters=NA, use.binary=FALSE) {
  # run multiple linear programming problems in parallel or in a single core
  #
  # Args:
  # list.load: list with demand curves for different periods being optimized.
  #             Usually each element of the list is a single 24 h day.
  # run.parallel: (boolean) If TRUE LPs are executed in parallel using ncores-1
  #
  
  n.days <- length(list.load)
  
  if(run.parallell) {
    if(is.na(n.clusters)) n.clusters <- parallel::detectCores()-1
    # overwrite old file and write number of clusters
    cat(n.clusters,'\n', file='~/out_parallel.txt', sep="", append=FALSE)
    parallelCluster <- parallel::makeCluster(n.clusters,
                                             outfile="~/out_parallel.txt")
    #print(parallelCluster)
    
    tryCatch({
      #start parallel processing. write start time
      cat(floor(proc.time()[3]), '\n', file = "~/out_parallel.txt", sep = "", 
          append = TRUE)
      cat(n.days,'\n', file = "~/out_parallel.txt", sep = "", append = TRUE)
      parallel::clusterExport(cl=parallelCluster, 
                              varlist=c("run.LP", "compute.solution", 
                                        "parse.solution"))
      list.out <- parallel::parLapply(cl=parallelCluster, 
                                      X=list.load[1:n.days], fun=wrapper.runLP,
                                      utility.data=utility.data,
                                      use.binary=use.binary)}, 
      error = function(e) print(e))
    
    # end of parallel processing. write end time and 'DONE!'
    cat(floor(proc.time()[3]), '\n', file = "~/out_parallel.txt", sep = "", 
        append = TRUE)
    cat('DONE!\n', file = "~/out_parallel.txt", sep = "", append = TRUE)
    # Shutdown cluster neatly
    if(!is.null(parallelCluster)) {
      parallel::stopCluster(parallelCluster)
      parallelCluster <- c()
    }
    
  } else {
    list.out <- lapply(X=list.load[1:n.days], FUN=wrapper.runLP,
                       utility.data=utility.data)
  }
  
  df.results <- consolidate.daily.results(list.out)
  
  return(df.results)
}

check.constraint.plant <- function() {
  
  x.sol <- sol$solution
  
  name.plant <- "Colbert"
  idx.plant <- which(utility.data$PNAME == name.plant)
  
  idx.var.plant <- seq(from=idx.plant, by=P, length.out = H)
  x.sol[idx.var.plant]
  
  #  H + P*H + H*length(idx.tdc) + H*length(idx.tdc) + num_hyd
  
  # capacity constraint
  idx.zero <-  H
  idx.const.plant <- idx.zero + seq(from=idx.plant, by=(P), length.out = H)
  
  # minimum dispatch
  idx.zero <- H + P*H
  idx2 <- which(idx.tdc == idx.plant)
  idx.ini <- (idx2-1)*H + 1
  seq.1 <- seq(from=idx.ini, by = 1, length.out=H)
  idx.const.plant <- c(idx.const.plant, (idx.zero + seq.1))

  # on and off
  idx.zero <- H + P*H + H*length(idx.tdc)
  idx2 <- which(idx.tdc == idx.plant)
  idx.ini <- (idx2-1)*H + 1
  seq.1 <- seq(from=idx.ini, by=1, length.out=H)
  idx.const.plant <- c(idx.const.plant, (idx.zero + seq.1))
    
  rr <- A.con %*% x.sol - b.rhs
  rr <- rr[idx.const.plant]
  
  
  
}
