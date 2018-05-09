# ************************************************************************
#           functions to create specific plots
# ************************************************************************

source("auxiliary_functions.R")

build.ldc <- function(load.data, time=c('year', 'month', 'day')) {
  # Computes load duration curve (LDC)
  #
  # Args:
  #   time: (string) time interval used to compute the LDC
  #   
  #
  # Returns:
  #   
  
  require(plyr)
  require(dplyr)
  require(ggplot2)
  require(lubridate)
  
  ## evaluate choices
  time <- match.arg(time)
  
  if (time == 'year') {
    
    load.data$gcm <- as.character(load.data$gcm)
    load.data$gcm[load.data$gcm != 'present'] <- 'gcm'
    
    out <- load.data %>% mutate(hour.year = (yday(day)-1)*24 + hour) %>%
      group_by(gcm, hour.year) %>% summarize(ldc=mean(gen, na.rm=TRUE)) %>%
      arrange(-ldc) %>% mutate(hour.year=seq(1, n())) %>% as.data.frame()
  }
  
  return(out)
}

compute.continuous.hours <- function(nf, value) {
  # Computes distribution of continous intervals (in hours) of demand above 
  # a certain value
  #
  # Args:
  #   nf: (string) name of file with simulation results
  #   value: (numeric) reference value of hourly demand (in MW)
  #   
  #
  # Returns:
  #   
  
  require(plyr)
  require(dplyr)
  require(ggplot2)
  require(lubridate)
  
  # reading data
  cat('Reading RDS file ', nf, '. Please wait...\n')
  rcp <- readRDS(file=nf)
  
  # filter demand
  cat('Get load values...\n')
  load <- lapply(X=rcp, 
                 FUN={function(x){
                   x %>% filter(plant == 'Load') %>% arrange(day, hour) %>%
                     select(season, day, hour, gen)}})
  cat('End of', nf, '! Freeing memory...\n')
  rm(rcp)
  gc()
  
  cat('Computing duration of events...\n')
  out <- lapply(X=load, 
                FUN={function(x, y){
                  cat('Computing durations for one case...')
                  
                  # create column with 1 if demand in that hour is > than value
                  z <- x %>% mutate(above.value = (gen >= y)*1)
                  
                  # sweep data frame to compute duration of events with demand > value
                  z <- z[!is.na(z$above.value), ]
                  dur <- rep(NA, nrow(z))
                  k <- 1
                  i <- 1
                  while (i <= nrow(z)){
                    if (z$above.value[i] > 0) {
                      j <- 0
                      while ((i+j <= nrow(z)) & (z$above.value[i+j] > 0)){
                        j <- j+1
                      }
                      dur[k] <- j
                      k <- k+1
                      i <- i + j
                    } else {
                      i <- i + 1
                    }
                  }
                  dur <- dur[!is.na(dur)]
                  
                  cat('Done!\n')
                  return(dur)
                }}, y=value)
  
  out <- ldply(.data=out, .fun=data.frame, .id = 'case')
  names(out)[2] <- 'duration'
  out$case<- as.character(out$case)
  out$case[out$case != 'present'] <- '2089-2099'
  
  ggplot(out, aes(x=duration)) + geom_density(aes(group=case, colour=case, fill=case), alpha=0.3)
  
  out %>% group_by(case) %>% summarise(x95=quantile(duration, probs=c(0.99)), mean=mean(duration)) 
}

plot.ldc <- function(load.data) {
  require(ggplot2)
  
  ldc.list <- lapply(load.data, build.ldc)
  
  ldc <- ldply(.data=ldc.list, .fun=data.frame, .id='rcp')
  ldc$rcp <- as.character(ldc$rcp)
  
  ldc$gcm[ldc$gcm != 'present'] <- ldc$rcp[ldc$gcm != 'present']
  ldc$rcp <- NULL
  
  g <- ggplot(data=ldc) + geom_line(aes(x=hour.year, y=ldc, linetype=gcm))
  g <- g + xlab('hour of year') + ylab('Hourly Load (MW)') + theme_bw()
  g <- g + theme(legend.position=c(0.95,0.95), legend.justification=c(1,1)) + 
    guides(linetype=guide_legend(title=NULL))
  print(g)
  
  return(g)
}

plot.temp.load <- function(reg.model,
                           data.set,
                           temp.breaks = c(-10, 0, 10, 20, 30),
                           temp = c(-15:40),
                           inter.term = NULL,
                           dir.out = "./out/",
                           add.marg.hist=FALSE,
                           normalize.values=FALSE){
  # Creates plot of temperature in X axis and estimated load in Y axis
  # (complete function)
  #
  # Args:
  #   reg.model: linear regression model object
  #   data.set: original data set used in regression model (before computing
  #             temperature components)
  #   temp.breaks:  breakpoints of temperature intervals in stepwise linear 
  #                 model
  #   temp: temperature values in X axis
  #   inter.term: (string) name of variable in interaction ("dew.point", "rh") 
  #               or NULL if no interaction
  #   add.marg.hist: (boolean) if TRUE we add marginal histograms to both axis
  #   normalize.values: (boolean) if TRUE we "normalize" the load values by 
  #                     subtracting the fixed effects from them
  #
  # Returns:
  #   nothing
  
  require(gridExtra)
  
  # compute temperature components and create X matrix  
  t.c <- createTempComponents(temp, temp.breaks = temp.breaks)
  
  X.2 <- as.matrix(cbind(rep(1, length(temp)), t.c))

  # finds positions of coefficients of temperature components
  idx.tc <- grep("tc.\\d", names(coef(reg.model)), perl = TRUE)

  # finds positions of coefficients of interaction terms
  idx.inter <- grep("tc.\\d:", names(coef(reg.model)), perl = TRUE)
  
  # exclude interaction terms to get only temperature components
  idx.tc <- setdiff(idx.tc, idx.inter)
  
  # gets coefficient values for temp. componentes (tc)
  coef.values <- c(coef(reg.model)[c(1, idx.tc)])
  
  if(!is.null(inter.term)) {
    # create tc * (mean(dew.point | tc))
    mean.value <- mean.value.by.tc(data.set, inter.term, temp.breaks)$mean
    interaction.columns <- t(apply(X=t.c, MARGIN = 1, 
                                   FUN = function(x,y) return(x * y),
                                   mean.value))
    colnames(interaction.columns) <- paste0(colnames(interaction.columns),
                                            ":", inter.term) 
    X.2 <- cbind(X.2, interaction.columns)
    coef.values <- c(coef.values, coef(reg.model)[c(idx.inter)])
  }
  
  coef.values[is.na(coef.values)] <- 0
  
  # computes matrix multiplication y = X * beta
  y.2 <- X.2 %*% coef.values
  
  if (!normalize.values) {
    # shift curve vertically by average value of annual fixed effects
    idx.annual <- grep("factor", names(coef(reg.model)))
    # shift curve vertically by average value of hourly fixed effects
    idx.hour <- grep("hour.of.day", names(coef(reg.model)))
    y.2 <- y.2 + mean(coef(reg.model)[idx.annual], na.rm = TRUE) + 
      mean(coef(reg.model)[idx.hour], na.rm = TRUE)
  }else{
    a <- data.frame(type=names(coef(reg.model)), 
                    values=coef(reg.model), row.names = NULL, 
                    stringsAsFactors = FALSE)
    annual.a <- a[grep('factor', a$type), ]
    names(annual.a) <- c('annual.type', 'annual.value')
    hour.a <- a[grep('hour.of.day', a$type), ]
    names(hour.a) <- c('hour.type', 'hour.value')

    data.set <- data.set %>% 
      mutate(hour.type=paste0('hour.of.day',hour.of.day,
                              ':type.day',type.day,
                              ':season', season),
             annual.type=paste0('factor(',year,')')) %>%
      left_join(annual.a, by='annual.type') %>%
      left_join(hour.a, by='hour.type') %>%
      mutate(norm.value=load-ifelse(is.na(hour.value), 0, hour.value)-
               ifelse(is.na(annual.value), 0, annual.value)) %>%
      mutate(temp=temp, load=norm.value)
    #data.set$load <- data.set$load - coef(reg.model)[1]
  }
  
  # **** plots resulting curve of load vs temperature ---
  
  df.plot.out1 <- data.frame(model=character(), x=numeric(), y=numeric())
  df.plot.out1 <- rbind(df.plot.out1, data.frame(model = "model", 
                                                 x = temp, 
                                                 y = y.2/1e3))
  g <- ggplot()
  g <- g + geom_point(data=data.set, aes(x=temp, y=load/1e3),
                      col = rgb(211, 211, 211, 20,
                                maxColorValue = 255), size=0.4)
  g <- g + theme_bw() + 
    xlab(expression(paste("Temperature (", degree, "C)")))
  
  if (normalize.values){
    g <- g + ylab("Residualized load (GW)")
  }else{
    g <- g + ylab("load (GW)")
  }
    
  g <- g + geom_line(data = df.plot.out1, aes(x=x, y=y, colour = model))
  g <- g + guides(colour=FALSE)
  g <- g + theme()
  #panel.border=element_blank(),
  #axis.line=element_line()
  #g <- g + geom_point(data = df.plot.out1, aes(x=x, y=y, colour = model)) 
  #g <- g + guides(colour=guide_legend(title=NULL))
  #g <- g + theme(legend.position=c(0.98,0.98), legend.justification=c(1,0))
  
  if (add.marg.hist) {
    
    # change margins of plot
#    g <- g + theme(plot.margin=unit(c(0,0,1,1), "cm"))
    
    #get info from main plot
    gb <- ggplot_build(g)
    
    y.range <- gb$layout$panel_ranges[[1]]$y.range
    x.range <- gb$layout$panel_ranges[[1]]$x.range
    
    y.major <- gb$layout$panel_ranges[[1]]$y.major_source
    y.minor <- gb$layout$panel_ranges[[1]]$y.minor_source
    
    x.major <- gb$layout$panel_ranges[[1]]$x.major_source
    x.minor <- gb$layout$panel_ranges[[1]]$x.minor_source
    
    # create marginal histograms to top and right
    x <- data.set$temp[!is.na(data.set$temp)]
    breaks <- pretty(x.range, n = nclass.FD(x), min.n = 1)
    bwidth <- breaks[2]-breaks[1]
    pTop <- ggplot(data.set, aes(x = temp)) +
      geom_histogram(binwidth=bwidth, na.rm = TRUE, 
                     fill='gray75', colour='gray40') +
      theme_bw() +
      scale_x_continuous(breaks = x.major, minor_breaks = x.minor,
                         limits = x.range, expand = c(0, 0)) +
      scale_y_continuous(breaks=NULL, minor_breaks = NULL) +
      theme(axis.title = element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            plot.margin=unit(c(5.5, 5.5, -1.5, 5.5), "points"),
            panel.border = element_blank()
            )
    
    x <- data.set$load[!is.na(data.set$load)]/1e3
    breaks <- pretty(y.range, n = nclass.FD(x), min.n = 1)
    bwidth <- breaks[2]-breaks[1]
    pRight <- ggplot(data.set, aes(x = load/1e3)) +
      geom_histogram(binwidth=bwidth, na.rm = TRUE,
                     fill='gray75', colour='gray40') + 
      scale_x_continuous(breaks = y.major, minor_breaks = y.minor,
                         limits = y.range, expand = c(0,0)) +
      scale_y_continuous(breaks=NULL, minor_breaks = NULL) +
      coord_flip() + theme_bw() +
      theme(axis.title = element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            plot.margin=unit(c(5.5, 5.5, 5.5, -2.75), "points"),
            panel.border = element_blank()
            )
  
    pEmpty <- textGrob("")
    
    gmain <- ggplotGrob(g)
    gtop <- ggplotGrob(pTop)
    gright <- ggplotGrob(pRight)
    
    maxWidth <- grid::unit.pmax(gmain$widths[2:5], gtop$widths[2:5])
    maxHeight <- grid::unit.pmax(gmain$heights[2:9], gright$heights[2:9])
    
    gmain$widths[2:5] <- as.list(maxWidth)
    gtop$widths[2:5] <- as.list(maxWidth)

    gmain$heights[2:9] <- as.list(maxHeight)
    gright$heights[2:9] <- as.list(maxHeight)
    
    g <- arrangeGrob(gtop, pEmpty, gmain, gright, ncol=2, nrow=2, 
                     widths = c(7, 1), heights = c(1, 7))
    
    #g <- arrangeGrob(gtop, pEmpty, gmain, gright, ncol=2, nrow=2, 
    #                 widths = c(7, 1), heights = c(1, 7))
    
    #grid.arrange(g)
  }
  
  if(!dir.exists(dir.out)) {
    dir.create(dir.out)
  }
  
  name.out <- "plotTempLoad.png"
  if (substring(dir.out, first = nchar(dir.out)) != "/") {
    dir.out <- paste0(dir.out, "/")
  }
  
  png(file = paste0(dir.out, name.out), width = 1500, 
      height = 1500, res = 300)
  grid.arrange(g)
  dev.off()
  
  return(df.plot.out1)
}

plot3d.temp.load <- function(reg.model, data.set,
                             temp.breaks = c(-10, 0, 10, 20, 30),
                             temp = c(-15:40)){
  # Creates 3D plot of temperature, humidity, estimated load
  # (complete function)
  #
  # Args:
  #   reg.model: linear regression model object
  #   data.set: original data set used in regression model (before computing
  #             temperature components)
  #   temp.breaks:  breakpoints of temperature intervals in stepwise linear 
  #                 model
  #   temp: temperature values in X axis
  #
  # Returns:
  #   nothing
  
  library(rgl)
  
  # compute temperature components and create X matrix  
  t.c <- createTempComponents(temp, temp.breaks = temp.breaks)
  
  t.c.aux <- cbind(temp,t.c)
  
  temp.ext <- rep(temp,rep(length(temp), length(temp)))
  dp.ext <- rep(temp, length(temp))
  
  data.aux <- data.frame(temp=temp.ext, dew.point=dp.ext)
  rm(temp.ext, dp.ext)
  
  data.aux$dew.point[data.aux$dew.point > data.aux$temp] <- NA
  
  data.aux <- merge(x=data.aux, y=t.c.aux, by="temp")
  
  data.aux <- cbind(data.aux,
                    t(apply(X=data.aux, MARGIN = 1, 
                            FUN = {function(x) {
                              nc = length(x)
                              return(x[3:nc]*x[2])}})))
  
  # fix names of interaction columns in data.aux
  n1 <- ncol(data.aux) # number of columns in data aux
  n2 <- ncol(t.c) # number of columns in t.c
  names(data.aux)[(n1-n2+1):n1] <- paste0("dp:", names(t.c))
  
  X.2 <- as.matrix(cbind(rep(1, nrow(data.aux)), 
                         data.aux[ ,c(3:n1)]))
  
  # finds positions of coefficients of temperature components
  idx.tc <- grep("tc.\\d", names(coef(reg.model)), perl = TRUE)
  
  # finds positions of coefficients of interaction terms
  idx.inter <- grep("tc.\\d:", names(coef(reg.model)), perl = TRUE)
  
  # exclude interaction terms to get only temperature components
  idx.tc <- setdiff(idx.tc, idx.inter)
  
  # gets coefficient values for temp. componentes (tc) & interactions
  coef.values <- c(coef(reg.model)[c(1, idx.tc)])
  
  coef.values <- c(coef.values, coef(reg.model)[c(idx.inter)])
  
  coef.values[is.na(coef.values)] <- 0
  
  # computes matrix multiplication y = X * beta
  y.2 <- X.2 %*% coef.values
  
  # **** plots resulting curve of load vs temperature ---
  
  df.plot.out1 <- data.frame(temp=data.aux$temp, 
                             dew.point=data.aux$dew.point, 
                             load=y.2)
  m <- df2mat(df.plot.out1)
  
  sample.data <- data.set[sample(nrow(data.set), 20000), ]
  
  # set size of windows (512 x 512)
  par3d("windowRect" = c(0, 0, 512, 512))
  
  # plot prediction surface
  plot3d(df.plot.out1$temp, df.plot.out1$dew.point, df.plot.out1$load, 
         type = "p", size = 0.2, alpha = 0.2, lit=FALSE, 
         xlab = "Temperature", ylab = "Dew Point", zlab = "Load")
  surface3d(m$temp, m$dew.point, t(m$load), alpha=0.4, 
            front="line", back="line")
  
  # Add some observed points
  points3d(sample.data$temperature, sample.data$dew.point, sample.data$load, 
           alpha=0.3, col="red", size=1)
  
  movie3d(spin3d(axis=c(0,0,1), rpm=4), duration=10, fps=24)
  
  # Add line segments showing the error
  # segments3d(data.set$temperature, data.set$temperature),
  #            interleave(data.set$dew.point, data.set$dew.point),
  #            interleave(data.set$load, m$pred_mpg),
  #            alpha=0.4, col="red")
  
  
  #png(file = "./out/3dSurfaceplot.png", width = 1500, height = 1500, res = 300)
  #wireframe(load ~ temp*dew.point, 
  #          data = df.plot.out1, colorkey = TRUE, drape = TRUE)
  #dev.off()
  
  
}

scatter.plot.temp.load <- function(data.set) {
  # Creates scatter plot of temperature in X axis and load in Y axis
  # (only observed values)
  #
  # Args:
  #   data.set: original data set
  #
  # Returns:
  #   ggplot object
  
  g <- ggplot() + 
    geom_point(data=data.set, aes(x=temp, y=load/1e3), 
               col = rgb(100, 100, 100, 20, maxColorValue = 255), size=0.4) + 
    theme_bw() + xlab(expression(paste("Temperature (", degree, "C)"))) +
    ylab("Hourly load (GW)")
  return(g)
}

grid.plot.temp.load <- function(data.set) {
  # Creates grid of scatter plots temp x load in each season
  # (only observed values)
  #
  # Args:
  #   data.set: original data set
  #
  # Returns:
  #   ggplot object
  
  library(grid)
  library(gridExtra)
  
  seasons <- unique(data.set$season)
  gg.list <- list()
  n.seasons <- length(seasons)
  x.range <- range(data.set$temperature, na.rm = TRUE)
  y.range <- range(data.set$load/1e3, na.rm = TRUE)
  
  for (i in 1:n.seasons) {
    data.set.season <- filter(data.set, season == seasons[i])
    g <- ggplot(data = data.set.season)
    g <- g + geom_point(aes(x=temperature, y=load/1e3), 
                        col = rgb(100, 100, 100, 20, maxColorValue = 255), 
                        size=0.4) + theme_bw()
    g <- g + annotate("text", 
                      x=mean(range(x.range)), 
                      y=Inf,
                      label=seasons[i], vjust=1.5, size=4)
    g <- g + xlim(x.range) + ylim(y.range)
    if (i %in% c(1,3)) {
      g <- g + ylab("load (GW)") 
    } else {
      g <- g + theme(axis.title.y=element_blank())
    }
    if (i %in% c(3, 4)) {
      g <- g + xlab(expression(paste("Temperature (", degree, "C)")))
    } else {
      g <- g + theme(axis.title.x=element_blank())
    }
    gg.list[[seasons[i]]] <- g
  }
  gg.out <- arrangeGrob(grobs = gg.list)
  grid.arrange(gg.out)
  return(gg.out)
}

df2mat <- function(p, xvar = NULL, yvar = NULL, zvar = NULL) {
  # Convert long-style data frame with x, y, and z vars into a list
  # with x and y as row/column values, and z as a matrix.
  # (@ R Graphics Cookbook by Winston Chang - page 286)
  #
  # Args:
  #   p: (data frame) data set (long format)
  #   xvar, yvar, zvar: (string) names of variables
  #
  # Returns:
  #   list
  
  if (is.null(xvar)) xvar <- names(p)[1]
  if (is.null(yvar)) yvar <- names(p)[2]
  if (is.null(zvar)) zvar <- names(p)[3]
  x <- unique(p[[xvar]])
  x <- x[!is.na(x)] # drop NA
  
  y <- unique(p[[yvar]])
  y <- y[!is.na(y)] # drop NA
  
  z <- matrix(p[[zvar]], nrow = length(y), ncol = length(x))
  m <- list(x, y, z)
  
  names(m) <- c(xvar, yvar, zvar)
  m
}

interleave <- function(v1, v2) {
  # Function to interleave the elements of two vectors
  # (@ R Graphics Cookbook by Winston Chang - page 286)
  return(as.vector(rbind(v1,v2)))
} 

predictgrid <- function(model, xvar, yvar, zvar, res = 16, type = NULL) {
  # Given a model, predict zvar from xvar and yvar
  # Defaults to range of x and y variables, and a 16x16 grid
  # (@ R Graphics Cookbook by Winston Chang - page 286)
  #
  # Args:
  #   model: (model object)
  #   xvar, yvar, zvar: (string) names of variables
  #   res:
  #   type:
  #
  # Returns:
  #   list
  
  # Find the range of the predictor variable. This works for lm and glm
  # and some others, but may require customization for others.
  xrange <- range(model$model[[xvar]])
  yrange <- range(model$model[[yvar]])
  newdata <- expand.grid(x = seq(xrange[1], xrange[2], length.out = res),
                         y = seq(yrange[1], yrange[2], length.out = res))
  names(newdata) <- c(xvar, yvar)
  newdata[[zvar]] <- predict(model, newdata = newdata, type = type)
  return(newdata)
}

plot.intraday.curve <- function(lm.out, out.of.sample = TRUE) {

  if(out.of.sample) {
    if ("performance" %in% names(lm.out)) {
      df.1 <- data.frame(hour.of.day=lm.out$performance$test.data$hour.of.day,
                         season=lm.out$performance$test.data$season,
                         load=lm.out$performance$test.data$load/1e3)
      df.2 <- data.frame(hour.of.day=lm.out$performance$test.data$hour.of.day,
                         season=lm.out$performance$test.data$season,
                         load=lm.out$performance$y.hat$y.hat/1e3)
    }else {
      cat("\nNo out of sample prediction!!! Creating plot with fitted data.\n")
      df.1 <- data.frame(hour.of.day=lm.out$model$model$hour.of.day,
                         season=lm.out$mode$modell$season,
                         load=lm.out$model$model$load/1e3)
      df.2 <- data.frame(hour.of.day=lm.out$model$model$hour.of.day,
                         season=lm.out$model$model$season,
                         load=lm.out$model$fitted.values/1e3)
    }
  } else {
    df.1 <- data.frame(hour.of.day=lm.out$model$model$hour.of.day,
                       season=lm.out$model$model$season,
                       load=lm.out$model$model$load/1e3)
    df.2 <- data.frame(hour.of.day=lm.out$model$model$hour.of.day,
                       season=lm.out$model$model$season,
                       load=lm.out$model$fitted.values/1e3)
  }
    
  
  # plot estimated and actual mean intraday load curve for each season
  df.real <- as.data.frame(df.1 %>% group_by(season, hour.of.day) %>%
    summarize(mean.value = mean(load, na.rm=TRUE)))
  df.real <- df.real %>% filter(season %in% c("Summer"))
  
  df.real$hour.of.day <- as.numeric(as.character(df.real$hour.of.day))
  df.real$type <- factor(x = "Observed", levels = c("Observed", "Predicted"))
  
  df.pred <- df.2
  df.pred <- df.pred %>% filter(season %in% c("Summer", "Winter"))
  
  df.pred <- as.data.frame(df.pred %>% group_by(season, hour.of.day) %>%
                             summarize(mean.value = mean(load, na.rm=TRUE)))
  df.pred$hour.of.day <- as.numeric(as.character(df.pred$hour.of.day))
  df.pred$type <- factor(x = "Predicted", levels = c("Observed", "Predicted"))
  
  df.total <- rbind(df.real, df.pred)
  
  df.range <- df.2 %>% group_by(season, hour.of.day) %>%
    summarize(qt1 = quantile(load, probs = c(0.025) ,na.rm=TRUE),
              qt2 = quantile(load, probs = c(0.975) ,na.rm=TRUE)) %>%
    filter(season == "Summer")
  
  df.range$hour.of.day <- as.numeric(as.character(df.range$hour.of.day))
  
  g <- ggplot()
  g <- g + geom_line(data=df.total, aes(x=hour.of.day, y=mean.value, 
                                        linetype=type))
  g <- g + scale_linetype_manual(values = c("solid", "dashed"))
  g <- g + scale_colour_manual(values = c("black", "gray"))
  
  g <- g + theme_bw() + xlab("hour of day") + ylab("Load (GW)")
  g <- g + theme(legend.position=c(0.95,0.05), legend.justification=c(1,0))
  g <- g + guides(colour=guide_legend(title=NULL), 
                  linetype=guide_legend(title=NULL))
  
  g <- g + geom_ribbon(data=df.range, 
                       aes(x=hour.of.day, ymin=qt1, ymax=qt2),
                       alpha=0.2)
  
  g <- g + facet_wrap(~ season, ncol=2)
  
  y.range <- layer_scales(g)$y$range$range
  g <- g + ylim(0, y.range[2])
  
  # if (substring(dir.out, first = nchar(dir.out)) != "/") {
  #   dir.out <- paste0(dir.out, "/")
  # }
  # png(file = paste0(dir.out, "plotintraday.png"), width = 2000, 
  #     height = 1333, res = 300)
  # print(g)
  # dev.off()
  
  return(g)
}

plot.intraday.curve.2 <- function(lm.out) {
  # lm.out: list with resulting linear model
  
  df.1 <- data.frame(hour.of.day=lm.out$performance$test.data$hour.of.day,
                     season=lm.out$performance$test.data$season,
                     load=lm.out$performance$test.data$load/1e3)
  df.2 <- data.frame(hour.of.day=lm.out$performance$test.data$hour.of.day,
                     season=lm.out$performance$test.data$season)
  
  y.hat <- predict.lm(lm.out$model, newdata = lm.out$performance$test.data, 
                      interval = "prediction")
  y.hat <- y.hat/1e3
  
  df.2 <- cbind(df.2, y.hat)
  
  # plot estimated and actual mean intraday load curve for each season
  df.real <- as.data.frame(df.1 %>% group_by(season, hour.of.day) %>%
                             summarize(mean.value = mean(load, na.rm=TRUE)))
  df.real <- df.real %>% filter(season %in% c("Summer"))
  
  df.real$hour.of.day <- as.numeric(as.character(df.real$hour.of.day))
  df.real$type <- factor(x = "Observed", levels = c("Observed", "Predicted"))
  
  df.pred <- df.2
  df.pred <- df.pred %>% filter(season %in% c("Summer"))
  
  df.pred <- as.data.frame(df.pred %>% group_by(season, hour.of.day) %>%
                             summarize(mean.value = mean(fit, na.rm=TRUE)))
  df.pred$hour.of.day <- as.numeric(as.character(df.pred$hour.of.day))
  df.pred$type <- factor(x = "Predicted", levels = c("Observed", "Predicted"))

  df.pred.2 <- as.data.frame(df.2 %>% filter(season %in% c("Summer")) %>% 
                               group_by(season, hour.of.day) %>%
                               summarize(ub = mean(upr, na.rm=TRUE),
                                         lb = mean(lwr, na.rm=TRUE)))
  df.pred.2$hour.of.day <- as.numeric(as.character(df.pred.2$hour.of.day))
  df.pred.2$type <- factor(x = "Predicted", levels = c("Observed", "Predicted"))
  
  df.total <- rbind(df.real, df.pred)
  
  # df.range <- df.2 %>% group_by(season, hour.of.day) %>%
  #   summarize(qt1 = quantile(load, probs = c(0.025) ,na.rm=TRUE),
  #             qt2 = quantile(load, probs = c(0.975) ,na.rm=TRUE)) %>%
  #   filter(season == "Summer")
  # df.range$hour.of.day <- as.numeric(as.character(df.range$hour.of.day))
  
  g <- ggplot()
  g <- g + geom_line(data=df.total, aes(x=hour.of.day, y=mean.value, 
                                        linetype=type))
  g <- g + scale_linetype_manual(values = c("solid", "dashed"))
  g <- g + scale_colour_manual(values = c("black", "gray"))
  
  g <- g + theme_bw() + xlab("hour of day") + ylab("Load (GW)")
  g <- g + theme(legend.position=c(0.95,0.05), legend.justification=c(1,0))
  g <- g + guides(colour=guide_legend(title=NULL), 
                  linetype=guide_legend(title=NULL))
#  g <- g + geom_point(data=df.pred, 
#                         aes(x=hour.of.day, y=mean.value), size=1)
  g <- g + geom_ribbon(data=df.pred.2, 
                         aes(ymin=lb, ymax=ub, x=hour.of.day), alpha=0.2)
#  g <- g + geom_ribbon(data=df.range, 
#                       aes(x=hour.of.day, ymin=qt1, ymax=qt2),
#                       alpha=0.2)
  
  y.range <- layer_scales(g)$y$range$range
  g <- g + ylim(0, y.range[2])
  #g
  return(g)
}
