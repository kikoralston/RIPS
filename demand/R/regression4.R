
library(car)
library(ggplot2)

compute.performance <- function(test.data, reg.model, annual.model) {
  # Computes the performance of forecasting power of regression model
  #
  # Args:
  #   test.data:  data frame with testing data (subset of total data)
  #   reg.model: fitted regression hourly model
  #   annual.model: fitted regression annual model
  #
  # Returns:
  #   list.out: list with prediction values, RMSE and MAPE
  
  y.hat <- predict.lm(reg.model, newdata = test.data)

  if (!is.null(annual.model)) {
    # subtract intercept
    y.hat <- y.hat - coef(reg.model)[1]
    
    annual.test.df <- get.annual.df(test.data)
    annual.fe.hat <- predict(annual.model, newdata = annual.test.df)
    
    name.response <- as.character(attr(annual.model$terms, "variables"))[2]
    if (name.response == "log(fix.eff)") {
      # make log adjustment
      MSE <- summary(annual.model)$sigma^2
      annual.fe.hat <- exp(MSE/2)*exp(annual.fe.hat)
    }
    
    y.hat <- y.hat + annual.fe.hat
  }
  
  name.response <- as.character(attr(reg.model$terms, "variables"))[2]
  if (name.response == "log(load)") {
    # make log adjustment
    MSE <- summary(reg.model)$sigma^2
    y.hat <- exp(MSE/2)*exp(y.hat)
  }
  
  rmse <- sqrt(mean((test.data$load - y.hat)^2, na.rm = TRUE))
  mape <- mean(abs((test.data$load - y.hat)/test.data$load), na.rm = TRUE)
  
  y.hat <- data.frame(test.data, y.hat=y.hat)
  
  list.out <- list()
  list.out$rmse <- rmse
  list.out$mape <- mape
  list.out$test.data <- test.data
  list.out$y.hat <- y.hat
  
  return(list.out)
}

parseRegressionFormula <- function(lm.formula) {
  
  formula.char <- deparse(lm.formula)
  
  x <- regexpr("\\^\\d\\)", formula.char, perl=TRUE)
  pos <- as.numeric(x)
  exp.string <- substring(formula.char, pos, pos+1)
    
  formula.char <- gsub("I\\(", "", formula.char, perl = TRUE)
  formula.char <- gsub("\\^\\d\\)", paste0("$", exp.string, "$"), 
                       formula.char, perl = TRUE)
  formula.char <- gsub("~", "$\\\\sim$", formula.char, perl = TRUE)
  
  #formula.char <- paste0("$", formula.char, "$")
  return(formula.char)
}

writeRegressionSummary <- function(reg.model, idx.coeff=NULL){
  # Creates table with output of regression (as in Anglin and Gencay (1996))
  #
  # Args:
  #   reg.model: linear regression object
  #   idx.coef:  vector if indexes of coefficients to show
  #
  # Returns:
  #   Matrix object with output table to be exported to csv

  # summary
  reg.sum <- summary(reg.model)

  if(!is.null(idx.coeff)) {
    # coefficients
    out.matrix <- signif(reg.sum$coefficients[idx.coeff, ], 3)
  } else {
    # coefficients
    out.matrix <- signif(reg.sum$coefficients, 3)
  }
  
  out.matrix[(out.matrix[ , 4] < 1e-10), 4] <- "$<10^{-10}$"
  
  # names of variables
  names.row <- row.names(out.matrix)
  
  # check for polynomial terms and correct them
  idx.poly <- grep("I\\(", names.row, perl = TRUE)
  for (i in idx.poly) {
    x <- regexpr("\\^\\d",row.names(out.matrix)[i], perl=TRUE)
    pos <- as.numeric(x)
    exp.string <- substring(row.names(out.matrix)[i], pos, pos+1)
    
    names.row[i] <- gsub("\\^\\d", "",names.row[i], perl = TRUE)
    names.row[i] <- gsub("I\\(", "",names.row[i], perl = TRUE)
    names.row[i] <- gsub("\\)", "",names.row[i], perl = TRUE)
    
    names.row[i] <- paste0("$\\text{", names.row[i], "}",exp.string, "$")
  }
  row.names(out.matrix) <- names.row
  
  # add row with SSR
  #ssr <- sum(reg.sum$residuals^2)
  #new.row <- matrix(data=c(round(ssr,2), "", "", ""), 
  #                  nrow=1, ncol=4)
  #rownames(new.row) <- "SSR"
  #out.matrix <- rbind(out.matrix, new.row)
  
  # add row with R^2
  new.row <- matrix(data=c(round(reg.sum$r.squared,2), "", "", ""), 
                    nrow=1, ncol=4)
  rownames(new.row) <- "$R^2$"
  out.matrix <- rbind(out.matrix, new.row)
  
  # add row with adjusted R^2
  new.row <- matrix(data=c(round(reg.sum$adj.r.squared,2), "", "", ""), 
                    nrow=1, ncol=4)
  rownames(new.row) <- "$R^2_{adj}$"
  out.matrix <- rbind(out.matrix, new.row)
  
  # add row with Number of observations
  new.row <- matrix(data=c(length(reg.sum$residuals), "", "", ""), 
                    nrow=1, ncol=4)
  rownames(new.row) <- "Number of observations"
  out.matrix <- rbind(out.matrix, new.row)
  
  return(out.matrix)
}

getAnnualFixedEffects <- function(reg.model, flag.plot = FALSE, 
                                  folder.out = ".out/") {
  
  # get indexes of annual fixed effects in vector of betas
  idx.year <- grep("year", names(coef(reg.model)))
  
  # annual fixed effects (\beta_0 + \gamma_y)
  annual.fe <- c(coef(reg.model)[1], 
                 coef(reg.model)[1] + coef(reg.model)[idx.year])
  
  years.data <- as.numeric(levels(reg.model$model$`factor(year)`))
  
  df.fe <- data.frame(year=years.data, fix.eff=annual.fe)
  row.names(df.fe) <- NULL
  
  # plot time series of annual fixed effects
  if (flag.plot) {
    png(filename = paste0(folder.out, "annual_fe.png"), width=3000, 
        height=2000, res=300)
    par(mar = c(5, 4, 1, 1) + 0.1)
    plot(df.fe, type="b", ylab="Annual Fixed Effects (GW)",
         xlab="Year")
    dev.off()
  }
  
  return(df.fe)
}

regression.model.1 <- function(data.set, 
                               temp.breaks = c(-10, 0, 10, 20, 30),
                               test.year = NULL,
                               type.humidity = NULL,
                               log.load = FALSE,
                               idx.annual = 1,
                               formula.regression=NULL) {
  # Fits regression model to training set and makes prediction using test set
  #
  # Args:
  #   data.set:  data frame with data set to be used in regression
  #   temp.breaks:  breakpoints of temperature intervals in stepwise linear 
  #                 model  
  #   test.year:  years that should be used as test set. If NULL whole set is 
  #               used as traning set and there is no prediction
  #   type.humidity:  (string) Type of humidity metric to use. 
  #                   Accepted Values: "dp" (dew point) or 
  #                   "rh" (relative humidity) or NULL (use none). 
  #   log.load: (boolean) TRUE if log transformation should be applied to load
  #   idx.annual: (integer) index of annual model to be used.
  #               (see function "list.annual.models()")
  #   formula.regression: (string) sets regression formula
  #               directly (overrides previous arguments != NULL)
  #
  # Returns:
  #   list.out: list with results
  
  # prepare data and fit model

  years.data <- unique(data.set$year)
  years.data <- years.data[order(years.data)]
  
  # temperature components
  temp.components <- createTempComponents(data.set$temp, temp.breaks)
  data.set <- cbind(data.set, temp.components)
  rm(temp.components)

  if (is.null(formula.regression)) {
    # test inputs
    if (!is.null(type.humidity)) {
      if (!(type.humidity %in% c("dp", "rh"))) {
        stop("Argument type.humidity invalid! Must be \"dp\" or \"rh\" ")
      }
    }
    
    # creates regression equation
    if(log.load) {
      equation.string <- "log(load) ~ hour.of.day:type.day:season + factor(year)"
    } else {
      equation.string <- "load ~ hour.of.day:type.day:season + factor(year)"
    }
    
    n.cp <- length(temp.breaks) + 1  # number of components
    names.cp <- paste0("tc.", c(1:n.cp))
    
    equation.string <- paste0(equation.string, " + ",
                              paste0(names.cp, collapse = " + "))
    
    #name.humidity <- NULL
    if (!is.null(type.humidity)) {
      if (type.humidity == "dp") {
        names.interaction <- paste0(names.cp,":dp")
        #name.humidity <- "dp"
      } else {
        names.interaction <- paste0(names.cp,":rh")
        #name.humidity <- "rh"
      }
      # add interaction term between temperature components and dew point
      equation.string <- paste0(equation.string, " + ",
                                paste0(names.interaction, collapse = " + "))
    }
  } else {
    equation.string <- formula.regression
  }
  
  cat("\n--------------------------------------------\n")
  cat(equation.string)
  cat("\n--------------------------------------------\n")

  # separate training data 
  if(!is.null(test.year)) {
    df.train <- dplyr::filter(data.set, !(year %in% test.year))
  } else {
    df.train <- data.set
  }
  
  # fits hourly model to training data set
  reg.model <- lm(formula = as.formula(equation.string), 
                  data = df.train)
  
  if (!is.null(idx.annual)) {
    # fits annual model
    annual.model <- fit.annual.model(df.train, reg.model, idx.annual)
  } else {
    # don't use annual model (repeat previous value of \gamma_a)
    annual.model <- NULL
  }
  
  # computes performance of prediction (only if a test set is defined)
  pred.perf <- NULL
  gg.plot <- NULL
  if (!is.null(test.year)) {
    df.test <- dplyr::filter(data.set, year %in% test.year)
    
    if (is.null(idx.annual)) {
      # repeat previous year
      df.test$year <- years.data[which(years.data == min(test.year)) - 1]
    } else { 
      # set to initial year (\gamma_a == 0)
      df.test$year <- min(df.train$year, na.rm = TRUE)
    }
    pred.perf <- compute.performance(df.test, reg.model, annual.model)
  }

  # prepare list with results to output
  
  list.out <- list()
  list.out$model <- reg.model
  list.out$temp.bins <- c(paste0('(-Inf,', min(temp.breaks), ']'), 
                          levels(cut(1, temp.breaks)), 
                          paste0('(', max(temp.breaks), ',+Inf)'))
  list.out$annual.model <- annual.model
  list.out$performance <- pred.perf
  # list.out$gg.plot <- gg.plot
  
  return(list.out)
}

leave.one.year.out.cv <- function(data.set,
                                  temp.breaks = c(-10, 0, 10, 20, 30),
                                  type.humidity = NULL,
                                  log.load = FALSE,
                                  idx.annual=1) {
  # leave one year out cross validation
  #
  # Args:
  #   data.set:  data frame with data set to be used in regression
  #   temp.breaks:  breakpoints of temperature intervals in stepwise linear 
  #                 model  
  #   type.humidity:  (string) Type of humidity metric to use. 
  #                   Accepted Values: "dp" (dew point) or 
  #                   "rh" (relative humidity) or NULL (use none). 
  #   log.load: (boolean) TRUE if log transformation should be applied to load
  #
  # Returns:
  #   list.out: list with results
  
  years.data <- unique(data.set$year)
  years.data <- years.data[!is.na(years.data)]
  
  cat("\n------------------------------\n")
  cat("    Running Cross Validation    ")
  cat("\n------------------------------\n")
  #cat("\n\n")
  
  cv.results <- list()
  
  for (i in 1:length(years.data)) {
    cat("\nTest year: ", years.data[i], "\n")
    reg.result <- regression.model.1(data.set = data.set,
                                     test.year = years.data[i],
                                     temp.breaks = temp.breaks,
                                     type.humidity = type.humidity,
                                     log.load = log.load,
                                     idx.annual = idx.annual)
    
    cv.results[[i]] <- list(rmse=reg.result$performance$rmse,
                            mape=reg.result$performance$mape)
  }
  
  return(cv.results)
}

get.annual.df <- function(data.set) {
  
  annual.data <- data.set %>% group_by(year) %>% 
    summarise(TN.GDP=mean(TN.GDP, na.rm=TRUE), 
              TN.POP=mean(TN.POP, na.rm=TRUE),
              us.energyprod=mean(us.energyprod,na.rm=TRUE),
              ng.price.real=mean(ng.price.real,na.rm=TRUE))
  annual.data <- data.frame(year=annual.data$year,
                            fix.eff=NA,
                            gdp.tn=annual.data$TN.GDP,
                            pop.tn=annual.data$TN.POP,
                            gdp.cap.tn=annual.data$TN.GDP*1e3/annual.data$TN.POP,
                            us.energyprod=annual.data$us.energyprod,
                            ng.price.real=annual.data$ng.price.real)
  return(annual.data)
}

fit.annual.model <- function(data.set, hourly.model, idx.annual = 1) {
  
  annual.data <- get.annual.df(data.set)
  annual.fe <- getAnnualFixedEffects(hourly.model)
  
  for (irow in 1:nrow(annual.data)) {
    irow2 <- which(annual.fe$year == annual.data$year[irow])
    if(length(irow2) > 0) {
      annual.data$fix.eff[irow] <- annual.fe$fix.eff[irow2] 
    }
  }
  
  annual.model <- list.annual.models()[[idx.annual]]
  
  lm.annual <- lm(annual.model, data=annual.data)
  
  return(lm.annual)
}

list.annual.models <- function() {
  
  annual.models <- list()
  annual.models[[1]] <- as.formula("fix.eff ~ gdp.cap.tn + I(gdp.cap.tn^2)")
  annual.models[[2]] <- as.formula("log(fix.eff) ~ log(gdp.cap.tn)")
  annual.models[[3]] <- as.formula("fix.eff ~ log(gdp.cap.tn)")
  annual.models[[4]] <- as.formula("fix.eff ~ gdp.cap.tn + I(gdp.cap.tn^2) + gdp.tn")
  annual.models[[5]] <- as.formula("fix.eff ~ gdp.cap.tn + 
                                   I(gdp.cap.tn^2) + pop.tn")
  annual.models[[6]] <- as.formula("fix.eff ~ gdp.cap.tn + I(gdp.cap.tn^2) + 
                                   us.energyprod")
  annual.models[[7]] <- as.formula("fix.eff ~ gdp.cap.tn + I(gdp.cap.tn^2) + 
                                   log(us.energyprod)")
  annual.models[[8]] <- as.formula("log(fix.eff) ~ log(gdp.cap.tn) + 
                                   log(us.energyprod)")
  annual.models[[9]] <- as.formula("fix.eff ~ gdp.cap.tn + ng.price.real")
  annual.models[[10]] <- as.formula("fix.eff ~ log(gdp.cap.tn) + log(ng.price.real)")
  annual.models[[11]] <- as.formula("fix.eff ~ gdp.cap.tn + 
                                   I(gdp.cap.tn^2) + ng.price.real")
  return(annual.models)
}

# ------- GARBAGE -----
#
#
# seasons <- c("Summer", "Fall", "Winter", "Spring")
# 
# models <- list()
# 
# for (i in 1:length(seasons)) {
#   df.train.season <- subset(df.final, season == seasons[i] & year(time) < 2014)
#   models[[i]] <- lm(load ~ hour.of.day + type.day + factor(year(time)) +
#                       tc.1 + tc.2 + tc.3 + tc.4 + tc.5 + tc.6 +
#                       tc.1:dew.point + tc.2:dew.point + tc.3:dew.point + 
#                       tc.4:dew.point + tc.5:dew.point + tc.6:dew.point, 
#                data = df.train.season)
# }
# 
# names(models) <- seasons
#df.train.season <- subset(df.final, season == seasons[i] & year(time) < 2014)
#mean.dp <- mean(df.train.season$dew.point, na.rm = TRUE)
# idx.year <- grep("year", names(coef(reg.model)))
# mean.year <- mean(coef(reg.model)[idx.year])

# i <- 4
# df.test.season <- subset(df.final, season == seasons[i] & year(time) == 2014)
# year(df.test.season$time) <- 2013
# 
# 
# # plot actual hourly load in test data set
# g <- ggplot()
# g <- g + geom_line(data=df.test.season, aes(x=time, y=load), 
#                    colour = "darkgray")
# g <- g + theme_bw() + scale_x_datetime(breaks=ticks.plot, date_labels = "%b")
# g <- g + ylim(10000, 35000) + xlab("Time") + 
#   ylab("Actual Hourly load (MW)")
# 
# #png(filename = "test_demand_actual.png", width = 3000, height = 2000, 
# #    res = 300)
# print(g)
# #dev.off()
# 
# # plot estimated hourly demand in test data set
# xx <- data.frame(x=df.test.season$time, y=exp(y.hat)) 
# g <- ggplot()
# g <- g + geom_line(data=xx, aes(x=x, y=y), 
#                    colour = "darkgray")
# g <- g + theme_bw() + scale_x_datetime(breaks=ticks.plot, date_labels = "%b")
# g <- g + ylim(10000, 35000) + xlab("Time") + 
#   ylab("Estimated Hourly load (MW)")
# 
# #png(filename = "test_demand_estimate.png", width = 3000, height = 2000, 
# #    res = 300)
# print(g)
# #dev.off()
# 
# 
# mm <- df.test.season %>% group_by(hour.of.day) %>% select(load) %>%
#   summarize(mean.value = mean(load, na.rm=TRUE))
# 
# df.xx <- data.frame(time=df.test.season$time, 
#                     hour.of.day=df.test.season$hour.of.day, 
#                     load=y.hat)
# mm.xx <- as.data.frame(df.xx %>% group_by(hour.of.day) %>% 
#                          select(load) %>%
#                          summarize(mean.value = mean(load, na.rm=TRUE)))
# 
# plot(mm$mean.value, type="l", col="blue")
# lines(mm.xx$mean.value, col="red")
#
# regression.model.2 <- function(data.set, 
#                                temp.breaks = c(-10, 0, 10, 20, 30),
#                                dew.point.inter = TRUE) {
#   # Fits regression model to complete data set
#   #
#   # Args:
#   #   data.set:  data frame with data set to be used in regression
#   #   temp.breaks:  
#   #   dew.point.inter:  
#   #
#   # Returns:
#   #   reg.model: regression model
#   
#   temp.components <- createTempComponents(data.set$temperature)
#   data.set <- cbind(data.set, temp.components)
#   rm(temp.components)
#   
#   # creates regression equation
#   equation.string <- "load ~ hour.of.day:type.day:season + factor(year)"
#   n.cp <- length(temp.breaks) + 1  # number of components
#   names.cp <- paste0("tc.", c(1:n.cp))
#   names.interaction <- paste0(names.cp,":dew.point")
#   equation.string <- paste0(equation.string, " + ",
#                             paste0(names.cp, collapse = " + "))
#   
#   if (dew.point.inter) {
#     # add interaction term between temperature components and dew point
#     equation.string <- paste0(equation.string, " + ",
#                               paste0(names.interaction, collapse = " + "))
#   }
#   
#   cat("\n--------------------------------------------\n")
#   cat(equation.string)
#   cat("\n--------------------------------------------\n")
#   
#   # fits model to training data set
#   reg.model <- lm(formula = as.formula(equation.string), 
#                   data = data.set)
#   
#   return(reg.model)
#   
# }