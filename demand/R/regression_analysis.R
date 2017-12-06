# ********** Regression Analysis *********************

read.from.file <- TRUE

setwd("~/GoogleDrive/CMU/RIPS/R")

source("auxiliary_functions.R")

# prepares data frame for regression
source("prepare_data_frame.R")
if(read.from.file) {
  df.final <- readDataFromFile()
} else {
  df.final <- readAllData()
}

# df.final <- df.final[complete.cases(df.final$temperature), ]

# loads function with regression analysis
source("regression4.R")

# df.final <- df.final[df.final$year >= 2006, ]
#m1 <- regression.model.1(df.final)
#m2 <- regression.model.1(df.final, dew.point.inter = FALSE)

lm.base <- regression.model.1(df.final, dew.point.inter = FALSE)

lm.dp.inter <- regression.model.2(df.final)

fix.eff.year <- m1$coefficients[grep("Intercept", names(m1$coefficients))] + 
  m1$coefficients[grep("year", names(m1$coefficients))]

fix.eff.year <- c(m1$coefficients[grep("Intercept", names(m1$coefficients))],
                  fix.eff.year)

fix.eff.year <- data.frame(year=c(c(1993:2003), c(2006:2014)), 
                                  fix.eff = fix.eff.year)
row.names(fix.eff.year) <- NULL

xx <- merge(fix.eff.year,
            select(df.final, c(1,11:26)),
            by="year")

# drop repeated rows
xx <- unique(xx)

xx$GDP.per.capita.US <- xx$GDP.US / xx$POP.US

# Fixed effects vs US GDP ----
g <- ggplot(data=xx)
g <- g + geom_point(aes(x=GDP.US, y=fix.eff)) + theme_bw()
g <- g + geom_text(aes(x=GDP.US, y=fix.eff, label=year), 
                   hjust=0, nudge_x = 1e4) 
coef.line <- coef(lm(fix.eff ~ GDP.US, data = xx))
g <- g + geom_abline(intercept = coef.line[1], slope = coef.line[2])
g
summary(lm(fix.eff ~ GDP.US, data = xx))

# Fixed effects vs US Population ----
g <- ggplot(data=xx)
g <- g + geom_point(aes(x=POP.US, y=fix.eff)) + theme_bw()
g <- g + geom_text(aes(x=POP.US, y=fix.eff, label=year), 
                   hjust=0, nudge_x = 1e4) 
coef.line <- coef(lm(fix.eff ~ POP.US, data = xx))
g <- g + geom_abline(intercept = coef.line[1], slope = coef.line[2])
g
summary(lm(fix.eff ~ POP.US, data = xx))

# Fixed effects vs US Population + GDP ----
summary(lm(fix.eff ~ GDP.US + POP.US, data = xx))

# Fixed effects vs GDP per capita ----
g <- ggplot(data=xx)
g <- g + geom_point(aes(x=GDP.per.capita.US, y=fix.eff)) + theme_bw()
g <- g + geom_text(aes(x=GDP.per.capita.US, y=fix.eff, label=year), 
                   hjust=0) 
coef.line <- coef(lm(fix.eff ~ GDP.per.capita.US, data = xx))
g <- g + geom_abline(intercept = coef.line[1], slope = coef.line[2])
g
summary(lm(fix.eff ~ GDP.per.capita.US, data = xx))

# Fixed effects vs TN Population + GDP ----
g <- ggplot(data=xx)
g <- g + geom_point(aes(x=GDP.Tennessee, y=fix.eff)) + theme_bw()
g <- g + geom_text(aes(x=GDP.Tennessee, y=fix.eff, label=year), 
                   hjust=0, nudge_x = 5e2) 
coef.line <- coef(lm(fix.eff ~ GDP.Tennessee, data = xx))
g <- g + geom_abline(intercept = coef.line[1], slope = coef.line[2])
g
summary(lm(fix.eff ~ GDP.Tennessee, data = xx))

# Linear regression including all GDPs and POPs of TVA states ----
formula.lm <- paste0("fix.eff ~ ", paste0(names(xx)[4:10], 
                                          collapse = " + "))
formula.lm <- paste0(formula.lm, " + ", paste0(names(xx)[12:18], 
                                               collapse = " + "))
xx.train <- xx[-nrow(xx), ]
lm1 <- lm(formula = formula.lm, data = xx.train)
summary(lm1)

y <- predict(lm1, newdata = xx[nrow(xx), ])
abs(y - xx$fix.eff[nrow(xx)])/xx$fix.eff[nrow(xx)]

# Comparison GDP growth per state

gdp.df <- data.frame(year=df.final$year, 
                     select(df.final, starts_with("GDP.")))
gdp.df <- unique(gdp.df)